#' Bayesian Factor Analysis and Covariance Estimation with a CUSP prior
#'
#' @param X Data matrix
#' @param iter_warmup Number of warmup iterations to run per chain
#' @param iter_sampling Number of post-warmup iterations to run per chain
#' @param adapt (logical) Adapt number of factors across samples?
#' @param k_init If `adapt=FALSE`, the persistent number of factors. If
#'     `adapt=TRUE`, the initial number of factors in the sampler. If `NULL`
#'     the number of factors is selected automatically as `3*log(p)`, where `p`
#'     is the number of columns of `X`
#' @param prior_params Parameters of the CUSP prior
#' @param verbose (logical) Show progress bar?
#'
#' @returns List with samples for each parameter
#'
#' @export
#' @examples
#' NULL
bfa_cusp <- function(
  X,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt = TRUE,
  k_init = NULL,
  prior_params = list(
    a_sigma = NULL,
    b_sigma = NULL,
    alpha = NULL,
    a_theta = NULL,
    b_theta = NULL,
    theta_inf = NULL
  ),
  verbose = TRUE
) {
  # 0. Input checks --------------------------------------------------------
  if (!is.matrix(X)) {
    X <- as.matrix(X)
    if (!is.numeric(X)) {
      stop("'X' must be a numeric matrix.")
    }
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop("'X' must not contain NA, NaN, or Inf values.")
  }
  if (nrow(X) <= 1 || ncol(X) < 2) {
    stop("'X' must have at least 2 rows and 2 columns.")
  }
  if (any(apply(X, 2, stats::sd) == 0)) {
    stop("'X' has constant columns; remove or impute them before fitting.")
  }

  if (!is.null(k_init)) {
    if (
      !is.numeric(k_init) ||
        length(k_init) != 1 ||
        k_init != as.integer(k_init) ||
        k_init < 1
    ) {
      stop("'k_init' must be a positive integer.")
    }
    if (k_init >= min(nrow(X), ncol(X))) {
      stop("'k_init' must be less than min(n, p).")
    }
  }

  if (
    !is.numeric(iter_warmup) ||
      length(iter_warmup) != 1 ||
      iter_warmup != as.integer(iter_warmup) ||
      iter_warmup < 1
  ) {
    stop("'iter_warmup' must be a positive integer.")
  }
  if (
    !is.numeric(iter_sampling) ||
      length(iter_sampling) != 1 ||
      iter_sampling != as.integer(iter_sampling) ||
      iter_sampling < 1
  ) {
    stop("'iter_sampling' must be a positive integer.")
  }

  if (!is.logical(adapt) || length(adapt) != 1) {
    stop("'adapt' must be a single logical value.")
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be a single logical value.")
  }

  a_sigma_check <- prior_params$a_sigma %||% 1
  b_sigma_check <- prior_params$b_sigma %||% 1
  alpha_check <- prior_params$alpha %||% 5
  a_theta_check <- prior_params$a_theta %||% 2
  b_theta_check <- prior_params$b_theta %||% 2
  theta_inf_check <- prior_params$theta_inf %||% 0.05
  if (a_sigma_check <= 0 || b_sigma_check <= 0) {
    stop("'a_sigma' and 'b_sigma' must be positive.")
  }
  if (alpha_check <= 0) {
    stop("'alpha' must be positive.")
  }
  if (a_theta_check <= 0 || b_theta_check <= 0) {
    stop("'a_theta' and 'b_theta' must be positive.")
  }
  if (theta_inf_check <= 0) {
    stop("'theta_inf' must be positive.")
  }
  if (adapt) {
    alpha1_check <- prior_params$alpha1 %||% -0.0005
    if (alpha1_check >= 0) {
      stop(
        "'alpha1' must be negative so that the adaptation probability decays over iterations."
      )
    }
  }

  # 1. Extract dimensions and hyperparameters ------------------------------
  n <- nrow(X)
  p <- ncol(X)
  sd_X <- apply(X, 2, stats::sd)
  scale_mat <- outer(sd_X, sd_X)
  X <- scale(X)

  k <- k_init %||% as.integer(3 * log(p))
  a_sigma <- prior_params$a_sigma %||% 1
  b_sigma <- prior_params$b_sigma %||% 1
  alpha <- prior_params$alpha %||% 5
  a_theta <- prior_params$a_theta %||% 2
  b_theta <- prior_params$b_theta %||% 2
  theta_inf <- prior_params$theta_inf %||% 0.05
  if (adapt) {
    alpha0 <- prior_params$alpha0 %||% -1
    alpha1 <- prior_params$alpha1 %||% -0.0005
  }

  # 2. Set-up storage of samples -------------------------------------------
  if (!adapt) {
    samples <- list(
      Lambda = array(dim = c(p, k, iter_sampling)),
      Eta = array(dim = c(n, k, iter_sampling)),
      sigma2_inv = matrix(nrow = p, ncol = iter_sampling),
      theta = matrix(nrow = k, ncol = iter_sampling),
      Sigma_X = array(dim = c(p, p, iter_sampling))
    )
  } else {
    samples <- list(
      k = numeric(iter_sampling),
      Sigma_X = array(dim = c(p, p, iter_sampling))
    )
  }

  # 3. Initialize parameters -----------------------------------------------
  # The update will add z necessary for adaptation
  cusp_params <- list(
    Lambda = matrix(stats::rnorm(p * k), p, k),
    sigma2_inv = stats::rgamma(p, a_sigma, b_sigma),
    theta = stats::rgamma(k, a_theta, b_theta),
    omega = stick_break(c(stats::rbeta(k - 1, 1, alpha), 1))
  )

  # 4. Gibbs sampler -------------------------------------------------------
  total_iter <- iter_warmup + iter_sampling
  shape_sigma <- a_sigma + n / 2
  shape_theta <- a_theta + p / 2
  p_bar <- utils::txtProgressBar(max = total_iter, style = 3)
  for (iter in 1:total_iter) {
    Eta <- latent_update_bfa(
      X,
      cusp_params$Lambda,
      cusp_params$sigma2_inv,
      n,
      k
    )
    cusp_params <- CUSP_update(
      cusp_params,
      X,
      Eta,
      n,
      p,
      k,
      shape_sigma,
      b_sigma,
      alpha,
      shape_theta,
      a_theta,
      b_theta,
      theta_inf
    )

    if (adapt) {
      p_adapt <- exp(alpha0 + alpha1 * iter)
      if (stats::runif(1) < p_adapt) {
        active <- cusp_params$z > seq_len(k)
        k_star <- sum(active)
        # k - 1 since last element is always from the spike
        if (k_star < k - 1) {
          k <- k_star + 1
          cusp_params$Lambda <- cbind(
            cusp_params$Lambda[, active, drop = FALSE],
            stats::rnorm(p, 0, sqrt(theta_inf))
          )
          cusp_params$theta <- c(cusp_params$theta[active], theta_inf)
          Eta <- cbind(Eta[, active, drop = FALSE], stats::rnorm(n))
          v_new <- c(stats::rbeta(k_star, 1, alpha), 1)
          cusp_params$omega <- stick_break(v_new)
        } else {
          k <- k + 1
          Eta <- cbind(Eta, stats::rnorm(n))
          new_v <- stats::rbeta(1, 1, alpha)
          cusp_params$omega <- c(
            cusp_params$omega[seq_len(k - 2)],
            cusp_params$omega[k - 1] * new_v
          )
          cusp_params$omega <- c(cusp_params$omega, 1 - sum(cusp_params$omega))
          cusp_params$theta <- c(cusp_params$theta[seq_len(k - 1)], theta_inf)
          cusp_params$Lambda <- cbind(
            cusp_params$Lambda,
            stats::rnorm(p, 0, sqrt(theta_inf))
          )
        }
      }
    }

    if (iter > iter_warmup) {
      c_iter <- iter - iter_warmup
      samples$Sigma_X[,, c_iter] <- (cusp_params$Lambda %*%
        t(cusp_params$Lambda) +
        diag(1 / cusp_params$sigma2_inv)) *
        scale_mat
      if (!adapt) {
        samples$Lambda[,, c_iter] <- cusp_params$Lambda
        samples$Eta[,, c_iter] <- Eta
        samples$sigma2_inv[, c_iter] <- cusp_params$sigma2_inv
        samples$theta[, c_iter] <- cusp_params$theta
      } else {
        samples$k[c_iter] <- k
      }
    }

    if (verbose) utils::setTxtProgressBar(p_bar, iter)
  }

  close(p_bar)
  samples
}
