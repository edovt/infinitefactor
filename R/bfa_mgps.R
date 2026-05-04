#' Bayesian Factor Analysis and Covariance Estimation with an MGPS prior
#'
#' @param X Data matrix
#' @param iter_warmup Number of warmup iterations to run per chain.
#' @param iter_sampling Number of post-warmup iterations to run per chain.
#' @param adapt (logical) Adapt number of factors across samples?
#' @param k_init If `adapt = FALSE`, the persistent number of factors. If
#'     `adapt = TRUE`, the initial number of factors in the sampler. If `NULL`
#'     the number of factors is selected automatically as `3*log(p)` where `p`
#'     is the number of columns of `X`.
#' @param prior_params Parameters of the MGPS prior.
#' @param eps If `adapt = TRUE`, tolerance for column sums of `Lambda` in the
#'     adaptation step.
#' @param verbose (logical) Show progress bar?
#'
#' @returns List with samples for each parameter
#'
#'
#' @export
#' @examples
#' NULL
bfa_mgps <- function(
  X,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt = TRUE,
  k_init = NULL,
  prior_params = list(
    a_sigma = NULL,
    b_sigma = NULL,
    nu = NULL,
    a1 = NULL,
    a2 = NULL,
    alpha0 = NULL,
    alpha1 = NULL
  ),
  eps = 1e-2,
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
      stop(
        "'k_init' must be less than min(n, p)."
      )
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
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0) {
    stop("'eps' must be a single positive number.")
  }

  a_sigma_check <- prior_params$a_sigma %||% 1
  b_sigma_check <- prior_params$b_sigma %||% 1
  nu_check <- prior_params$nu %||% 3
  a1_check <- prior_params$a1 %||% 2.1
  a2_check <- prior_params$a2 %||% 3.1
  if (a_sigma_check <= 0 || b_sigma_check <= 0) {
    stop("'a_sigma' and 'b_sigma' must be positive.")
  }
  if (nu_check <= 0) {
    stop("'nu' must be positive.")
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
  nu <- prior_params$nu %||% 3
  a1 <- prior_params$a1 %||% 2.1
  a2 <- prior_params$a2 %||% 3.1
  if (adapt) {
    alpha0 <- prior_params$alpha0 %||% -1
    alpha1 <- prior_params$alpha1 %||% -0.0005
  }

  # 2. Set-up storage of samples -------------------------------------------
  if (!adapt) {
    samples <- list(
      Lambda = array(NA, dim = c(p, k, iter_sampling)),
      Eta = array(NA, dim = c(n, k, iter_sampling)),
      sigma2_inv = matrix(NA, nrow = p, ncol = iter_sampling),
      Phi = array(NA, dim = c(p, k, iter_sampling)),
      delta = matrix(NA, nrow = k, ncol = iter_sampling),
      tau = matrix(NA, nrow = k, ncol = iter_sampling),
      Sigma_X = array(NA, dim = c(p, p, iter_sampling))
    )
  } else {
    samples <- list(
      k = numeric(iter_sampling),
      Sigma_X = array(NA, dim = c(p, p, iter_sampling))
    )
  }

  # 3. Initialize parameters -----------------------------------------------
  # Eta not needed since it's the first to get updated
  delta_init <- c(stats::rgamma(1, a1, 1), stats::rgamma(k - 1, a2, 1))
  mgps_params <- list(
    Lambda = matrix(stats::rnorm(p * k), nrow = p, ncol = k),
    sigma2_inv = stats::rgamma(p, shape = a_sigma, rate = b_sigma),
    Phi = matrix(stats::rgamma(p * k, nu / 2, nu / 2), nrow = p, ncol = k),
    delta = delta_init,
    tau = cumprod(delta_init)
  )

  # 4. Gibbs sampler -------------------------------------------------------
  total_iter <- iter_warmup + iter_sampling
  shape_sigma <- a_sigma + n / 2 # constant across iterations
  shape_phi <- (nu + 1) / 2 # constant across iterations
  p_bar <- utils::txtProgressBar(max = total_iter, style = 3)

  for (iter in 1:total_iter) {
    Eta <- latent_update_bfa(
      X,
      mgps_params$Lambda,
      mgps_params$sigma2_inv,
      n,
      k
    )
    mgps_params <- MGPS_update(
      mgps_params,
      X,
      Eta,
      n,
      p,
      k,
      shape_sigma,
      b_sigma,
      shape_phi,
      nu,
      a1,
      a2
    )

    # 4.2 If k is not fixed, update
    if (adapt) {
      p_adapt <- exp(alpha0 + alpha1 * iter)
      unif <- stats::runif(1)
      if (unif < p_adapt) {
        redundant <- (colSums(abs(mgps_params$Lambda) < eps) == p)
        n_redundant <- sum(redundant)
        if (n_redundant == 0) {
          k <- k + 1
          Eta <- cbind(Eta, stats::rnorm(n))
          mgps_params$Phi <- cbind(
            mgps_params$Phi,
            stats::rgamma(p, nu / 2, nu / 2)
          )
          mgps_params$delta[k] <- stats::rgamma(1, a2, 1)
          mgps_params$tau <- cumprod(mgps_params$delta)
          mgps_params$Lambda <- cbind(
            mgps_params$Lambda,
            stats::rnorm(
              p,
              0,
              1 / sqrt(mgps_params$Phi[, k] * mgps_params$tau[k])
            )
          )
        } else {
          # TODO: there is a bug here if n_redundant = k
          k <- max(k - n_redundant, 1)
          mgps_params$Lambda <- mgps_params$Lambda[, !redundant, drop = F]
          mgps_params$Phi <- mgps_params$Phi[, !redundant, drop = F]
          Eta <- Eta[, !redundant, drop = F]
          mgps_params$delta <- mgps_params$delta[!redundant]
          mgps_params$tau <- cumprod(mgps_params$delta)
        }
      }
    }

    # 4.3 Save samples
    if (iter > iter_warmup) {
      c_iter <- iter - iter_warmup
      samples$Sigma_X[,, c_iter] <- (mgps_params$Lambda %*%
        t(mgps_params$Lambda) +
        diag(1 / mgps_params$sigma2_inv)) *
        scale_mat
      if (!adapt) {
        samples$Lambda[,, c_iter] <- mgps_params$Lambda
        samples$Eta[,, c_iter] <- Eta
        samples$sigma2_inv[, c_iter] <- mgps_params$sigma2_inv
        samples$Phi[,, c_iter] <- mgps_params$Phi
        samples$delta[, c_iter] <- mgps_params$delta
        samples$tau[, c_iter] <- mgps_params$tau
      } else {
        samples$k[c_iter] <- k
      }
    }

    if (verbose) utils::setTxtProgressBar(p_bar, iter)
  }

  close(p_bar)
  samples
}
