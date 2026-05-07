#' Bayesian Latent Factor Regression with a CUSP prior
#'
#' @param y Response vector
#' @param X Predictors matrix
#' @param iter_warmup Number of warmup iterations to run per chain.
#' @param iter_sampling Number of post-warmup iterations to run per chain.
#' @param k Number of factors
#' @param induced (logical) return induced effects of original predictors?
#' @param interactions (logical) include interactions between latent factors?
#' @param prior_params Parameters of the CUSP prior
#' @param mala_eps If `interactions=TRUE` step-size of the MALA step for
#'     latent factors. If `NULL`, uses `k^(-1/3)` as default.
#' @param adapt_mala_eps (logical) Should `mala_eps` be adaptive?
#' @param window_mala IF `adapt_mala_eps=TRUE`, number of iterations during
#'     warmup that are considered for the adaptation of `mala_eps`
#' @param verbose (logical) Show progress bar?
#'
#' @returns List with samples for each parameter
#'
#' @export
#' @examples
#' NULL
blfr_cusp <- function(
  y,
  X,
  iter_warmup = 1000,
  iter_sampling = 1000,
  k = NULL,
  induced = TRUE,
  interactions = FALSE,
  prior_params = list(
    nu_y = NULL,
    beta_var = NULL,
    Omega_var = NULL,
    a_sigma = NULL,
    b_sigma = NULL,
    alpha = NULL,
    a_theta = NULL,
    b_theta = NULL,
    theta_inf = NULL
  ),
  mala_eps = NULL,
  adapt_mala_eps = TRUE,
  window_mala = floor(iter_warmup / 10),
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

  if (!is.numeric(y) || !is.null(dim(y))) {
    stop("'y' must be a numeric vector.")
  }
  if (length(y) != nrow(X)) {
    stop("'y' must have length equal to nrow(X).")
  }
  if (anyNA(y) || any(!is.finite(y))) {
    stop("'y' must not contain NA, NaN, or Inf values.")
  }

  if (!is.null(k)) {
    if (
      !is.numeric(k) ||
        length(k) != 1 ||
        k != as.integer(k) ||
        k < 1
    ) {
      stop("'k' must be a positive integer.")
    }
    if (k >= min(nrow(X), ncol(X))) {
      stop("'k' must be less than min(n, p).")
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

  if (!is.logical(induced) || length(induced) != 1) {
    stop("'induced' must be a single logical value.")
  }
  if (!is.logical(interactions) || length(interactions) != 1) {
    stop("'interactions' must be a single logical value.")
  }
  if (!is.logical(adapt_mala_eps) || length(adapt_mala_eps) != 1) {
    stop("'adapt_mala_eps' must be a single logical value.")
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be a single logical value.")
  }

  if (!is.null(mala_eps)) {
    if (!is.numeric(mala_eps) || length(mala_eps) != 1 || mala_eps <= 0) {
      stop("'mala_eps' must be a single positive number.")
    }
  }
  if (interactions && adapt_mala_eps) {
    if (
      !is.numeric(window_mala) ||
        length(window_mala) != 1 ||
        window_mala != as.integer(window_mala) ||
        window_mala < 1
    ) {
      stop("'window_mala' must be a positive integer.")
    }
  }

  nu_y_check <- prior_params$nu_y %||% 1
  beta_var_check <- prior_params$beta_var %||% 100
  Omega_var_check <- prior_params$Omega_var %||% 100
  a_sigma_check <- prior_params$a_sigma %||% 1
  b_sigma_check <- prior_params$b_sigma %||% 1
  alpha_check <- prior_params$alpha %||% 5
  a_theta_check <- prior_params$a_theta %||% 2
  b_theta_check <- prior_params$b_theta %||% 2
  theta_inf_check <- prior_params$theta_inf %||% 0.05
  if (nu_y_check <= 0) {
    stop("'nu_y' must be positive.")
  }
  if (beta_var_check <= 0) {
    stop("'beta_var' must be positive.")
  }
  if (interactions && Omega_var_check <= 0) {
    stop("'Omega_var' must be positive.")
  }
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

  # 1. Extract dimensions and hyperparameters ------------------------------
  n <- nrow(X)
  p <- ncol(X)
  sd_X <- apply(X, 2, stats::sd)
  scale_mat <- outer(sd_X, sd_X)
  X <- scale(X)

  k <- k %||% as.integer(3 * log(p))
  mala_eps <- mala_eps %||% k^(-1 / 3)
  nu_y <- prior_params$nu_y %||% 1
  beta_var <- prior_params$beta_var %||% 100
  Omega_var <- prior_params$Omega_var %||% 100
  a_sigma <- prior_params$a_sigma %||% 1
  b_sigma <- prior_params$b_sigma %||% 1
  alpha <- prior_params$alpha %||% 5
  a_theta <- prior_params$a_theta %||% 2
  b_theta <- prior_params$b_theta %||% 2
  theta_inf <- prior_params$theta_inf %||% 0.05

  # 2. Set-up storage of samples -------------------------------------------
  samples <- list(
    sigma2_y = numeric(iter_sampling),
    beta = matrix(NA, k, iter_sampling),
    Lambda = array(NA, c(p, k, iter_sampling)),
    Eta = array(NA, c(n, k, iter_sampling)),
    sigma2_inv = matrix(NA, p, iter_sampling),
    theta = matrix(NA, k, iter_sampling)
  )
  if (induced) {
    samples$beta_X <- matrix(NA, p, iter_sampling)
  }
  if (interactions) {
    samples$Omega <- array(NA, c(k, k, iter_sampling))
  }
  if (induced && interactions) {
    samples$Omega_X <- array(NA, c(p, p, iter_sampling))
  }

  # 3. Initialize parameters -----------------------------------------------
  cusp_params <- list(
    Lambda = matrix(stats::rnorm(p * k), p, k),
    sigma2_inv = stats::rgamma(p, a_sigma, b_sigma),
    theta = stats::rgamma(k, a_theta, b_theta),
    omega = stick_break(c(stats::rbeta(k - 1, 1, alpha), 1))
  )
  Eta <- matrix(stats::rnorm(n * k), n, k)
  reg_params <- list(
    beta = stats::rnorm(k, sd = sqrt(beta_var)),
    sigma2_y = 1 / stats::rgamma(1, nu_y / 2, nu_y / 2)
  )
  if (interactions) {
    Omega_upper <- matrix(0, k, k)
    Omega_upper[upper.tri(Omega_upper, diag = TRUE)] <- stats::rnorm(
      k * (k + 1) / 2,
      sd = sqrt(Omega_var)
    )
    reg_params$Omega <- (Omega_upper + t(Omega_upper)) / 2
  }

  # 4. Gibbs sampler -------------------------------------------------------
  total_iter <- iter_warmup + iter_sampling
  p_bar <- utils::txtProgressBar(max = total_iter, style = 3)
  n_acceptances <- numeric(n)

  # Constants across iterations
  shape_sigmay <- (nu_y + n) / 2
  shape_sigma <- a_sigma + n / 2
  shape_theta <- a_theta + p / 2
  norm_covariance <- diag(theta_inf, p)
  t_covariance <- diag(b_theta / a_theta, p)

  for (iter in 1:total_iter) {
    # 4.1 Update parameters
    latent_update <- latent_update_blfr(
      X,
      y,
      Eta,
      reg_params,
      cusp_params$Lambda,
      cusp_params$sigma2_inv,
      n,
      k,
      mala_eps,
      n_acceptances
    )
    Eta <- latent_update$Eta
    n_acceptances <- latent_update$n_acceptances

    reg_params <- update_reg_params(
      reg_params,
      y,
      Eta,
      beta_var,
      Omega_var,
      shape_sigmay,
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
      theta_inf,
      norm_covariance,
      t_covariance
    )

    # 4.2 Adapt stepsize if wanted during warmup and after window
    if (adapt_mala_eps && (iter %% window_mala == 0) && iter <= iter_warmup) {
      acceptance_mean <- mean(n_acceptances) / window_mala
      mala_eps <- exp(log(mala_eps) + acceptance_mean - 0.574)
      n_acceptances <- numeric(n)
    }

    # 4.3 Save samples
    if (iter > iter_warmup) {
      c_iter <- iter - iter_warmup
      samples$sigma2_y[c_iter] <- reg_params$sigma2_y
      samples$beta[, c_iter] <- reg_params$beta
      samples$Lambda[,, c_iter] <- cusp_params$Lambda
      samples$Eta[,, c_iter] <- Eta
      samples$sigma2_inv[, c_iter] <- cusp_params$sigma2_inv
      samples$theta[, c_iter] <- cusp_params$theta

      if (interactions) {
        samples$Omega[,, c_iter] <- reg_params$Omega
      }
      if (induced) {
        L <- cusp_params$Lambda
        V <- solve(t(L) %*% diag(cusp_params$sigma2_inv) %*% L + diag(k))
        A <- V %*% t(L) %*% diag(cusp_params$sigma2_inv)
        samples$beta_X[, c_iter] <- t(A) %*% reg_params$beta

        if (interactions) {
          samples$Omega_X[,, c_iter] <- t(A) %*% reg_params$Omega %*% A
        }
      }
    }

    if (verbose) utils::setTxtProgressBar(p_bar, iter)
  }

  close(p_bar)
  samples
}
