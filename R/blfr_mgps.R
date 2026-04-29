#' Bayesian Latent Factor Regression with an MGPS prior
#'
#' @param y Response vector
#' @param X Predictors matrix
#' @param iter_warmup Number of warmup iterations to run per chain.
#' @param iter_sampling Number of post-warmup iterations to run per chain.
#' @param k Number of factors
#' @param induced (logical) return induced effects of original predictors?
#' @param interactions (logical) include interactions between latent factors?
#' @param prior_params Parameters of the MGPS prior
#' @param mala_eps If `interactions=TRUE`, step-size of the MALA step for
#'     latent factors
#' @param verbose (logical) Show progress bar?
#'
#' @returns List with samples for each parameter
#'
#' @export
#' @examples
#' NULL
blfr_mgps <- function(
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
    nu = NULL,
    a1 = NULL,
    a2 = NULL
  ),
  mala_eps = 10,
  verbose = TRUE
) {
  # 0. Input checks --------------------------------------------------------

  # 1. Extract dimensions and hyperparameters ------------------------------
  n <- nrow(X)
  p <- ncol(X)
  sd_X <- apply(X, 2, stats::sd)
  scale_mat <- outer(sd_X, sd_X)
  X <- scale(X)

  k <- k %||% as.integer(3 * log(p))
  nu_y <- prior_params$nu_y %||% 1
  beta_var <- prior_params$beta_var %||% 100
  Omega_var <- prior_params$Omega_var %||% 100
  a_sigma <- prior_params$a_sigma %||% 1
  b_sigma <- prior_params$b_sigma %||% 1
  nu <- prior_params$nu %||% 3
  a1 <- prior_params$a1 %||% 2.1
  a2 <- prior_params$a2 %||% 3.1

  # 2. Set-up storage of samples -------------------------------------------
  samples <- list(
    sigma2_y = numeric(iter_sampling),
    beta = matrix(NA, nrow = k, ncol = iter_sampling),
    Lambda = array(NA, dim = c(p, k, iter_sampling)),
    Eta = array(NA, dim = c(n, k, iter_sampling)),
    sigma2_inv = matrix(NA, nrow = p, ncol = iter_sampling),
    Phi = array(NA, dim = c(p, k, iter_sampling)),
    delta = matrix(NA, nrow = k, ncol = iter_sampling),
    tau = matrix(NA, nrow = k, ncol = iter_sampling)
  )
  if (induced) {
    samples$beta_X = matrix(NA, nrow = p, ncol = iter_sampling)
  }
  if (interactions) {
    samples$Omega = array(NA, dim = c(k, k, iter_sampling))
  }
  if (induced && interactions) {
    samples$Omega_X = array(NA, dim = c(p, p, iter_sampling))
  }

  # 3. Initialize parameters -----------------------------------------------
  delta_init <- c(stats::rgamma(1, a1, 1), stats::rgamma(k - 1, a2, 1))
  mgps_params <- list(
    Lambda = matrix(stats::rnorm(p * k), nrow = p, ncol = k),
    sigma2_inv = stats::rgamma(p, shape = a_sigma, rate = b_sigma),
    Phi = matrix(stats::rgamma(p * k, nu / 2, nu / 2), nrow = p, ncol = k),
    delta = delta_init,
    tau = cumprod(delta_init)
  )
  Eta <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  reg_params <- list(
    beta = stats::rnorm(k, sd = sqrt(beta_var)),
    sigma2_y = 1 / stats::rgamma(1, nu_y / 2, nu_y / 2)
  )
  if (interactions) {
    Omega_upper <- matrix(
      stats::rnorm(k * (k + 1) / 2, sd = sqrt(Omega_var)),
      nrow = k,
      ncol = k
    )
    reg_params$Omega = (Omega_upper + t(Omega_upper)) / 2
  }

  # 4. Gibbs sampler -------------------------------------------------------
  total_iter <- iter_warmup + iter_sampling
  p_bar <- utils::txtProgressBar(max = total_iter, style = 3)

  # constants across iterations
  shape_sigmay <- (nu_y + n) / 2
  shape_sigma <- a_sigma + n / 2
  shape_phi <- (nu + 1) / 2

  for (iter in 1:total_iter) {
    Eta <- latent_update_blfr(
      X,
      y,
      Eta,
      reg_params,
      mgps_params$Lambda,
      mgps_params$sigma2_inv,
      n,
      k,
      mala_eps
    )

    reg_params <- update_reg_params(
      reg_params,
      y,
      Eta,
      shape_sigmay,
      beta_var,
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

    # 4.6 Save samples
    if (iter > iter_warmup) {
      c_iter <- iter - iter_warmup
      samples$Lambda[,, c_iter] <- mgps_params$Lambda
      samples$Eta[,, c_iter] <- Eta
      samples$sigma2_inv[, c_iter] <- mgps_params$sigma2_inv
      samples$Phi[,, c_iter] <- mgps_params$Phi
      samples$delta[, c_iter] <- mgps_params$delta
      samples$tau[, c_iter] <- mgps_params$tau

      if (interactions) {
        samples$Omega[,, c_iter] <- reg_params$Omega
      }
      if (induced) {
        L <- mgps_params$Lambda
        V <- solve(t(L) %*% diag(mgps_params$sigma2_inv) %*% L + diag(k))
        A <- V %*% t(L) %*% diag(mgps_params$sigma2_inv)
        samples$beta_X[, c_iter] <- t(A) %*% reg_params$beta
      }
      if (induced && interactions) {
        samples$Omega_X[,, c_iter] <- t(A) %*% reg_params$Omega %*% A
      }
    }

    if (verbose) utils::setTxtProgressBar(p_bar, iter)
  }

  close(p_bar)
  samples
}
