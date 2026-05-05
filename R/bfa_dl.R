#' Bayesian Factor Analysis and Covariance Estimation with a row-wise DL prior
#'
#' @param X Data matrix
#' @param iter_warmup Number of warmup iterations to run per chain.
#' @param iter_sampling Number of post-warmup iterations to run per chain.
#' @param k Number of factors (fixed)
#' @param prior_params Parameters of the DL prior.
#' @param verbose (logical) Show progress bar?
#'
#' @returns
#'
#' @export
#' @examples
#' NULL
bfa_dl <- function(
  X,
  iter_warmup = 1000,
  iter_sampling = 1000L,
  k = NULL,
  prior_params = list(
    a_sigma = NULL,
    b_sigma = NULL,
    a = NULL
  ),
  verbose = TRUE
) {
  # 0. Input checks ---------------------------------------------------------

  # 1. Extract dimensions and hyperparameters -------------------------------
  n <- nrow(X)
  p <- ncol(X)
  sd_X <- apply(X, 2, stats::sd)
  scale_mat <- outer(sd_X, sd_X)
  X <- scale(X)

  k <- k %||% as.integer(3 * log(p))
  a_sigma <- prior_params$a_sigma %||% 1
  b_sigma <- prior_params$b_sigma %||% 1
  a <- prior_params$a %||% 1 / k

  # 2. Set-up storage of samples --------------------------------------------
  samples <- list(
    Lambda = array(NA, dim = c(p, k, iter_sampling)),
    Eta = array(NA, dim = c(n, k, iter_sampling)),
    sigma2_inv = matrix(NA, p, iter_sampling),
    Delta = array(NA, dim = c(p, k, iter_sampling)),
    Psi = array(NA, dim = c(p, k, iter_sampling)),
    Sigma_X = array(NA, dim = c(p, p, iter_sampling))
  )

  # 3. Initialize parameters ------------------------------------------------
  # Eta not needed since it's the first to get updated
  dl_params <- list(
    Lambda = matrix(stats::rnorm(p * k), p, k),
    sigma2_inv = stats::rgamma(p, a_sigma, b_sigma),
    Delta = matrix(stats::rgamma(p * k, a, 1 / 2), p, k),
    Psi = matrix(stats::rexp(p * k, 1 / 2))
  )

  # 4. Gibbs sampler --------------------------------------------------------
  total_iter <- iter_warmup + iter_sampling
  p_bar <- utils::txtProgressBar(max = total_iter, style = 3)

  # Constants across iterations
  shape_sigma <- a_sigma + n / 2

  for (iter in 1:total_iter) {
    # 4.1 All updates
    Eta <- latent_update_bfa(X, dl_params$Lambda, dl_params$sigma2_inv, n, k)
    dl_params <- DL_update(dl_params, X, Eta, n, p, k, shape_sigma, b_sigma, a)

    # 4.2 Save samples
    if (iter > iter_warmup) {
      c_iter <- iter - iter_warmup
      samples$Lambda[,, c_iter] <- dl_params$Lambda
      samples$Eta[,, c_iter] <- Eta
      samples$sigma2_inv[, c_iter] <- dl_params$sigma2_inv
      samples$Delta[,, c_iter] <- dl_params$Delta
      samples$Psi[,, c_iter] <- dl_params$Psi
      samples$Sigma_X[,, c_iter] <- (dl_params$Lambda %*%
        t(dl_params$Lambda) +
        diag(1 / dl_params$sigma2_inv)) *
        scale_mat
    }

    if (verbose) utils::setTxtProgressBar(p_bar, iter)
  }

  close(p_bar)
  samples
}
