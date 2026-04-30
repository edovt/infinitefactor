bfa_dl <- function(
  X,
  way = "rows",
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
  # 0. Input checks

  # 1. Extract dimensions and hyperparameters
  n <- nrow(X)
  p <- ncol(X)
  sd_X <- apply(X, 2, stats::sd)
  scale_mat <- outer(sd_X, sd_X)
  X <- scale(X)

  k <- k %||% as.integer(3 * log(p))
  a_sigma <- prior_params$a_sigma %||% 1
  b_sigma <- prior_params$b_sigma %||% 1
  a <- prior_params$a %||% 1

  # 2. Set-up storage of samples
  samples <- list(
    Lambda = array(NA, dim = c(p, k, iter_sampling)),
    Eta = array(NA, dim = c(n, k, iter_sampling)),
    sigma2_inv = matrix(NA, nrow = p, ncol = iter_sampling),
    Phi = array(NA, dim = c(p, k, iter_sampling)),
    tau = matrix(NA, nrow = p, ncol = iter_sampling),
    Omega = array(NA, dim = c(p, p, iter_sampling))
  )

  # 3. Initialize parameters
  dir_gamma_matrix <- matrix(stats::rgamma(p * k, a), nrow = p, ncol = k)
  dl_params <- list(
    Lambda = matrix(stats::rnorm(p * k), nrow = p, ncol = k),
    sigma2_inv = stats::rgamma(p, shape = a_sigma, rate = b_sigma),
    Phi = dir_gamma_matrix / rowSums(dir_gamma_matrix),
    Psi = matrix(rexp(p * k, 1 / 2), nrow = p, ncol = k),
    tau = stats::rgamma(p, n * a, 1 / 2),
  )

  # 4. Gibbs sampler
  total_iter <- iter_warmup + iter_sampling
  shape_sigma <- a_sigma + n / 2
  p_bar <- utils::txtProgressBar(max = total_iter, style = 3)
  for (iter in 1:total_iter) {
    # 4.1 All updates
    Eta <- latent_update_bfa(X, dl_params$Lambda, dl_params$sigma2_inv, n, k)
    dl_params <- DL_update(dl_params, X, Eta, n, p, k, shape_sigma, b_sigma)

    # 4.2 Save samples
    if (iter > iter_warmup) {
      c_iter <- iter - iter_warmup
      samples$Lambda[,, c_iter] <- dl_params$Lambda
      samples$Eta[,, c_iter] <- Eta
      samples$sigma2_inv[, c_iter] <- dl_params$sigma2_inv
      samples$Phi[,, c_iter] <- dl_params$Phi
      samples$tau[, c_iter] <- dl_params$tau
      samples$Omega[,, c_iter] <- (dl_params$Lambda %*%
        t(dl_params$Lambda) +
        diag(1 / dl_params$sigma2_inv)) *
        scale_mat
    }

    if (verbose) utils::setTxtProgressBar(p_bar, iter)
  }

  close(p_bar)
  samples
}
