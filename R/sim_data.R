sim_data_bfa <- function(n, p, k) {
  # Number of non-zeros per column, linearly decreasing from 2k down to k+1
  n_nonzero <- round(seq(2 * k, k + 1, length.out = k))

  Lambda <- matrix(0, nrow = p, ncol = k)
  for (h in seq_len(k)) {
    rows <- sample(p, n_nonzero[h])
    Lambda[rows, h] <- stats::rnorm(n_nonzero[h], mean = 0, sd = 3)
  }

  sigma2 <- 1 / stats::rgamma(p, shape = 1, rate = 0.25)
  Sigma_X_true <- tcrossprod(Lambda) + diag(sigma2)
  X <- mvnfast::rmvn(n, rep(0, p), Sigma_X_true)

  list(
    X = X,
    Lambda = Lambda,
    sigma2 = sigma2,
    Sigma_X_true = Sigma_X_true,
    k_true = k
  )
}

sim_data_blfr <- function() {}
