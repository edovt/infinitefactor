DL_update <- function(
  dl_params,
  X,
  Eta,
  n,
  p,
  k,
  shape_sigma,
  b_sigma,
  a
) {
  Lambda <- dl_params$Lambda
  sigma2_inv <- dl_params$sigma2_inv
  Delta <- dl_params$Delta
  Psi <- dl_params$Psi

  # 4.1 Update Lambda row-wise
  Eta_cross <- crossprod(Eta)
  EtaTX <- crossprod(Eta, X)
  for (j in 1:p) {
    D_j_inv <- diag(1 / (Psi[j, ] * Delta[j, ]^2), nrow = k)
    V_j <- solve(D_j_inv + sigma2_inv[j] * Eta_cross)
    mu_j <- sigma2_inv[j] * (V_j %*% EtaTX[, j])
    Lambda[j, ] <- mvnfast::rmvn(1, mu_j, V_j)
  }

  # 4.2 Update sigma_j
  rate_sigma <- b_sigma + .5 * colSums((X - Eta %*% t(Lambda))^2)
  sigma2_inv <- stats::rgamma(p, shape_sigma, rate_sigma)

  # 4.3 Update Delta
  for (j in 1:p) {
    for (h in 1:k) {
      Delta[j, h] <- GIGrvg::rgig(
        1,
        a - 1,
        chi = 2 * abs(Lambda[j, h]),
        psi = 1
      )
    }
  }

  # 4.4 Update Psi
  ## For the inverse Gaussian, order of parameters change
  for (j in 1:p) {
    for (h in 1:k) {
      Psi[j, h] <- GIGrvg::rgig(
        1,
        1 / 2,
        chi = (Lambda[j, h] / Delta[j, h])^2,
        psi = 1
      )
    }
  }

  list(
    Lambda = Lambda,
    sigma2_inv = sigma2_inv,
    Delta = Delta,
    Psi = Psi
  )
}
