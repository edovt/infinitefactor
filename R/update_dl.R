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
  Phi <- dl_params$Phi
  Psi <- dl_params$Psi
  tau <- dl_params$tau

  # 4.1 Update Lambda row-wise
  Eta_cross <- crossprod(Eta)
  EtaTX <- crossprod(Eta, X)
  for (j in 1:p) {
    D_j <- diag(1 / (Psi[j, ] * Phi[j, ] * tau^2), nrow = k)
    V_j <- solve(D_j + sigma2_inv[j] * Eta_cross)
    mu_j <- sigma2_inv[j] * (V_j %*% EtaTX[, j])
    Lambda[j, ] <- mvtnorm::rmvnorm(1, mu_j, V_j)
  }

  # 4.2 Update sigma_j
  rate_sigma <- b_sigma + .5 * colSums((X - Eta %*% t(Lambda))^2)
  sigma2_inv <- stats::rgamma(p, shape_sigma, rate_sigma)

  # 4.3 Update Psi
  for (j in 1:p) {
    for (h in 1:k) {
      Psi[j, h] <- GIGrvg::rgig(1, -1 / 2, chi = 1, psi = tau[j] * Phi[j, h])
    }
  }

  # 4.4 Update tau
  chi_tau <- 2 * rowSums(abs(Lambda) / Phi)
  for (j in 1:p) {
    tau[j] <- GIGrvg::rgig(1, 1 - k, chi = chi_tau[j], psi = 1)
  }

  # 4.5 Update Phi
  T_m <- matrix(NA, nrow = p, ncol = k)
  for (j in 1:p) {
    for (h in 1:k) {
      T_m[j, h] <- GIGrvg::rgig(1, a - 1, chi = 2 * abs(Lambda[j, h]), psi = 1)
    }
  }
  Psi <- T_m / rowSums(T_m)

  list(
    Lambda = Lambda,
    sigma2_inv = sigma2_inv,
    Phi = Phi,
    Psi = Psi,
    tau = tau
  )
}
