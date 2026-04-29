MGPS_update <- function(
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
) {
  list2env(mgps_params, envir = environment()) # Phi, tau, sigma2_inv, Lambda

  # 4.1 Update Lambda row-wise
  Eta_cross <- crossprod(Eta)
  EtaTX <- crossprod(Eta, X)
  for (j in 1:p) {
    D_j <- diag(Phi[j, ] * tau, nrow = k)
    V_j <- solve(D_j + sigma2_inv[j] * Eta_cross)
    mu_j <- sigma2_inv[j] * (V_j %*% EtaTX[, j])
    Lambda[j, ] <- mvtnorm::rmvnorm(1, mu_j, V_j)
  }

  # 4.2 Update sigma_j
  rate_sigma <- b_sigma + .5 * colSums((X - Eta %*% t(Lambda))^2)
  sigma2_inv <- stats::rgamma(p, shape_sigma, rate_sigma)

  # 4.3 Update Phi
  rate_phi <- (nu + t(t(Lambda^2) * tau)) / 2
  Phi <- matrix(
    stats::rgamma(p * k, shape_phi, rate_phi),
    nrow = p,
    ncol = k
  )

  # 4.4 Update delta and tau
  sum_phi_lambda2 <- colSums(Phi * Lambda^2)
  shape_d1 <- a1 + (p * k) / 2
  rate_d1 <- 1 + 0.5 * sum((tau / delta[1]) * sum_phi_lambda2)
  delta[1] <- stats::rgamma(1, shape_d1, rate_d1)
  tau <- cumprod(delta)

  if (k >= 2) {
    for (h in 2:k) {
      shape_dh <- a2 + p * (k - h + 1) / 2
      tau_l_h <- tau[h:k] / delta[h]
      rate_dh <- 1 + 0.5 * sum(tau_l_h * sum_phi_lambda2[h:k])
      delta[h] <- stats::rgamma(1, shape_dh, rate_dh)
      tau <- cumprod(delta)
    }
  }

  list(
    Lambda = Lambda,
    sigma2_inv = sigma2_inv,
    Phi = Phi,
    delta = delta,
    tau = tau
  )
}
