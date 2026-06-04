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
  Delta_dl <- dl_params$Delta_dl
  Psi <- dl_params$Psi

  # 4.1 Update Delta_dl
  for (j in 1:p) {
    for (h in 1:k) {
      Delta_dl[j, h] <- GIGrvg::rgig(
        1,
        a - 1,
        chi = 2 * abs(Lambda[j, h]),
        psi = 1
      )
    }
  }

  # 4.2 Update Psi
  ## For the inverse Gaussian, order of parameters change
  for (j in 1:p) {
    for (h in 1:k) {
      Psi[j, h] <- GIGrvg::rgig(
        1,
        1 / 2,
        chi = (Lambda[j, h] / Delta_dl[j, h])^2,
        psi = 1
      )
    }
  }

  # 4.3 Update Lambda row-wise
  Lambda <- update_Lambda_dl(Eta, X, sigma2_inv, Psi, Delta_dl)

  # 4.4 Update sigma_j
  rate_sigma <- b_sigma + .5 * colSums((X - Eta %*% t(Lambda))^2)
  sigma2_inv <- stats::rgamma(p, shape_sigma, rate_sigma)

  list(
    Lambda = Lambda,
    sigma2_inv = sigma2_inv,
    Delta_dl = Delta_dl,
    Psi = Psi
  )
}
