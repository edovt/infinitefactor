CUSP_update <- function(
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
) {
  Lambda <- cusp_params$Lambda
  sigma2_inv <- cusp_params$sigma2_inv
  theta <- cusp_params$theta
  omega <- cusp_params$omega

  # 1. Update Lambda row-wise
  Eta_cross <- crossprod(Eta)
  EtaTX <- crossprod(Eta, X)
  for (j in 1:p) {
    V_j <- solve(diag(theta, nrow = k) + sigma2_inv[j] * Eta_cross)
    mu_j <- sigma2_inv[j] * (V_j %*% EtaTX[, j])
    Lambda[j, ] <- mvnfast::rmvn(1, mu_j, V_j)
  }

  # 2. Update sigma_j
  rate_sigma <- b_sigma + .5 * colSums((X - Eta %*% t(Lambda))^2)
  sigma2_inv <- stats::rgamma(p, shape_sigma, rate_sigma)

  # 3. Update z
  z <- integer(k)
  log_norm <- mvnfast::dmvn(
    t(Lambda),
    rep(0, p),
    norm_covariance,
    log = TRUE
  )
  log_t <- mvnfast::dmvt(
    t(Lambda),
    rep(0, p),
    t_covariance,
    df = 2 * a_theta,
    log = TRUE
  )
  log_omega <- log(omega)

  hh <- seq_len(k)
  for (h in hh) {
    log_p <- log_omega + ifelse(hh <= h, log_norm[h], log_t[h])
    z[h] <- sample.int(k, 1, prob = exp(log_p - max(log_p)))
  }

  # 4. Update v and omega
  v <- numeric(k)
  for (l in seq_len(k - 1)) {
    v[l] <- stats::rbeta(1, 1 + sum(z == l), alpha + sum(z > l))
  }
  v[k] <- 1
  omega <- stick_break(v)

  # 5. Update theta
  lambda2_sum <- colSums(Lambda^2)
  for (h in seq_len(k)) {
    theta[h] <- ifelse(
      z[h] <= h,
      theta_inf,
      1 / stats::rgamma(1, shape_theta, b_theta + .5 * lambda2_sum[h])
    )
  }

  list(
    Lambda = Lambda,
    sigma2_inv = sigma2_inv,
    theta = theta,
    omega = omega,
    z = z
  )
}
