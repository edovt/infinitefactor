update_reg_params <- function(
  reg_params,
  y,
  Z,
  Eta,
  beta_var,
  Omega_var,
  alpha_var,
  Delta_var,
  shape_sigmay,
  dup_matrix,
  n,
  q,
  k
) {
  beta <- reg_params$beta
  sigma2_y <- reg_params$sigma2_y

  Omega <- reg_params$Omega
  inter <- !is.null(Omega)
  quad_term_Omega <- if (inter) rowSums((Eta %*% Omega) * Eta) else 0

  alpha <- reg_params$alpha
  covs <- !is.null(alpha)
  term_alpha <- if (covs) Z %*% alpha else 0

  Delta <- reg_params$Delta
  covs_inter <- !is.null(Delta)
  quad_term_Delta <- if (covs_inter) rowSums((Eta %*% Delta) * Z) else 0

  # 1. Update beta
  r_beta <- y - quad_term_Omega - term_alpha - quad_term_Delta
  V_beta <- solve(crossprod(Eta) / sigma2_y + diag(1 / beta_var, nrow = k))
  mu_beta <- V_beta %*% (crossprod(Eta, r_beta) / sigma2_y)
  beta <- as.numeric(mvnfast::rmvn(1, mu_beta, V_beta))
  term_beta <- Eta %*% beta

  # 2. Update alpha
  if (covs) {
    r_alpha <- y - term_beta - quad_term_Omega - quad_term_Delta
    V_alpha <- solve(crossprod(Z) / sigma2_y + diag(1 / alpha_var, nrow = q))
    mu_alpha <- V_alpha %*% (crossprod(Z, r_alpha) / sigma2_y)
    alpha <- as.numeric(mvnfast::rmvn(1, mu_alpha, V_alpha))
    term_alpha <- Z %*% alpha
  }

  # 3. Update Delta
  if (covs_inter) {
    r_Delta <- y - term_beta - term_alpha - quad_term_Omega
    W <- t(sapply(1:n, \(i) t(kronecker(Z[i, ], Eta[i, ]))))
    V_Delta <- solve(
      crossprod(W) / sigma2_y + diag(1 / Delta_var, nrow = k * q)
    )
    mu_Delta <- V_Delta %*% (crossprod(W, r_Delta) / sigma2_y)
    Delta <- matrix(as.numeric(mvnfast::rmvn(1, mu_Delta, V_Delta)), k, q)
    quad_term_Delta <- rowSums((Eta %*% Delta) * Z)
  }

  # 4. Update Omega
  if (inter) {
    r_Omega <- y - term_beta - term_alpha - quad_term_Delta
    H <- t(sapply(1:n, \(i) t(kronecker(Eta[i, ], Eta[i, ]))))
    HD <- H %*% dup_matrix
    V_Omega <- solve(
      crossprod(HD) / sigma2_y + diag(1 / Omega_var, nrow = k * (k + 1) / 2)
    )
    mu_Omega <- V_Omega %*% (crossprod(HD, r_Omega) / sigma2_y)
    vech_Omega <- as.numeric(mvnfast::rmvn(1, mu_Omega, V_Omega))
    Omega <- vech2sym(vech_Omega)
    quad_term_Omega <- rowSums((Eta %*% Omega) * Eta)
  }

  # 5. Update sigma2_y
  r <- y - term_beta - quad_term_Omega - term_alpha - quad_term_Delta
  rate_sigmay <- 0.5 + 0.5 * as.numeric(crossprod(r))
  sigma2_y <- 1 / stats::rgamma(1, shape_sigmay, rate_sigmay)

  # 6. Return
  out <- list(sigma2_y = sigma2_y, beta = beta)
  if (inter) {
    out$Omega <- Omega
  }
  if (covs) {
    out$alpha <- alpha
  }
  if (covs_inter) {
    out$Delta <- Delta
  }
  out
}
