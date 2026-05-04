eta_star <- function(Eta, k) {
  m <- k * (k + 1L) / 2L
  Eta_star <- matrix(0, nrow(Eta), m)
  idx <- 1L
  for (l in seq_len(k)) {
    for (h in seq_len(l)) {
      Eta_star[, idx] <- if (h == l) Eta[, h]^2 else 2 * Eta[, h] * Eta[, l]
      idx <- idx + 1L
    }
  }
  Eta_star
}

vech_to_sym <- function(omega_vec, k) {
  Omega <- matrix(0, k, k)
  idx <- 1L
  for (l in seq_len(k)) {
    for (h in seq_len(l)) {
      Omega[h, l] <- omega_vec[idx]
      Omega[l, h] <- omega_vec[idx]
      idx <- idx + 1L
    }
  }
  Omega
}

update_reg_params <- function(
  reg_params,
  y,
  Eta,
  beta_var,
  Omega_var,
  shape_sigmay,
  n,
  k
) {
  beta <- reg_params$beta
  sigma2_y <- reg_params$sigma2_y
  Omega <- reg_params$Omega
  inter <- !is.null(Omega)
  quad_term <- if (inter) rowSums((Eta %*% Omega) * Eta) else 0

  # 1. Update beta
  resid_beta <- y - quad_term
  V_beta <- solve(crossprod(Eta) / sigma2_y + diag(1 / beta_var, nrow = k))
  mu_beta <- V_beta %*% (crossprod(Eta, resid_beta) / sigma2_y)
  beta <- mvnfast::rmvn(1, mu_beta, V_beta)

  # 2. Update Omega
  if (inter) {
    Eta_star <- eta_star(Eta, k) # n x k(k+1)/2
    resid_Omega <- y - Eta %*% beta
    m <- k * (k + 1L) / 2L
    P_Omega <- crossprod(Eta_star) / sigma2_y + diag(1 / Omega_var, nrow = m)
    V_Omega <- solve(P_Omega)
    mu_Omega <- V_Omega %*% (crossprod(Eta_star, resid_Omega) / sigma2_y)
    omega_vec <- mvnfast::rmvn(1, mu_Omega, V_Omega)
    Omega <- vech_to_sym(omega_vec, k)
  }

  # 3. Update sigma2_y
  resid_all <- y - Eta %*% beta - quad_term
  rate_sigmay <- 0.5 + 0.5 * as.numeric(crossprod(resid_all))
  sigma2_y <- 1 / stats::rgamma(1, shape_sigmay, rate_sigmay)

  # 4. Return
  out <- list(sigma2_y = sigma2_y, beta = beta)
  if (inter) {
    out$Omega <- Omega
  }
  out
}
