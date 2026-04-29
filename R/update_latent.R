latent_update_bfa <- function(X, Lambda, sigma2_inv, n, k) {
  LaTSinv <- t(Lambda * sigma2_inv)
  V_eta <- solve(diag(k) + LaTSinv %*% Lambda)
  mu_eta <- V_eta %*% LaTSinv %*% t(X)
  Eta <- matrix(0, nrow = n, ncol = k)
  for (i in 1:n) {
    Eta[i, ] <- mvtnorm::rmvnorm(1, mu_eta[, i], V_eta)
  }
  Eta
}

mala_logdens <- function(Eta, f_r, r_r, sigma2_inv, sigma2_y) {
  c1 <- -r_r^2 / (2 * sigma2_y)
  c2 <- -.5 * sum(sigma2_inv * f_r^2)
  c3 <- -.5 * sum(Eta^2)

  c1 + c2 + c3
}

mala_gradient <- function(
  Eta_i,
  OEta_i,
  fr_i,
  rr_i,
  LaTSinv,
  sigma2_y,
  beta
) {
  c1 <- LaTSinv %*% fr_i
  c2 <- (rr_i / sigma2_y) * (beta + 2 * OEta_i)
  c1 - Eta_i + c2
}

latent_update_blfr <- function(
  X,
  y,
  Eta,
  reg_params,
  Lambda,
  sigma2_inv,
  n,
  k,
  eps
) {
  beta <- reg_params$beta
  sigma2_y <- reg_params$sigma2_y
  Omega <- reg_params$Omega
  Eta_new <- matrix(0, nrow = n, ncol = k)

  if (!is.null(Omega)) {
    # MALA step
    LaTSinv <- t(Lambda * sigma2_inv)
    OEta <- tcrossprod(Omega, Eta)
    LEta <- tcrossprod(Lambda, Eta)

    for (i in 1:n) {
      Eta_i <- Eta[i, ]
      fr_i <- X[i, ] - LEta[, i]
      rr_i <- y[i] - crossprod(beta, Eta_i) - t(Eta_i) %*% OEta[, i]
      gradient_i <- mala_gradient(
        Eta_i,
        OEta[, i],
        fr_i,
        rr_i,
        LaTSinv,
        sigma2_y,
        beta
      )
      mu_prop <- Eta_i + .5 * gradient_i
      Eta_prop <- mu_prop + sqrt(eps) * stats::rnorm(k)
      fr_p <- X[i, ] - Lambda %*% Eta_prop
      rr_p <- y[i] -
        crossprod(beta, Eta_prop) -
        t(Eta_prop) %*% Omega %*% Eta_prop
      gradient_p <- mala_gradient(
        Eta_prop,
        Omega %*% Eta_prop,
        fr_p,
        rr_p,
        LaTSinv,
        sigma2_y,
        beta
      )

      log_ratio <- mala_logdens(Eta_prop, fr_p, rr_p, sigma2_inv, sigma2_y) -
        mala_logdens(Eta_i, fr_i, rr_i, sigma2_inv, sigma2_y) +
        (-.5 / eps) * sum((Eta_i - Eta_prop + .5 * gradient_p)^2) -
        (-.5 / eps) * sum((Eta_prop - Eta_i + .5 * gradient_i)^2)

      Eta_new[i, ] <- ifelse(
        log(stats::runif(1)) < log_ratio,
        Eta_prop,
        Eta[i, ]
      )
    }
  } else {
    # Conjugate Gibbs
    LaTSinv <- t(Lambda * sigma2_inv)
    V_eta <- solve(tcrossprod(beta) / sigma2_y + LaTSinv %*% Lambda + diag(k))
    mu_eta <- V_eta %*% (tcrossprod(beta, y / sigma2_y) + LaTSinv %*% t(X))
    for (i in 1:n) {
      Eta_new[i, ] <- mvtnorm::rmvnorm(1, mu_eta[, i], V_eta)
    }
  }
  Eta_new
}
