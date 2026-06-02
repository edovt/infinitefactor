latent_update_bfa <- function(X, Lambda, sigma2_inv, n, k) {
  LaTSinv <- t(Lambda * sigma2_inv)
  V_eta <- solve(diag(k) + LaTSinv %*% Lambda)
  mu_eta <- t(V_eta %*% LaTSinv %*% t(X))
  U <- chol(V_eta)
  mu_eta + matrix(stats::rnorm(n * k), n, k) %*% U
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
  beta_i
) {
  c1 <- LaTSinv %*% fr_i
  c2 <- (rr_i / sigma2_y) * (beta_i + 2 * OEta_i)
  c1 - Eta_i + c2
}

latent_update_blfr <- function(
  X,
  y,
  Z,
  Eta,
  reg_params,
  Lambda,
  sigma2_inv,
  n,
  k,
  mala_eps,
  n_acceptances
) {
  beta <- reg_params$beta
  Omega <- reg_params$Omega
  alpha <- reg_params$alpha
  Delta <- reg_params$Delta
  sigma2_y <- reg_params$sigma2_y
  Eta_new <- matrix(0, nrow = n, ncol = k)

  if (!is.null(Omega)) {
    # MALA step
    LaTSinv <- t(Lambda * sigma2_inv)
    OEta <- tcrossprod(Omega, Eta)
    LEta <- tcrossprod(Lambda, Eta)
    y_adj <- if (is.null(alpha)) y else y - Z %*% alpha

    for (i in 1:n) {
      Eta_i <- Eta[i, ]
      beta_i <- if (is.null(Delta)) beta else beta + Delta %*% Z[i, ]

      fr_i <- X[i, ] - LEta[, i]
      rr_i <- as.numeric(
        y_adj[i] - crossprod(beta_i, Eta_i) - t(Eta_i) %*% OEta[, i]
      )
      gradient_i <- mala_gradient(
        Eta_i,
        OEta[, i],
        fr_i,
        rr_i,
        LaTSinv,
        sigma2_y,
        beta_i
      )

      mu_prop <- Eta_i + .5 * mala_eps * gradient_i
      Eta_prop <- mu_prop + sqrt(mala_eps) * stats::rnorm(k)
      fr_p <- X[i, ] - Lambda %*% Eta_prop
      rr_p <- as.numeric(
        y_adj[i] -
          crossprod(beta_i, Eta_prop) -
          t(Eta_prop) %*% Omega %*% Eta_prop
      )
      gradient_p <- mala_gradient(
        Eta_prop,
        Omega %*% Eta_prop,
        fr_p,
        rr_p,
        LaTSinv,
        sigma2_y,
        beta_i
      )

      log_ratio <- mala_logdens(Eta_prop, fr_p, rr_p, sigma2_inv, sigma2_y) -
        mala_logdens(Eta_i, fr_i, rr_i, sigma2_inv, sigma2_y) +
        (-.5 / mala_eps) *
          sum((Eta_i - Eta_prop - .5 * mala_eps * gradient_p)^2) -
        (-.5 / mala_eps) *
          sum((Eta_prop - Eta_i - .5 * mala_eps * gradient_i)^2)
      if (log(stats::runif(1)) < log_ratio) {
        Eta_new[i, ] <- Eta_prop
        n_acceptances[i] <- n_acceptances[i] + 1
      } else {
        Eta_new[i, ] <- Eta[i, ]
      }
    }
  } else {
    # Conjugate Gibbs: Various cases
    LaTSinv <- t(Lambda * sigma2_inv)

    if (is.null(alpha) && is.null(Delta)) {
      V_eta <- solve(tcrossprod(beta) / sigma2_y + LaTSinv %*% Lambda + diag(k))
      mu_eta <- t(V_eta %*% (tcrossprod(beta, y / sigma2_y) + LaTSinv %*% t(X))) # n × k
      U_eta <- chol(V_eta)
      Eta_new <- mu_eta + matrix(stats::rnorm(n * k), n, k) %*% U_eta
    } else if (!is.null(alpha) && is.null(Delta)) {
      V_eta <- solve(tcrossprod(beta) / sigma2_y + LaTSinv %*% Lambda + diag(k))
      y_adj <- y - Z %*% alpha
      mu_eta <- t(
        V_eta %*% (tcrossprod(beta, y_adj / sigma2_y) + LaTSinv %*% t(X))
      ) # n × k
      U_eta <- chol(V_eta)
      Eta_new <- mu_eta + matrix(stats::rnorm(n * k), n, k) %*% U_eta
    } else if (!is.null(alpha) && !is.null(Delta)) {
      # Mean and variance changes for every observation
      for (i in seq_len(n)) {
        beta_i <- beta + Delta %*% Z[i, ]
        V_eta_i <- solve(
          tcrossprod(beta_i) / sigma2_y + LaTSinv %*% Lambda + diag(k)
        )
        mu_eta_i <- V_eta_i %*%
          ((y[i] - crossprod(alpha, Z[i, ])) /
            sigma2_y *
            beta_i +
            LaTSinv %*% X[i, ])
        Eta_new[i, ] <- mvnfast::rmvn(1, mu_eta_i, V_eta_i)
      }
    } else {
      stop("Case of null alpha but non-null Delta not implemented yet")
    }
  }

  list(Eta = Eta_new, n_acceptances = n_acceptances)
}
