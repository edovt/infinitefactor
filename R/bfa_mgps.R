#' Bayesian Factor Analysis with an MGPS prior
#'
#' @param X Data matrix
#' @param iter_warmup Number of warmup iterations to run per chain.
#' @param iter_sampling Number of post-warmup iterations to run per chain.
#' @param adapt (logical) Adapt number of factors across samples?
#' @param k_init If `adapt = FALSE`, the persistent number of factors. If
#'     `adapt = TRUE`, the initial number of factors in the sampler. If `NULL`
#'     the number of factors is selected automatically as `3*log(p)` where `p`
#'     is the number of columns of `X`.
#' @param prior_params Parameters of the MGPS prior.
#' @param eps If `adapt = TRUE`, tolerance for column sums of `Lambda` in the
#'     adaptation step.
#' @param verbose (logical) Show progress bar?
#'
#' @returns List with samples for each parameter
#'
#' @export
#' @examples
#' NULL
bfa_mgps <- function(
  X,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt = TRUE,
  k_init = NULL,
  prior_params = list(
    a_sigma = NULL,
    b_sigma = NULL,
    nu = NULL,
    a1 = NULL,
    a2 = NULL,
    alpha0 = NULL,
    alpha1 = NULL
  ),
  eps = 1e-2,
  verbose = TRUE
) {
  # 0. Input checks --------------------------------------------------------

  # 1. Extract dimensions and hyperparameters ------------------------------
  n <- nrow(X)
  p <- ncol(X)
  sd_X <- apply(X, 2, stats::sd)
  scale_mat <- outer(sd_X, sd_X)
  X <- scale(X)

  k <- ifelse(is.null(k_init), as.integer(3 * log(p)), k_init)
  a_sigma <- ifelse(is.null(prior_params$a_sigma), 1, prior_params$a_sigma)
  b_sigma <- ifelse(is.null(prior_params$b_sigma), 1, prior_params$b_sigma)
  nu <- ifelse(is.null(prior_params$nu), 3, prior_params$nu)
  a1 <- ifelse(is.null(prior_params$a1), 2.1, prior_params$a1)
  a2 <- ifelse(is.null(prior_params$a2), 3.1, prior_params$a2)
  if (adapt) {
    alpha0 <- ifelse(is.null(prior_params$alpha0), -1, prior_params$alpha0)
    alpha1 <- ifelse(is.null(prior_params$alpha1), -0.0005, prior_params$alpha1)
  }

  # 2. Set-up storage of samples -------------------------------------------
  if (!adapt) {
    samples <- list(
      Lambda = array(NA, dim = c(p, k, iter_sampling)),
      Eta = array(NA, dim = c(n, k, iter_sampling)),
      sigma2_inv = matrix(NA, nrow = p, ncol = iter_sampling),
      Phi = array(NA, dim = c(p, k, iter_sampling)),
      delta = matrix(NA, nrow = k, ncol = iter_sampling),
      tau = matrix(NA, nrow = k, ncol = iter_sampling)
    )
  } else {
    samples <- list(
      k = numeric(iter_sampling),
      Omega = array(NA, dim = c(p, p, iter_sampling))
    )
  }

  # 3. Initialize parameters -----------------------------------------------
  Lambda <- matrix(stats::rnorm(p * k), nrow = p, ncol = k)
  Eta <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  sigma2_inv <- stats::rgamma(p, shape = a_sigma, rate = b_sigma)
  Phi <- matrix(
    stats::rgamma(p * k, shape = nu / 2, rate = nu / 2),
    nrow = p,
    ncol = k
  )
  delta <- c(stats::rgamma(1, a1, 1), stats::rgamma(k - 1, a2, 1))
  tau <- cumprod(delta)

  # 4. Gibbs sampler -------------------------------------------------------
  total_iter <- iter_warmup + iter_sampling
  shape_sigma <- a_sigma + n / 2 # constant across iterations
  shape_phi <- (nu + 1) / 2 # constant across iterations
  p_bar <- utils::txtProgressBar(max = total_iter, style = 3)

  for (iter in 1:total_iter) {
    # 4.1 Sample Lambda row-wise
    Eta_cross <- crossprod(Eta) # k×k: η'η
    EtaTX <- crossprod(Eta, X) # k×p: η'X
    for (j in 1:p) {
      D_j <- diag(Phi[j, ] * tau, nrow = k)
      V_j <- solve(D_j + sigma2_inv[j] * Eta_cross)
      mu_j <- sigma2_inv[j] * (V_j %*% EtaTX[, j])
      Lambda[j, ] <- mvtnorm::rmvnorm(1, mu_j, V_j)
    }

    # 4.2 Sample sigma_j
    rate_sigma <- b_sigma + .5 * colSums((X - Eta %*% t(Lambda))^2)
    sigma2_inv <- stats::rgamma(p, shape_sigma, rate_sigma)

    # 4.3 Sample Eta row-wise
    LaTSinv <- t(Lambda * sigma2_inv)
    V_eta <- solve(diag(k) + LaTSinv %*% Lambda)
    mu_eta <- V_eta %*% LaTSinv %*% t(X)
    for (i in 1:n) {
      Eta[i, ] <- mvtnorm::rmvnorm(1, mu_eta[, i], V_eta)
    }

    # 4.4 Sample Phi
    rate_phi <- (nu + t(t(Lambda^2) * tau)) / 2
    Phi <- matrix(
      stats::rgamma(p * k, shape_phi, rate_phi),
      nrow = p,
      ncol = k
    )

    # 4.5 Sample delta
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

    # 4.6 If k is not fixed, update
    if (adapt) {
      p_adapt <- exp(alpha0 + alpha1 * iter)
      unif <- stats::runif(1)
      if (unif < p_adapt) {
        redundant <- (colSums(abs(Lambda) < eps) == p)
        n_redundant <- sum(redundant)
        if (n_redundant == 0) {
          k <- k + 1
          Eta <- cbind(Eta, stats::rnorm(n))
          Phi <- cbind(Phi, stats::rgamma(p, nu / 2, nu / 2))
          delta[k] <- stats::rgamma(1, a2, 1)
          tau <- cumprod(delta)
          Lambda <- cbind(
            Lambda,
            stats::rnorm(p, 0, 1 / sqrt(Phi[, k] * tau[k]))
          )
        } else {
          k <- max(k - n_redundant, 1)
          Lambda <- Lambda[, !redundant, drop = F]
          Phi <- Phi[, !redundant, drop = F]
          Eta <- Eta[, !redundant, drop = F]
          delta <- delta[!redundant]
          tau <- cumprod(delta)
        }
      }
    }

    # 4.7 Save samples
    if (iter > iter_warmup) {
      c_iter <- iter - iter_warmup
      if (!adapt) {
        samples$Lambda[,, c_iter] <- Lambda
        samples$Eta[,, c_iter] <- Eta
        samples$sigma2_inv[, c_iter] <- sigma2_inv
        samples$Phi[,, c_iter] <- Phi
        samples$delta[, c_iter] <- delta
        samples$tau[, c_iter] <- tau
      } else {
        samples$k[c_iter] <- k
        samples$Omega[,, c_iter] <- (Lambda %*%
          t(Lambda) +
          diag(1 / sigma2_inv)) *
          scale_mat
      }
    }

    if (verbose) utils::setTxtProgressBar(p_bar, iter)
  }

  close(p_bar)
  samples
}
