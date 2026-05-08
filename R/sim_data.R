#' Generate Data for Factor Analysis Examples
#'
#' @description
#' Data generation as in Bhattacharya & Dunson (2011)
#'
#' @details
#' Data is generated as in Bhattacharya & Dunson (2011): number of non-zero
#' columns of the matrix of loadings is linearly decreasing from 2k down to
#' k + 1 and selected at random. Then non-zero elements are sampled
#' from a N(0, 9) distribution. Note that this is a simulation from any
#' "true" model.
#'
#' @param n Number of observations
#' @param p Number of covariates
#' @param k Number of latent factors
#'
#' @returns Named list with the following elements:
#'
#' * `X`: \eqn{n \times p} matrix of observations
#' * `Lambda`: Simulated factor loadings
#' * `sigma2`: Vector of intrinsic variances for each covariate
#' * `Sigma_X_true`: True covariance matrix of `X`, i.e.
#'     \eqn{\Lambda^T\lambda + \text{diag}(\sigma^2)}
#' * `k_true`: True number of latent factors
#'
#' @export
#' @examples
#' set.seed(219)
#'
#' # 0. Data simulation as in Bhattacharya & Dunson (2011)
#' sim <- sim_data_bfa(n = 200, p = 20, k = 3)
#'
#' # 1. Adaptive number of factors, starting at floor(3*log(p)) = 8
#' fit_mgps <- bfa_mgps(sim$X)
sim_data_bfa <- function(n, p, k) {
  # Number of non-zeros per column, linearly decreasing from 2k down to k+1
  n_nonzero <- round(seq(2 * k, k + 1, length.out = k))
  Lambda <- matrix(0, nrow = p, ncol = k)
  for (h in seq_len(k)) {
    rows <- sample(p, n_nonzero[h])
    Lambda[rows, h] <- stats::rnorm(n_nonzero[h], mean = 0, sd = 3)
  }

  sigma2 <- 1 / stats::rgamma(p, shape = 1, rate = 0.25)
  Sigma_X_true <- tcrossprod(Lambda) + diag(sigma2)
  X <- mvnfast::rmvn(n, rep(0, p), Sigma_X_true)

  list(
    X = X,
    Lambda = Lambda,
    sigma2 = sigma2,
    Sigma_X_true = Sigma_X_true,
    k_true = k
  )
}

#' Generate Data for Latent Regression Examples
#'
#' @details
#' Data is generated as follows:
#' * Latent level as in Bhattacharya & Dunson (2011): number of non-zero
#' columns of the matrix of loadings is linearly decreasing from 2k down to
#' k + 1 and selected at random. Then non-zero elements are sampled
#' from a N(0, 9) distribution.
#' * Latent effects: `beta` is drawn from a Normal distribution with mean 0 and
#' standard deviation `beta_sd`. `Omega` elements are drawn from a Normal
#' distribution with mean 0 and sd `omega_sd`, then \eqn{\Omega = \frac{\Omega + \Omega^T}{2}}
#' * Observations are drawn with mean from the corresponding latent regression
#' model and variance `sigma2_y`.
#'
#' @param n Number of observations
#' @param n_test Number of test observations
#' @param p Number of covariates in X
#' @param k Number of latent factors
#' @param interactions (logical) Include interaction effects?
#' @param beta_sd standard deviation of latent main effects
#' @param omega_sd standard deviation of latent interaction effects
#' @param sigma2_y observation level noise in response y
#'
#' @returns Named list with the following elements:
#'
#' * `X`: \eqn{n \times p} matrix of predictors
#' * `y`: \eqn{n}-dimensional vector of the response
#' * `X_test`: \eqn{\text{n\_test} \times p} matrix of test predictors
#' * `y_test`: `n_test`-dimensional vector of test responses
#' * `Lambda`: Simulated factor loadings
#' * `sigma2`: Vector of intrinsic variances for each covariate
#' * `beta`: latent main effects
#' * `Omega`: latent interaction effects, `NULL` if `interactions=FALSE`
#' * `sigma2_y`: response level noise
#' * `Sigma_X_true`: True covariance matrix of `X`, i.e.
#'     \eqn{\Lambda^T\lambda + \text{diag}(\sigma^2)}
#' * `beta_X`: induced main effects in original X scale
#' * `Omega_X`: induced interaction effects in original X scale
#' * `intercept_X`: induced intercept, i.e. \eqn{tr(\Omega V)}
#' * `k_true`: True number of latent factors
#'
#' @export
#' @examples
#' set.seed(219)
#'
#' # 0. Data simulation
#' sim <- sim_data_blfr(n = 200, n_test = 20, p = 20, k = 3)
#'
#' # 1. Adaptive number of factors, starting at floor(3*log(p)) = 8
#' fit_mgps <- bfa_mgps(sim$X)
sim_data_blfr <- function(
  n,
  n_test,
  p,
  k,
  interactions = FALSE,
  beta_sd = 2,
  omega_sd = 1,
  sigma2_y = 0.5
) {
  # 1. First get Lambda as in the bfa case
  n_nonzero <- round(seq(2 * k, k + 1, length.out = k))
  Lambda <- matrix(0, nrow = p, ncol = k)
  for (h in seq_len(k)) {
    rows <- sample(p, n_nonzero[h])
    Lambda[rows, h] <- stats::rnorm(n_nonzero[h], mean = 0, sd = 3)
  }
  sigma2 <- 1 / stats::rgamma(p, shape = 1, rate = 0.25)

  # 2. Regression parameters for latent Eta
  beta <- stats::rnorm(
    k,
    mean = sample(c(-3, 3), k, replace = TRUE),
    sd = beta_sd
  )
  Omega <- if (interactions) {
    M <- matrix(stats::rnorm(k * k, sd = omega_sd), k, k)
    (M + t(M)) / 2
  } else {
    NULL
  }

  # 3. Generate train and test set
  generate <- function(N) {
    eta <- matrix(stats::rnorm(N * k), N, k)
    Xi <- matrix(stats::rnorm(N * p), N, p) %*% diag(sqrt(sigma2))
    X <- eta %*% t(Lambda) + Xi
    quad <- if (interactions) rowSums((eta %*% Omega) * eta) else 0
    y <- as.numeric(eta %*% beta) + quad + stats::rnorm(N, sd = sqrt(sigma2_y))
    list(X = X, y = y, eta = eta)
  }
  train <- generate(n)
  test <- generate(n_test)

  # 4. Induced quantities in original X scale
  Sigma_X_true <- tcrossprod(Lambda) + diag(sigma2)
  V_true <- solve(crossprod(Lambda, Lambda / sigma2) + diag(k))
  A_true <- V_true %*% t(Lambda / sigma2) # k x p, in original scale
  sd_X_true <- sqrt(diag(Sigma_X_true)) # marginal sd of x
  D_inv <- diag(1 / sd_X_true)

  beta_X <- as.numeric(D_inv %*% crossprod(A_true, beta))
  Omega_X <- if (interactions) {
    D_inv %*% crossprod(A_true, Omega) %*% A_true %*% D_inv
  } else {
    NULL
  }
  intercept_X <- if (interactions) sum(diag(Omega %*% V_true)) else 0

  list(
    X = train$X,
    y = train$y,
    X_test = test$X,
    y_test = test$y,
    Lambda = Lambda,
    sigma2 = sigma2,
    beta = beta,
    Omega = Omega,
    sigma2_y = sigma2_y,
    Sigma_X_true = Sigma_X_true,
    beta_X = beta_X, # original X scale
    Omega_X = Omega_X, # original X scale
    intercept_X = intercept_X,
    k_true = k
  )
}
