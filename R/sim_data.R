#' Generate data for Factor Analysis examples
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

sim_data_blfr <- function() {}
