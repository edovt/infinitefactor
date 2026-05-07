#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma;

// ----------------------------------------------------------------------------
// Row-wise Lambda update.
//
// Posterior is N_k(mu_j, V_j) with
//     V_j^{-1} = D_j^{-1} + sigma_j^{-2} eta^T eta,
//     mu_j     = V_j sigma_j^{-2} eta^T x^{(j)},
// where D_j is the prior covariance (diagonal) of lambda_j and varies by
// prior. We do one Cholesky per j and triangular solves for both the mean
// and the noise:
//     R^T R = V_j^{-1}        (R upper-tri)
//     mu = R^{-1} R^{-T} (sigma_j^{-2} eta^T x^{(j)})
//     lambda_j = mu + R^{-1} z,   z ~ N(0, I_k).
//
// What changes per prior is just the diagonal of D_j^{-1}:
//     MGPS:  D_j^{-1}[h] = Phi[j,h] * tau[h]
//     DL:    D_j^{-1}[h] = 1 / (Psi[j,h] * Delta[j,h]^2)
//     CUSP:  D^{-1}[h]   = 1 / theta[h]               (shared across j)
// ----------------------------------------------------------------------------

// Internal helper: sample one row of Lambda given the prior-precision
// diagonal d_inv.   Not exported to R.
//   d_inv  : prior precision diagonal, length k  (entries of D_j^{-1})
//   eta2   : eta^T eta, k x k, shared across j
//   EtaTX  : eta^T X, k x p, shared across j
//   j      : row index (0-based)
//   s_inv  : sigma_j^{-2} (scalar)
arma::rowvec sample_lambda_row(
    const arma::vec& d_inv,
    const arma::mat& eta2,
    const arma::mat& EtaTX,
    int j,
    double s_inv
) {
  int k = d_inv.n_elem;
  mat P  = diagmat(d_inv) + s_inv * eta2;          // posterior precision
  mat R  = chol(P);                                // R upper-tri, R^T R = P
  vec b  = s_inv * EtaTX.col(j);
  vec tmp = solve(trimatl(R.t()), b);              // R^T tmp = b
  vec mu  = solve(trimatu(R), tmp);                // R  mu  = tmp   =>   mu = P^{-1} b
  vec lam = mu + solve(trimatu(R), randn<vec>(k));  // lambda_j ~ N(mu, P^{-1})
  return lam.t();
}

// ----------------------------------------------------------------------------
// MGPS:  D_j^{-1}[h] = Phi[j,h] * tau[h]
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat update_Lambda_mgps(
    const arma::mat& Eta,
    const arma::mat& X,
    const arma::vec& sigma2_inv,
    const arma::mat& Phi,
    const arma::vec& tau
) {
  int p = X.n_cols;
  int k = Eta.n_cols;
  mat Lambda(p, k);
  mat eta2  = Eta.t() * Eta;
  mat EtaTX = Eta.t() * X;

  for (int j = 0; j < p; ++j) {
    vec d_inv = Phi.row(j).t() % tau;
    Lambda.row(j) = sample_lambda_row(d_inv, eta2, EtaTX, j, sigma2_inv(j));
  }
  return Lambda;
}

// ----------------------------------------------------------------------------
// DL (alternative parametrisation):  D_j^{-1}[h] = 1 / (Psi[j,h] * Delta[j,h]^2)
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat update_Lambda_dl(
    const arma::mat& Eta,
    const arma::mat& X,
    const arma::vec& sigma2_inv,
    const arma::mat& Psi,
    const arma::mat& Delta
) {
  int p = X.n_cols;
  int k = Eta.n_cols;
  mat Lambda(p, k);
  mat eta2  = Eta.t() * Eta;
  mat EtaTX = Eta.t() * X;

  for (int j = 0; j < p; ++j) {
    vec d_inv = 1.0 / (Psi.row(j).t() % square(Delta.row(j).t()));
    Lambda.row(j) = sample_lambda_row(d_inv, eta2, EtaTX, j, sigma2_inv(j));
  }
  return Lambda;
}

// ----------------------------------------------------------------------------
// CUSP:  D^{-1}[h] = 1 / theta[h]   (shared across j; built once)
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat update_Lambda_cusp(
    const arma::mat& Eta,
    const arma::mat& X,
    const arma::vec& sigma2_inv,
    const arma::vec& theta
) {
  int p = X.n_cols;
  int k = Eta.n_cols;
  mat Lambda(p, k);
  mat eta2  = Eta.t() * Eta;
  mat EtaTX = Eta.t() * X;
  vec d_inv = 1.0 / theta;          // shared across j

  for (int j = 0; j < p; ++j) {
    Lambda.row(j) = sample_lambda_row(d_inv, eta2, EtaTX, j, sigma2_inv(j));
  }
  return Lambda;
}
