#' Bayesian Latent Factor Regression with a Dirichlet-Laplace prior
#'
#' @param y Response vector
#' @param X Predictors matrix
#' @param Z optional matrix of covariates
#' @param iter_warmup Number of warmup iterations to run per chain
#' @param iter_sampling Number of post-warmup iterations to run per chain
#' @param k Number of factors
#' @param induced (logical) return induced effects of original predictors?
#' @param interactions (logical) include interactions between latent factors?
#' @param covariates_interactions (logical) include interactions between
#'     covariates and latent factors?
#' @param prior_params Parameters of the Dirichlet-Laplace prior
#' @param mala_eps If `interactions=TRUE`, step-size of the MALA step for
#'     latent factors. If `NULL`, uses `k^(-1/3)` as default.
#' @param adapt_mala_eps (logical) Should `mala_eps` be adaptive?
#' @param window_mala If `adapt_mala_eps = TRUE`, number of iterations during
#'     warmup that are considered for the adaptation of `mala_eps`
#' @param verbose (logical) Show progress bar?
#'
#' @returns List with samples for each parameter
#'
#' @export
#' @examples
#' set.seed(219)
#'
#' # 1. Example without interactions
#' sim <- sim_data_blfr(n = 200, n_test = 100, p = 20, k = 3)
#' fit_dl <- blfr_dl(sim$y, sim$X)
#'
#' # 1.1 Induced effect summaries
#' Sigma_X_plots <- plot_Sigma_X(fit_dl, real = sim$Sigma_X_true)
#' Sigma_X_plots$scatter
#' main_effects_plots <- plot_effects_X(
#'     fit_dl,
#'     real_beta = sim$beta_X
#' )
#' main_effects_plots$main
#' main_effects_plots$scatter
#'
#' # 1.2 Match Align
#' aligned <- match_align(fit_dl)
#' plot_match_align(aligned, "Lambda", "mean")
#'
#' # 1.3 Estimates on test set
#' preds <- predict_blfr(fit_dl, sim$X_test)
#' mean((preds - sim$y_test)^2)
#' plot(1:length(preds), preds - sim$y_test)
#'
#' # 2. For example with interactions, see help(blfr_mgps)
#'
blfr_dl <- function(
  y,
  X,
  Z = NULL,
  iter_warmup = 1000,
  iter_sampling = 1000,
  k = NULL,
  induced = TRUE,
  interactions = FALSE,
  covariates_interactions = FALSE,
  prior_params = list(
    nu_y = NULL,
    beta_var = NULL,
    Omega_var = NULL,
    alpha_var = NULL,
    Delta_var = NULL,
    a_sigma = NULL,
    b_sigma = NULL,
    a = NULL
  ),
  mala_eps = NULL,
  adapt_mala_eps = TRUE,
  window_mala = floor(iter_warmup / 10),
  verbose = TRUE
) {
  # 0. Input checks --------------------------------------------------------
  if (!is.matrix(X)) {
    X <- as.matrix(X)
    if (!is.numeric(X)) {
      stop("'X' must be a numeric matrix.")
    }
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop("'X' must not contain NA, NaN, or Inf values.")
  }
  if (nrow(X) <= 1 || ncol(X) < 2) {
    stop("'X' must have at least 2 rows and 2 columns.")
  }
  if (any(apply(X, 2, stats::sd) == 0)) {
    stop("'X' has constant columns; remove or impute them before fitting.")
  }

  if (!is.numeric(y) || !is.null(dim(y))) {
    stop("'y' must be a numeric vector.")
  }
  if (length(y) != nrow(X)) {
    stop("'y' must have length equal to nrow(X).")
  }
  if (anyNA(y) || any(!is.finite(y))) {
    stop("'y' must not contain NA, NaN, or Inf values.")
  }
  if (!is.null(Z)) {
    if (!is.matrix(Z)) {
      Z <- as.matrix(Z)
      if (!is.numeric(Z)) {
        stop("'Z' must be a numeric matrix.")
      }
    }
    if (anyNA(Z) || any(!is.finite(Z))) {
      stop("'Z' must not contain NA, NaN, or Inf values.")
    }
  }

  if (!is.null(k)) {
    if (
      !is.numeric(k) ||
        length(k) != 1 ||
        k != as.integer(k) ||
        k < 1
    ) {
      stop("'k' must be a positive integer.")
    }
    if (k >= min(nrow(X), ncol(X))) {
      stop("'k' must be less than min(n, p).")
    }
  }

  if (
    !is.numeric(iter_warmup) ||
      length(iter_warmup) != 1 ||
      iter_warmup != as.integer(iter_warmup) ||
      iter_warmup < 1
  ) {
    stop("'iter_warmup' must be a positive integer.")
  }
  if (
    !is.numeric(iter_sampling) ||
      length(iter_sampling) != 1 ||
      iter_sampling != as.integer(iter_sampling) ||
      iter_sampling < 1
  ) {
    stop("'iter_sampling' must be a positive integer.")
  }

  if (!is.logical(induced) || length(induced) != 1) {
    stop("'induced' must be a single logical value.")
  }
  if (!is.logical(interactions) || length(interactions) != 1) {
    stop("'interactions' must be a single logical value.")
  }
  if (
    !is.logical(covariates_interactions) || length(covariates_interactions) != 1
  ) {
    stop("'covariates_interactions' must be a single logical value.")
  }
  if (!is.logical(adapt_mala_eps) || length(adapt_mala_eps) != 1) {
    stop("'adapt_mala_eps' must be a single logical value.")
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be a single logical value.")
  }

  if (!is.null(mala_eps)) {
    if (!is.numeric(mala_eps) || length(mala_eps) != 1 || mala_eps <= 0) {
      stop("'mala_eps' must be a single positive number.")
    }
  }
  if (interactions && adapt_mala_eps) {
    if (
      !is.numeric(window_mala) ||
        length(window_mala) != 1 ||
        window_mala != as.integer(window_mala) ||
        window_mala < 1
    ) {
      stop("'window_mala' must be a positive integer.")
    }
  }

  nu_y_check <- prior_params$nu_y %||% 1
  beta_var_check <- prior_params$beta_var %||% 100
  Omega_var_check <- prior_params$Omega_var %||% 100
  alpha_var_check <- prior_params$alpha_var %||% 100
  Delta_var_check <- prior_params$Delta_var %||% 100
  a_sigma_check <- prior_params$a_sigma %||% 1
  b_sigma_check <- prior_params$b_sigma %||% 1
  a_check <- prior_params$a %||% 1
  if (nu_y_check <= 0) {
    stop("'nu_y' must be positive.")
  }
  if (beta_var_check <= 0) {
    stop("'beta_var' must be positive.")
  }
  if (interactions && Omega_var_check <= 0) {
    stop("'Omega_var' must be positive.")
  }
  if (!is.null(Z)) {
    if (alpha_var_check <= 0) {
      stop("'alpha_var' must be positive")
    }
    if (covariates_interactions && Delta_var_check <= 0) {
      stop("'Delta_var' must be positive.")
    }
  }
  if (a_sigma_check <= 0 || b_sigma_check <= 0) {
    stop("'a_sigma' and 'b_sigma' must be positive.")
  }
  if (a_check <= 0) {
    stop("'a' must be positive.")
  }

  # 1. Extract dimensions and hyperparameters ------------------------------
  n <- nrow(X)
  p <- ncol(X)
  sd_X <- apply(X, 2, stats::sd)
  scale_mat <- outer(sd_X, sd_X)
  X <- scale(X)
  covariates <- !is.null(Z)
  q <- if (covariates) ncol(Z) else NULL

  if (!covariates && covariates_interactions) {
    stop("Case of null alpha and non-null Delta not implemented")
  }

  k <- k %||% as.integer(3 * log(p))
  mala_eps <- mala_eps %||% k^(-1 / 3)
  nu_y <- prior_params$nu_y %||% 1
  beta_var <- prior_params$beta_var %||% 100
  Omega_var <- prior_params$Omega_var %||% 100
  alpha_var <- prior_params$alpha_var %||% 100
  Delta_var <- prior_params$Delta_var %||% 100
  a_sigma <- prior_params$a_sigma %||% 1
  b_sigma <- prior_params$b_sigma %||% 1
  a <- prior_params$a %||% (1 / 2)

  # 2. Set-up storage ------------------------------------------------------
  samples <- list(
    sigma2_y = numeric(iter_sampling),
    beta = matrix(NA, k, iter_sampling),
    Lambda = array(NA, c(p, k, iter_sampling)),
    Eta = array(NA, c(n, k, iter_sampling)),
    sigma2_inv = matrix(NA, p, iter_sampling),
    Delta_dl = array(NA, dim = c(p, k, iter_sampling)),
    Psi = array(NA, dim = c(p, k, iter_sampling)),
    Sigma_X = array(NA, dim = c(p, p, iter_sampling))
  )
  if (induced) {
    samples$beta_X <- matrix(NA, p, iter_sampling)
  }
  if (interactions) {
    samples$Omega <- array(NA, c(k, k, iter_sampling))
  }
  if (induced && interactions) {
    samples$Omega_X <- array(NA, c(p, p, iter_sampling))
    samples$intercept_X <- numeric(iter_sampling)
  }
  if (covariates) {
    samples$alpha <- matrix(NA, q, iter_sampling)
  }
  if (covariates && covariates_interactions) {
    samples$Delta <- array(NA, c(k, q, iter_sampling))
  }
  if (covariates && covariates_interactions && induced) {
    samples$Delta_X <- array(NA, c(p, q, iter_sampling))
  }

  # 3. Initialize parameters -----------------------------------------------
  dl_params <- list(
    Lambda = matrix(stats::rnorm(p * k), p, k),
    sigma2_inv = stats::rgamma(p, a_sigma, b_sigma),
    Delta_dl = matrix(stats::rgamma(p * k, a, 1 / 2), p, k),
    Psi = matrix(stats::rexp(p * k, 1 / 2), p, k)
  )
  Eta <- matrix(stats::rnorm(n * k), n, k)
  reg_params <- list(
    beta = stats::rnorm(k, sd = sqrt(beta_var)),
    sigma2_y = 1 / stats::rgamma(1, nu_y / 2, nu_y / 2)
  )
  if (interactions) {
    reg_params$Omega <- vech2sym(stats::rnorm(
      k * (k + 1) / 2,
      sd = sqrt(Omega_var)
    ))
  }
  if (covariates) {
    reg_params$alpha <- stats::rnorm(q, sd = sqrt(alpha_var))
    if (covariates_interactions) {
      reg_params$Delta <- matrix(
        stats::rnorm(q * k, sd = sqrt(Delta_var)),
        k,
        q
      )
    }
  }

  # 4. Gibbs sampler -------------------------------------------------------
  total_iter <- iter_warmup + iter_sampling
  n_acceptances <- numeric(n)
  if (verbose) {
    p_bar <- utils::txtProgressBar(max = total_iter, style = 3)
  }

  # Constants across iterations
  shape_sigmay <- (nu_y + n) / 2
  shape_sigma <- a_sigma + n / 2
  sd_X_inv <- diag(1 / sd_X)
  dup_matrix <- matrixcalc::duplication.matrix(k)

  for (iter in 1:total_iter) {
    # 4.1 Update parameters
    latent_update <- latent_update_blfr(
      X,
      y,
      Z,
      Eta,
      reg_params,
      dl_params$Lambda,
      dl_params$sigma2_inv,
      n,
      k,
      mala_eps,
      n_acceptances
    )
    Eta <- latent_update$Eta
    n_acceptances <- latent_update$n_acceptances

    reg_params <- update_reg_params(
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
    )
    dl_params <- DL_update(dl_params, X, Eta, n, p, k, shape_sigma, b_sigma, a)

    # 4.2 Adapt stepsize if wanted during warm-up and after window
    if (adapt_mala_eps && (iter %% window_mala == 0) && iter <= iter_warmup) {
      acceptance_mean <- mean(n_acceptances) / window_mala
      mala_eps <- exp(log(mala_eps) + acceptance_mean - 0.574)
      n_acceptances <- numeric(n)
    }

    # 4.3 Save samples
    if (iter > iter_warmup) {
      c_iter <- iter - iter_warmup
      samples$Sigma_X[,, c_iter] <- (dl_params$Lambda %*%
        t(dl_params$Lambda) +
        diag(1 / dl_params$sigma2_inv)) *
        scale_mat
      samples$sigma2_y[c_iter] <- reg_params$sigma2_y
      samples$beta[, c_iter] <- reg_params$beta
      samples$Lambda[,, c_iter] <- dl_params$Lambda
      samples$Eta[,, c_iter] <- Eta
      samples$sigma2_inv[, c_iter] <- dl_params$sigma2_inv
      samples$Delta_dl[,, c_iter] <- dl_params$Delta_dl
      samples$Psi[,, c_iter] <- dl_params$Psi

      if (interactions) {
        samples$Omega[,, c_iter] <- reg_params$Omega
      }
      if (covariates) {
        samples$alpha[, c_iter] <- reg_params$alpha
        if (covariates_interactions) {
          samples$Delta[,, c_iter] <- reg_params$Delta
        }
      }
      if (induced) {
        L <- dl_params$Lambda
        V <- solve(t(L) %*% diag(dl_params$sigma2_inv) %*% L + diag(k))
        A <- V %*% t(L) %*% diag(dl_params$sigma2_inv)
        samples$beta_X[, c_iter] <- sd_X_inv %*% t(A) %*% reg_params$beta
        if (interactions) {
          samples$Omega_X[,, c_iter] <- sd_X_inv %*%
            t(A) %*%
            reg_params$Omega %*%
            A %*%
            sd_X_inv
          samples$intercept_X[c_iter] <- sum(diag(reg_params$Omega %*% V))
        }
        if (covariates && covariates_interactions) {
          samples$Delta_X[,, c_iter] <- sd_X_inv %*% t(A) %*% reg_params$Delta
        }
      }
    }

    if (verbose) utils::setTxtProgressBar(p_bar, iter)
  }

  if (verbose) {
    close(p_bar)
  }
  samples
}
