#' Summary for adaptive number of factors
#'
#' @param fit Returned samples from any `bfa_*()` with `adapt=TRUE`.
#' @param real Real number of factors to include in summary, `NULL` otherwise.
#'
#' @returns Median number of factors in `fit$k`, also prints summary of the
#'     posterior distribution of `k`.
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
#' k_fit <- summary_k(fit_mgps, real = sim$k_true)
summary_k <- function(fit, real = NULL) {
  if (!("k" %in% names(fit))) {
    stop("Fit with non-adaptive number of factors")
  }

  if (is.null(real)) {
    cat(sprintf(
      "Estimated k -- median: %d, 95%% CI: [%d, %d]\n",
      as.integer(stats::median(fit$k)),
      as.integer(stats::quantile(fit$k, 0.025)),
      as.integer(stats::quantile(fit$k, 0.975)),
    ))
  } else {
    cat(sprintf(
      "Estimated k -- median: %d, 95%% CI: [%d, %d]  (true k = %d)\n",
      as.integer(stats::median(fit$k)),
      as.integer(stats::quantile(fit$k, 0.025)),
      as.integer(stats::quantile(fit$k, 0.975)),
      real
    ))
  }

  stats::median(fit$k)
}

#' Summary Plots of Posterior Covariance
#'
#' @param fit Returned samples from any `bfa_*()` or `blfr_*()` function.
#' @param real Real Sigma_X if known, `NULL` otherwise
#'
#' @returns Named list of plots. If `real=NULL`, `plots$main` is the plot of the
#' mean posterior covariance of X. Otherwise, three plots:
#' * `main`: compares mean posterior covariance with real covariance
#' * `residual`: plot of residuals between mean posterior and real covariance
#' * `scatter`: scatterplot of mean posterior elements vs real elements.
#'
#' @importFrom rlang .data
#' @export
#' @examples
#' set.seed(219)
#'
#' # 0. Data simulation as in Bhattacharya & Dunson (2011)
#' sim <- sim_data_bfa(n = 200, p = 20, k = 3)
#'
#' # 1. Adaptive number of factors, starting at floor(3*log(p)) = 8
#' fit_mgps <- bfa_mgps(sim$X)
#' plots_Sigma_X <- plot_Sigma_X(fit_mgps, real = sim$Sigma_X_true)
#' plots_Sigma_X$main
#' plots_Sigma_X$residual
#' plots_Sigma_X$scatter
plot_Sigma_X <- function(fit, real = NULL) {
  Sigma_X_est <- apply(fit$Sigma_X, c(1, 2), mean)
  p <- nrow(Sigma_X_est)
  to_tile <- function(df) {
    df$row <- p - df$row + 1
    df
  }
  heat_theme <- ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key.height = ggplot2::unit(0.8, "cm")
    )

  if (!is.null(real)) {
    resid <- Sigma_X_est - real
    clim <- max(abs(c(real, Sigma_X_est)))
    rlim <- max(abs(resid))

    df_main <- to_tile(rbind(
      mat_long(real, "True Sigma_X"),
      mat_long(Sigma_X_est, "Estimated Sigma_X")
    ))
    df_main$panel <- factor(
      df_main$panel,
      levels = c("True Sigma_X", "Estimated Sigma_X")
    )
    df_resid <- to_tile(mat_long(resid, "Residual"))

    g_main <- ggplot2::ggplot(
      df_main,
      ggplot2::aes(x = .data$col, y = .data$row, fill = .data$value)
    ) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~ .data$panel) +
      ggplot2::scale_fill_gradient2(
        low = "#3d52bf",
        mid = "white",
        high = "#800000",
        limits = c(-clim, clim)
      ) +
      ggplot2::coord_fixed() +
      heat_theme +
      ggplot2::labs(fill = NULL)

    g_resid <- ggplot2::ggplot(
      df_resid,
      ggplot2::aes(x = .data$col, y = .data$row, fill = .data$value)
    ) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~ .data$panel) +
      ggplot2::scale_fill_gradient2(
        low = "#3d52bf",
        mid = "white",
        high = "#800000",
        limits = c(-rlim, rlim)
      ) +
      ggplot2::coord_fixed() +
      heat_theme +
      ggplot2::labs(fill = NULL)

    idx <- which(upper.tri(real, diag = TRUE))
    lims <- range(c(real[idx], Sigma_X_est[idx]))
    g_scat <- ggplot2::ggplot(
      data.frame(true = real[idx], est = Sigma_X_est[idx]),
      ggplot2::aes(x = .data$true, y = .data$est)
    ) +
      ggplot2::geom_point(alpha = 0.4, size = 0.9) +
      ggplot2::geom_abline(colour = "red", linewidth = 0.8) +
      ggplot2::coord_fixed(xlim = lims, ylim = lims) +
      ggplot2::labs(
        x = "True Sigma_X entries",
        y = "Estimated Sigma_X entries",
        title = "Entry-wise comparison (upper triangle)"
      ) +
      ggplot2::theme_minimal(base_size = 15)

    return(list(main = g_main, residual = g_resid, scatter = g_scat))
  } else {
    clim <- max(abs(Sigma_X_est))
    df_main <- to_tile(rbind(mat_long(Sigma_X_est, "Estimated Sigma_X")))

    g_main <- ggplot2::ggplot(
      df_main,
      ggplot2::aes(x = .data$col, y = .data$row, fill = .data$value)
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        low = "#3d52bf",
        mid = "white",
        high = "#800000",
        limits = c(-clim, clim)
      ) +
      ggplot2::coord_fixed() +
      heat_theme +
      ggplot2::labs(fill = NULL)

    return(list(main = g_main))
  }
}

#' Summary Plots of Induced Main Effects
#'
#' @param fit_blfr Returned samples from any `blfr_*()` function
#' @param real_beta Real vector of main effects, if known
#' @param real_intercept Real intercept, if known
#'
#' @returns Named list of plots
#'
#' @export
#' @examples
#' NULL
plot_effects_X <- function(fit_blfr, real_beta = NULL, real_intercept = NULL) {
  if (!is.null(real_intercept) && is.null(real_beta)) {
    stop("real_beta cannot be NULL if real_intercept is supplied.")
  }
  intercept <- "intercept_X" %in% names(fit_blfr)
  if (intercept) {
    if (!is.null(real_beta) && is.null(real_intercept)) {
      stop(
        "Intercept in fit and real_beta provided, but real_intercept is NULL"
      )
    }
    real <- c(real_intercept, real_beta)
  } else {
    real <- real_beta
  }

  betaX_mean <- rowMeans(fit_blfr$beta_X)
  betaX_q05 <- apply(fit_blfr$beta_X, 1, stats::quantile, probs = 0.025)
  betaX_q95 <- apply(fit_blfr$beta_X, 1, stats::quantile, probs = 0.975)
  if (intercept) {
    betaX_mean <- c(mean(fit_blfr$intercept_X), betaX_mean)
    betaX_q05 <- c(stats::quantile(fit_blfr$intercept_X, 0.025), betaX_q05)
    betaX_q95 <- c(stats::quantile(fit_blfr$intercept_X, 0.975), betaX_q95)
  }
  df <- data.frame(
    idx = seq_along(betaX_mean),
    est = betaX_mean,
    lo = betaX_q05,
    hi = betaX_q95
  )

  if (is.null(real_beta)) {
    g_caterpillar <- ggplot2::ggplot(df, ggplot2::aes(x = .data$idx)) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
        width = 0,
        alpha = 0.5
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = .data$est),
        colour = "#800000",
        size = 1.5
      ) +
      ggplot2::labs(
        x = "Predictor index",
        y = expression(beta[X]),
        title = "Induced regression coefficients (95% CI)"
      ) +
      ggplot2::theme_minimal(base_size = 13)

    return(list(main = g_caterpillar))
  } else {
    df$true <- real
    g_caterpillar <- ggplot2::ggplot(df, ggplot2::aes(x = .data$idx)) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
        width = 0,
        alpha = 0.5
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = .data$est),
        colour = "#800000",
        size = 1.5
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = .data$true),
        colour = "#3d52bf",
        size = 1.5,
        shape = 4
      ) +
      ggplot2::labs(
        x = "Predictor index",
        y = expression(beta[X]),
        title = "Induced regression coefficients (95% CI)",
        caption = "x = true,   . = posterior mean"
      ) +
      ggplot2::theme_minimal(base_size = 13)

    lims <- range(c(df$true, df$lo, df$hi))
    g_scat <- ggplot2::ggplot(df, ggplot2::aes(x = .data$true, y = .data$est)) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
        width = 0,
        alpha = 0.4
      ) +
      ggplot2::geom_point(size = 1.2) +
      ggplot2::geom_abline(colour = "red") +
      ggplot2::coord_fixed(xlim = lims, ylim = lims) +
      ggplot2::labs(
        x = expression("True " * beta[X]),
        y = expression("Estimated " * beta[X])
      ) +
      ggplot2::theme_minimal(base_size = 13)

    return(list(main = g_caterpillar, scatter = g_scat))
  }
}

#' Summary Plots of Induced Interaction Effects
#'
#' @param fit_blfr Returned samples from any `blfr_*()` function with interactions
#' @param real Real matrix of induced interaction effects, if known
#'
#' @returns Named list of plots
#'
#' @export
#' @examples
#' NULL
plot_Omega_X <- function(fit_blfr, real = NULL) {
  if (!("Omega_X" %in% names(fit_blfr))) {
    stop("Model not fitted with interactions")
  }
  Omega_X_est <- apply(fit_blfr$Omega_X, c(1, 2), mean)
  p <- nrow(Omega_X_est)
  to_tile <- function(df) {
    df$row <- p - df$row + 1
    df
  }
  heat_theme <- ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key.height = ggplot2::unit(0.8, "cm")
    )
  if (!is.null(real)) {
    resid <- Omega_X_est - real
    clim <- max(abs(c(real, Omega_X_est)))
    rlim <- max(abs(resid))

    df_main <- to_tile(rbind(
      mat_long(real, "True Omega_X"),
      mat_long(Omega_X_est, "Estimated Omega_X")
    ))
    df_main$panel <- factor(
      df_main$panel,
      levels = c("True Omega_X", "Estimated Omega_X")
    )
    df_resid <- to_tile(mat_long(resid, "Residual"))

    g_main <- ggplot2::ggplot(
      df_main,
      ggplot2::aes(x = .data$col, y = .data$row, fill = .data$value)
    ) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~ .data$panel) +
      ggplot2::scale_fill_gradient2(
        low = "#3d52bf",
        mid = "white",
        high = "#800000",
        limits = c(-clim, clim)
      ) +
      ggplot2::coord_fixed() +
      heat_theme +
      ggplot2::labs(fill = NULL)

    g_resid <- ggplot2::ggplot(
      df_resid,
      ggplot2::aes(x = .data$col, y = .data$row, fill = .data$value)
    ) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~ .data$panel) +
      ggplot2::scale_fill_gradient2(
        low = "#3d52bf",
        mid = "white",
        high = "#800000",
        limits = c(-rlim, rlim)
      ) +
      ggplot2::coord_fixed() +
      heat_theme +
      ggplot2::labs(fill = NULL)

    idx <- which(upper.tri(real, diag = TRUE))
    lims <- range(c(real[idx], Omega_X_est[idx]))
    g_scat <- ggplot2::ggplot(
      data.frame(true = real[idx], est = Omega_X_est[idx]),
      ggplot2::aes(x = .data$true, y = .data$est)
    ) +
      ggplot2::geom_point(alpha = 0.4, size = 0.9) +
      ggplot2::geom_abline(colour = "red", linewidth = 0.8) +
      ggplot2::coord_fixed(xlim = lims, ylim = lims) +
      ggplot2::labs(
        x = "True Omega_X entries",
        y = "Estimated Omega_X entries",
        title = "Entry-wise comparison (upper triangle)"
      ) +
      ggplot2::theme_minimal(base_size = 15)

    return(list(main = g_main, residual = g_resid, scatter = g_scat))
  } else {
    clim <- max(abs(Omega_X_est))
    df_main <- to_tile(rbind(mat_long(Omega_X_est, "Estimated Omega_X")))

    g_main <- ggplot2::ggplot(
      df_main,
      ggplot2::aes(x = .data$col, y = .data$row, fill = .data$value)
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        low = "#3d52bf",
        mid = "white",
        high = "#800000",
        limits = c(-clim, clim)
      ) +
      ggplot2::coord_fixed() +
      heat_theme +
      ggplot2::labs(fill = NULL, title = "Estimated Omega_X entries")

    return(list(main = g_main))
  }
}

#' Make Predictions with a Latent Regression Model
#'
#' @param fit_blfr Returned samples from any `blfr_*()` function
#' @param X_test `n_test` times p matrix of covariates
#'
#' @returns `n_test` dimensional vector of predictions
#'
#' @export
#' @examples
#' NULL
predict_blfr <- function(fit_blfr, X_test) {
  if (!("beta_X" %in% names(fit_blfr))) {
    stop("No blfr output provided")
  }
  beta_X_hat <- rowMeans(fit_blfr$beta_X)
  yhat <- as.numeric(X_test %*% beta_X_hat)

  if ("Omega_X" %in% names(fit_blfr)) {
    Omega_X_hat <- apply(fit_blfr$Omega_X, c(1, 2), mean)
    intercept_X_hat <- mean(fit_blfr$intercept_X)
    yhat <- yhat + rowSums((X_test %*% Omega_X_hat) * X_test) + intercept_X_hat
  }

  yhat
}
