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
      as.integer(stats::quantile(fit$k, 0.975))
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
#' @param real Real matrix if known, `NULL` otherwise. When `covariance = TRUE`
#'   this is the real covariance matrix. When `covariance = FALSE` it is the real
#'   correlation matrix.
#' @param X_names Optional character vector of covariate names, of length equal
#'   to the number of variables in `X`. Used to label axes of plots.
#' @param covariance If `TRUE` (default), plots the posterior covariance matrix.
#'   If `FALSE`, plots the implied correlation matrix instead.
#'
#' @returns Named list of plots. If `real=NULL`, `plots$main` is the plot of the
#' mean posterior covariance (or correlation) of X. Otherwise, three plots:
#' * `main`: compares mean posterior estimate with the real matrix
#' * `residual`: plot of residuals between mean posterior and real matrix
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
plot_Sigma_X <- function(fit, real = NULL, X_names = NULL, covariance = TRUE) {
  Sigma_X_est <- apply(fit$Sigma_X, c(1, 2), mean)
  if (!covariance) {
    Sigma_X_est <- stats::cov2cor(Sigma_X_est)
  }
  mat_label <- if (covariance) "Sigma_X" else "Corr_X"
  p <- nrow(Sigma_X_est)
  if (!is.null(X_names) && length(X_names) != p) {
    stop(sprintf(
      "X_names has length %d but Sigma_X has %d rows/columns.",
      length(X_names),
      p
    ))
  }
  to_tile <- function(df) {
    df$row <- p - df$row + 1
    df
  }
  add_axis_names <- function(g) {
    if (is.null(X_names)) {
      return(g)
    }
    g +
      ggplot2::scale_x_continuous(breaks = seq_len(p), labels = X_names) +
      ggplot2::scale_y_continuous(breaks = seq_len(p), labels = rev(X_names)) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = ggplot2::element_text()
      )
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
      mat_long(real, sprintf("True %s", mat_label)),
      mat_long(Sigma_X_est, sprintf("Estimated %s", mat_label))
    ))
    df_main$panel <- factor(
      df_main$panel,
      levels = sprintf(c("True %s", "Estimated %s"), mat_label)
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
        x = sprintf("True %s entries", mat_label),
        y = sprintf("Estimated %s entries", mat_label),
        title = "Entry-wise comparison (upper triangle)"
      ) +
      ggplot2::theme_minimal(base_size = 15)

    return(list(
      main = add_axis_names(g_main),
      residual = add_axis_names(g_resid),
      scatter = g_scat
    ))
  } else {
    clim <- max(abs(Sigma_X_est))
    df_main <- to_tile(rbind(
      mat_long(Sigma_X_est, sprintf("Estimated %s", mat_label))
    ))

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
      ggplot2::labs(
        fill = NULL,
        title = if (covariance) {
          "Covariance estimator"
        } else {
          "Correlation estimator"
        }
      )

    return(list(main = add_axis_names(g_main)))
  }
}

#' Summary Plots of Induced Main Effects
#'
#' @param fit_blfr Returned samples from any `blfr_*()` function
#' @param real_beta Real vector of main effects, if known
#' @param real_intercept Real intercept, if known
#' @param X_names Optional character vector of covariate names, of length equal
#'   to the number of predictors `p`. If supplied, these label the x axis of the
#'   caterpillar plot (prefixed with "Intercept" when the fit has an intercept).
#'
#' @returns Named list of plots
#'
#' @export
#' @examples
#' NULL
plot_effects_X <- function(
  fit_blfr,
  real_beta = NULL,
  real_intercept = NULL,
  X_names = NULL
) {
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

  n_pred <- nrow(fit_blfr$beta_X)
  if (!is.null(X_names) && length(X_names) != n_pred) {
    stop(sprintf(
      "X_names has length %d but there are %d predictors.",
      length(X_names),
      n_pred
    ))
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

  add_x_names <- function(g) {
    if (is.null(X_names)) {
      return(g)
    }
    x_labels <- if (intercept) c("Intercept", X_names) else X_names
    g +
      ggplot2::scale_x_continuous(breaks = df$idx, labels = x_labels) +
      ggplot2::xlab(NULL) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
      )
  }

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

    return(list(main = add_x_names(g_caterpillar)))
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

    return(list(main = add_x_names(g_caterpillar), scatter = g_scat))
  }
}

#' Summary Plots of Induced Interaction Effects
#'
#' @param fit_blfr Returned samples from any `blfr_*()` function with interactions
#' @param real Real matrix of induced interaction effects, if known
#' @param X_names Optional character vector of covariate names, of length equal
#'   to the number of variables in `X`. If supplied, these are used to label both
#'   the x and y axes of the interaction-effect heatmaps.
#' @param mask_ci If `TRUE` (default), entries whose `level` credible interval
#'   contains 0 are faded out in the estimated-effect heatmap, so that only the
#'   credibly non-zero interactions stand out.
#' @param level Credible-interval level used for the mask. Defaults to `0.95`.
#'
#' @returns Named list of plots
#'
#' @export
#' @examples
#' NULL
plot_Omega_X <- function(
  fit_blfr,
  real = NULL,
  X_names = NULL,
  mask_ci = TRUE,
  level = 0.95
) {
  if (!("Omega_X" %in% names(fit_blfr))) {
    stop("Model not fitted with interactions")
  }
  Omega_X_est <- apply(fit_blfr$Omega_X, c(1, 2), mean)
  p <- nrow(Omega_X_est)
  if (!is.null(X_names) && length(X_names) != p) {
    stop(sprintf(
      "X_names has length %d but Omega_X has %d rows/columns.",
      length(X_names),
      p
    ))
  }
  # Per-entry credible interval mask: entries whose CI contains 0 are made
  # transparent so that only the credibly non-zero interactions are coloured.
  if (mask_ci) {
    probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
    Omega_X_ci <- apply(fit_blfr$Omega_X, c(1, 2), stats::quantile, probs = probs)
    signif_mat <- Omega_X_ci[1, , ] > 0 | Omega_X_ci[2, , ] < 0
    mask_vec <- !as.vector(signif_mat)
  } else {
    mask_vec <- rep(FALSE, p * p)
  }
  to_tile <- function(df) {
    df$row <- p - df$row + 1
    df
  }
  add_axis_names <- function(g) {
    if (is.null(X_names)) {
      return(g)
    }
    g +
      ggplot2::scale_x_continuous(breaks = seq_len(p), labels = X_names) +
      ggplot2::scale_y_continuous(breaks = seq_len(p), labels = rev(X_names)) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = ggplot2::element_text()
      )
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

    df_true <- mat_long(real, "True Omega_X")
    df_est <- mat_long(Omega_X_est, "Estimated Omega_X")
    df_est$value[mask_vec] <- NA
    df_main <- to_tile(rbind(df_true, df_est))
    df_main$panel <- factor(
      df_main$panel,
      levels = c("True Omega_X", "Estimated Omega_X")
    )
    df_resid <- to_tile(mat_long(resid, "Residual"))

    g_main <- ggplot2::ggplot(
      df_main,
      ggplot2::aes(x = .data$col, y = .data$row, fill = .data$value)
    ) +
      ggplot2::geom_tile(colour = "grey80", linewidth = 0.2) +
      ggplot2::facet_wrap(~ .data$panel) +
      ggplot2::scale_fill_gradient2(
        low = "#3d52bf",
        mid = "white",
        high = "#800000",
        limits = c(-clim, clim),
        na.value = "transparent"
      ) +
      ggplot2::coord_fixed() +
      heat_theme +
      ggplot2::labs(
        fill = NULL,
        caption = if (mask_ci) {
          sprintf("Blank: %g%% CI contains 0", 100 * level)
        } else {
          NULL
        }
      )

    g_resid <- ggplot2::ggplot(
      df_resid,
      ggplot2::aes(x = .data$col, y = .data$row, fill = .data$value)
    ) +
      ggplot2::geom_tile(colour = "grey80", linewidth = 0.2) +
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

    return(list(
      main = add_axis_names(g_main),
      residual = add_axis_names(g_resid),
      scatter = g_scat
    ))
  } else {
    clim <- max(abs(Omega_X_est))
    df_main <- mat_long(Omega_X_est, "Estimated Omega_X")
    df_main$value[mask_vec] <- NA
    df_main <- to_tile(df_main)

    g_main <- ggplot2::ggplot(
      df_main,
      ggplot2::aes(x = .data$col, y = .data$row, fill = .data$value)
    ) +
      ggplot2::geom_tile(colour = "grey80", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(
        low = "#3d52bf",
        mid = "white",
        high = "#800000",
        limits = c(-clim, clim),
        na.value = "transparent"
      ) +
      ggplot2::coord_fixed() +
      heat_theme +
      ggplot2::labs(
        fill = NULL,
        title = "Estimated Omega_X entries",
        caption = if (mask_ci) {
          sprintf("Blank: %g%% CI contains 0", 100 * level)
        } else {
          NULL
        }
      )

    return(list(main = add_axis_names(g_main)))
  }
}

#' Summary Plots of Covariate Effects
#' @param fit_blfr Returned samples from any `blfr_*()` function with covariates
#' @param real Real vector of covariate effects, if known
#' @param Z_names Optional character vector of covariate names, of length equal
#'   to the number of covariates `q`.
#'
#' @returns Named list of plots
#'
#' @export
#' @examples
#' NULL
plot_effects_Z <- function(fit_blfr, real = NULL, Z_names = NULL) {
  if (!("alpha" %in% names(fit_blfr))) {
    stop("Model not fitted with covariates")
  }
  alpha_est <- apply(fit_blfr$alpha, 1, mean)
  alpha_q05 <- apply(fit_blfr$alpha, 1, stats::quantile, probs = 0.025)
  alpha_q95 <- apply(fit_blfr$alpha, 1, stats::quantile, probs = 0.975)
  q <- length(alpha_est)
  if (!is.null(Z_names) && length(Z_names) != q) {
    stop(sprintf(
      "Z_names has length %d but alpha has %d elements.",
      length(Z_names),
      q
    ))
  }

  df <- data.frame(
    idx = seq_along(alpha_est),
    est = alpha_est,
    lo = alpha_q05,
    hi = alpha_q95
  )
  add_z_names <- function(g) {
    if (is.null(Z_names)) {
      return(g)
    }
    g +
      ggplot2::scale_x_continuous(breaks = df$idx, labels = Z_names) +
      ggplot2::xlab(NULL) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 20
        )
      )
  }

  if (is.null(real)) {
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
        x = "Covariate index",
        y = expression(alpha[Z]),
        title = "Covariates regression coefficients (95% CI)"
      ) +
      ggplot2::theme_minimal(base_size = 13)

    return(list(main = add_z_names(g_caterpillar)))
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
        x = "Covariate index",
        y = expression(alpha[Z]),
        title = "Covariates regression coefficients (95% CI)",
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
        x = expression("True " * alpha[Z]),
        y = expression("Estimated " * alpha[Z])
      ) +
      ggplot2::theme_minimal(base_size = 13)

    return(list(main = add_z_names(g_caterpillar), scatter = g_scat))
  }
}

#' Summary Plots of Induced Covariate and Predictor effects.
#' @param fit_blfr Returned samples from any `blfr_*()` function with covariates
#' @param real Real matrix of induced covariate and predictor interaction effects,
#'   if known
#' @param X_names Optional character vector of predictor names, of length equal
#'   to the number of variables in `X`.
#' @param Z_names Optional character vector of covariate names, of length equal
#'   to the number of covariates `q`
#'
#' @returns Named list of plots
#'
#' @export
#' @examples
#' NULL
plot_Delta_X <- function(
  fit_blfr,
  real = NULL,
  X_names = NULL,
  Z_names = NULL
) {
  if (!("Delta_X" %in% names(fit_blfr))) {
    stop("Model not fitted with covariate-predictor interactions")
  }
  # I take the transpose since in general q < p
  Delta_X_est <- t(apply(fit_blfr$Delta_X, c(1, 2), mean)) # q x p
  q <- nrow(Delta_X_est)
  p <- ncol(Delta_X_est)
  if (!is.null(X_names) && length(X_names) != p) {
    stop(sprintf(
      "X_names has length %d but Delta_X has %d rows.",
      length(X_names),
      p
    ))
  }
  if (!is.null(Z_names) && length(Z_names) != q) {
    stop(sprintf(
      "Z_names has length %d but Delta_X has %d columns.",
      length(Z_names),
      q
    ))
  }

  # Delta_X_est is q x p: row = covariate (y axis), col = predictor (x axis)
  to_long <- function(mat, label) {
    nr <- nrow(mat)
    nc <- ncol(mat)
    data.frame(
      row = rep(seq_len(nr), nc),
      col = rep(seq_len(nc), each = nr),
      value = as.vector(mat),
      panel = label
    )
  }
  to_tile <- function(df) {
    df$row <- q - df$row + 1
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
  add_axis_names <- function(g) {
    if (!is.null(X_names)) {
      g <- g +
        ggplot2::scale_x_continuous(breaks = seq_len(p), labels = X_names) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
        )
    }
    if (!is.null(Z_names)) {
      g <- g +
        ggplot2::scale_y_continuous(
          breaks = seq_len(q),
          labels = rev(Z_names)
        ) +
        ggplot2::theme(axis.text.y = ggplot2::element_text())
    }
    g
  }

  if (!is.null(real)) {
    if (nrow(real) != p || ncol(real) != q) {
      stop(sprintf(
        "real has dimension %d x %d but Delta_X is %d x %d (p x q).",
        nrow(real),
        ncol(real),
        p,
        q
      ))
    }
    real <- t(real) # match q x p orientation of Delta_X_est
    resid <- Delta_X_est - real
    clim <- max(abs(c(real, Delta_X_est)))
    rlim <- max(abs(resid))

    df_main <- to_tile(rbind(
      to_long(real, "True Delta_X"),
      to_long(Delta_X_est, "Estimated Delta_X")
    ))
    df_main$panel <- factor(
      df_main$panel,
      levels = c("True Delta_X", "Estimated Delta_X")
    )
    df_resid <- to_tile(to_long(resid, "Residual"))

    g_main <- ggplot2::ggplot(
      df_main,
      ggplot2::aes(x = .data$col, y = .data$row, fill = .data$value)
    ) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~ .data$panel, ncol = 1) +
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

    lims <- range(c(real, Delta_X_est))
    g_scat <- ggplot2::ggplot(
      data.frame(true = as.vector(real), est = as.vector(Delta_X_est)),
      ggplot2::aes(x = .data$true, y = .data$est)
    ) +
      ggplot2::geom_point(alpha = 0.4, size = 0.9) +
      ggplot2::geom_abline(colour = "red", linewidth = 0.8) +
      ggplot2::coord_fixed(xlim = lims, ylim = lims) +
      ggplot2::labs(
        x = "True Delta_X entries",
        y = "Estimated Delta_X entries",
        title = "Entry-wise comparison"
      ) +
      ggplot2::theme_minimal(base_size = 15)

    return(list(
      main = add_axis_names(g_main),
      residual = add_axis_names(g_resid),
      scatter = g_scat
    ))
  } else {
    clim <- max(abs(Delta_X_est))
    df_main <- to_tile(to_long(Delta_X_est, "Estimated Delta_X"))

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
      ggplot2::labs(fill = NULL, title = "Estimated Delta_X entries")

    return(list(main = add_axis_names(g_main)))
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
