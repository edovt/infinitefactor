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
      as.integer(median(fit$k)),
      as.integer(quantile(fit$k, 0.025)),
      as.integer(quantile(fit$k, 0.975)),
    ))
  } else {
    cat(sprintf(
      "Estimated k -- median: %d, 95%% CI: [%d, %d]  (true k = %d)\n",
      as.integer(median(fit$k)),
      as.integer(quantile(fit$k, 0.025)),
      as.integer(quantile(fit$k, 0.975)),
      real
    ))
  }

  median(fit$k)
}

#' Summary plots of the posterior covariance of X
#'
#' @param fit Returned samples from any `bfa_*()` or `blfr_*()` function.
#' @param real Real Sigma_X if known, `NULL` otherwise
#'
#' @returns List of plots. If `real=NULL`, `plots$main` is the plot of the
#'     mean posterior covariance of X. Otherwise, three plots `main`, `residual`
#'     and `scatter`: first compares mean posterior covariance with real
#'     covariance, second a plot of the residuals between the two, and third
#'     plot is a scatterplot of mean posterior elements vs real elements.
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
      ggplot2::facet_wrap(~panel) +
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
      ggplot2::facet_wrap(~panel) +
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
      ggplot2::aes(x = true, y = est)
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
