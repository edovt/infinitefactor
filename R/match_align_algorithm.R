#' Resolve rotational ambiguity with the MatchAlign algortithm
#'
#' Performs the varimax rotation on the factor loadings samples and
#' column-based matching to resolve resultant sign and label switching.
#' Rotates the factors along with the loadings to induce identifiability jointly.
#'
#' @param fit Returned samples from any `bfa_*()` and `blfr_*()` functions with fixed number of factors.
#' @param orth_procedure Orthogonalization procedure to use. Currently only supports VARIMAX.
#'
#' @returns List of rotationally aligned factor loadings and factors samples
#'
#' @export
#' @examples
#' set.seed(219)
#'
#' # 0. Data simulation as in Bhattacharya & Dunson (2011)
#' sim <- sim_data_bfa(n = 200, p = 20, k = 3)
#'
#' # 1. Adaptive number of factors, starting at floor(3*log(p)) = 8
#' fit_cusp <- bfa_cusp(sim$X)
#' k_fit <- summary_k(fit_cusp, real = sim$k_true) # returns median too
#'
#' # 2. Fixed number of factors -> MatchAlign
#' fit_cusp_fixed <- bfa_mgps(sim$X, adapt = FALSE, k_init = k_fit)
#' aligned <- match_align(fit_cusp_fixed)
#' plot_match_align(aligned, "Lambda", "mean")
#' plot_match_align(aligned, "Eta", "mean")
match_align <- function(fit, orth_procedure = "varimax") {
  if (!("Lambda" %in% names(fit))) {
    stop("Model fitted with changing number of factors")
  }
  # Coerce 3D arrays to lists (previous implementation)
  Lambda <- lapply(seq_len(dim(fit$Lambda)[3L]), function(s) fit$Lambda[,, s])
  Eta <- lapply(seq_len(dim(fit$Eta)[3L]), function(s) fit$Eta[,, s])

  # 1. Apply orthogonalization procedure to each sample
  if (orth_procedure == "varimax") {
    orth_output <- lapply(Lambda, stats::varimax, normalize = FALSE)
  } else {
    # orth_output <- lapply(Lambda, get(orth_procedure))
    stop("Only VARIMAX rotation can be used.")
  }

  # 2. Choose the pivot
  loads <- lapply(orth_output, `[[`, 1)
  rots <- lapply(orth_output, `[[`, 2)
  rotfact <- mapply(`%*%`, Eta, rots, SIMPLIFY = FALSE)
  norms <- sapply(loads, norm, "2")
  piv <- loads[order(norms)][[round(length(loads) / 2)]]

  # 3. Solve column-label / sign ambiguities
  matches <- lapply(loads, msfOUT, piv)
  lambda_out <- mapply(aplr, loads, matches, SIMPLIFY = FALSE)
  eta_out <- mapply(aplr, rotfact, matches, SIMPLIFY = FALSE)

  list(Lambda = lambda_out, Eta = eta_out)
}

#' Plots of results of MatchAlign
#'
#' Plots the summary matrix after performing a rotation of the factors and
#' loadings to induce joint identifiability.
#'
#' @param match_align_out Output of `match_align()`.
#' @param type Which samples to plot? One of `c("Lambda", "Eta")`.
#' @param summary Which summary metric to use? One of `c("mean", "median")`.
#' @param color Color scheme. One of `c("green", "red", "wes")`
#' @param title Optional plot title
#' @param X_names Optional character vector of covariate names, of length equal
#'   to the number of variables `p`. When supplied (and the plotted matrix has
#'   `p` rows, i.e. `type = "Lambda"`), these are used to label the y axis.
#'
#' @returns ggplot object
#'
#' @export
#' @examples
#' set.seed(219)
#'
#' # 0. Data simulation as in Bhattacharya & Dunson (2011)
#' sim <- sim_data_bfa(n = 200, p = 20, k = 3)
#'
#' # 1. Adaptive number of factors, starting at floor(3*log(p)) = 8
#' fit_cusp <- bfa_cusp(sim$X)
#' k_fit <- summary_k(fit_cusp, real = sim$k_true) # returns median too
#'
#' # 2. Fixed number of factors -> MatchAlign
#' fit_cusp_fixed <- bfa_mgps(sim$X, adapt = FALSE, k_init = k_fit)
#' aligned <- match_align(fit_cusp_fixed)
#' plot_match_align(aligned, "Lambda", "mean")
#' plot_match_align(aligned, "Eta", "mean")
plot_match_align <- function(
  match_align_out,
  type = c("Lambda", "Eta"),
  summary = c("mean", "median"),
  color = "green",
  title = NULL,
  X_names = NULL
) {
  type <- match.arg(type)
  summary <- match.arg(summary)

  samples <- match_align_out[[type]]
  arr <- simplify2array(samples)
  reducer <- if (summary == "mean") base::mean else stats::median
  mat <- apply(arr, c(1L, 2L), reducer)
  n_row <- nrow(mat)
  n_fac <- ncol(mat)
  if (!is.null(X_names) && length(X_names) != n_row) {
    stop(sprintf(
      "X_names has length %d but the %s matrix has %d rows.",
      length(X_names), type, n_row
    ))
  }

  mat <- apply(mat, 2, rev)
  longmat <- reshape2::melt(mat)
  Var1 <- Var2 <- value <- NULL # nolint
  p <- ggplot2::ggplot(longmat, ggplot2::aes(x = Var2, y = Var1)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value), colour = "grey20")

  if (color == "green") {
    p <- p +
      ggplot2::scale_fill_gradient2(
        low = "#3d52bf",
        high = "#33961b",
        mid = "white"
      )
  } else if (color == "red") {
    p <- p +
      ggplot2::scale_fill_gradient2(
        low = "#191970",
        high = "#800000",
        mid = "white"
      )
  } else if (color == "wes") {
    p <- p +
      ggplot2::scale_fill_gradient2(
        low = "#046C9A",
        high = "#D69C4E",
        mid = "white"
      )
  }

  p <- p +
    ggplot2::scale_x_continuous(
      breaks = seq_len(n_fac),
      labels = paste("Factor", seq_len(n_fac))
    )
  if (!is.null(X_names)) {
    p <- p +
      ggplot2::scale_y_continuous(
        breaks = seq_len(n_row),
        labels = rev(X_names)
      )
  }

  p <- p +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = if (is.null(X_names)) {
        ggplot2::element_blank()
      } else {
        ggplot2::element_text()
      },
      legend.title = ggplot2::element_text(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(fill = " ")
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }
  p
}
