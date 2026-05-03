`%||%` <- function(x, y) if (is.null(x)) y else x

stick_break <- function(v) {
  v * c(1, cumprod(1 - v[-length(v)]))
}
