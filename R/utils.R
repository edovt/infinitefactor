`%||%` <- function(x, y) if (is.null(x)) y else x

stick_break <- function(v) {
  v * c(1, cumprod(1 - v[-length(v)]))
}

mat_long <- function(mat, label) {
  p <- nrow(mat)
  data.frame(
    row = rep(seq_len(p), p),
    col = rep(seq_len(p), each = p),
    value = as.vector(mat),
    panel = label
  )
}
