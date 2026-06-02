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

## Half-vectorization and inverse, taken from OpenMx
vech <- function(X) {
  X[lower.tri(X, diag = TRUE)]
}

vech2sym <- function(X) {
  if (is.matrix(X)) {
    if (nrow(X) > 1 && ncol(X) > 1) {
      stop("Input to the full vech2full must be a (1 x n) or (n x 1) matrix.")
    }
    dimension <- max(dim(X))
  } else if (is.vector(X)) {
    dimension <- length(X)
  } else {
    stop("Input to the function vech2full must be either a matrix or a vector.")
  }
  k <- sqrt(2 * dimension + 0.25) - 0.5
  ret <- matrix(0, nrow = k, ncol = k)
  if (nrow(ret) != k) {
    stop(
      "Incorrect number of elements in vector to construct a matrix from a half-vectorization."
    )
  }
  ret[lower.tri(ret, diag = TRUE)] <- as.vector(X)
  ret[upper.tri(ret)] <- t(ret)[upper.tri(ret)]
  return(ret)
}
