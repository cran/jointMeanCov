#' Center Each Column By Subtracting Group or Global GLS Mean
#'
#' This function takes a data matrix, an inverse row covariance
#' matrix, group indices (i.e. row indices for membership in
#' groups one and two), and a subset of column indices indicating
#' which columns should be group centered.  It returns
#' a centered data matrix.  For each group centered column,
#' the two group means are estimated using GLS; then the group
#' one mean is subtracted from entries in group one, and the group
#' two mean is subtracted from entries in group two.  For each
#' globally centered column, a single global mean is estimated
#' using GLS and subtracted from each entry in the column.
#' In addition to returning the centered data matrix, this
#' function also returns the means estimated using GLS.
#'
#' \bold{Example}
#' \preformatted{n <- 4
#' m <- 3
#' X <- matrix(1:12, nrow=n, ncol=m)
#' # Group center the first two columns, globally center
#' # the last column.
#' out <- centerDataGLSModelSelection(
#'   X, B.inv=diag(n), group.one.indices=1:2,
#'   group.two.indices=3:4,
#'   group.cen.indices=1:2)
#' # Display the centered data matrix
#' print(out$X.cen)
#' }
#'
#' @param X a data matrix.
#' @param B.inv an inverse row covariance matrix used in GLS
#' @param group.one.indices indices of observations in group one.
#' @param group.two.indices indices of observations in group two.
#' @param group.cen.indices indices of columns to be group centered
#' @return Returns a centered data matrix of the same dimensions
#' as the original data matrix.
#' \item{X.cen}{Centered data matrix.}
#' \item{group.means.gls}{Group means estimated using GLS;
#' if all columns are globally centered, then \code{NULL}.}
#' \item{global.means.gls}{Global means estimated using GLS;
#' if all columns are group centered, then \code{NULL}.}
centerDataGLSModelSelection <- function(
  X, B.inv, group.one.indices, group.two.indices,
  group.cen.indices) {

  stopifnot(nrow(X) == nrow(B.inv))
  stopifnot(length(group.cen.indices) <= ncol(X))

  n <- nrow(X)
  m <- ncol(X)
  global.cen.indices <- setdiff(1:m, group.cen.indices)

  # Design matrix for two group means
  D.group <- twoGroupDesignMatrix(
    group.one.indices, group.two.indices)
  # Design matrix for global mean
  D.global <- matrix(1, n, 1)

  # Centered data matrix
  X.cen <- matrix(NA, nrow=n, ncol=m)

  # Global centering using GLS
  if (length(global.cen.indices) > 0) {
    global.means.gls <- GLSMeans(
      X[, global.cen.indices], D.global, B.inv)
    X.cen[, global.cen.indices] <- scale(
      X[, global.cen.indices], center=global.means.gls,
      scale=FALSE)
  }

  # Group centering using GLS
  if (length(group.cen.indices) > 0) {
    group.means.gls <- GLSMeans(
      X[, group.cen.indices], D.group, B.inv)
    X.cen[group.one.indices, group.cen.indices] <- scale(
      X[group.one.indices, group.cen.indices],
      center=group.means.gls[1, ],
      scale=FALSE)
    X.cen[group.two.indices, group.cen.indices] <- scale(
      X[group.two.indices, group.cen.indices],
      center=group.means.gls[2, ],
      scale=FALSE)
  }

  out <- list()
  out$X.cen <- X.cen

  if (length(global.cen.indices) > 0) {
    out$global.means.gls <- global.means.gls
  } else {
    out$global.means.gls <- NULL
  }

  if (length(group.cen.indices) > 0) {
    out$group.means.gls <- group.means.gls
  } else {
    out$group.means.gls <- NULL
  }

  return(out)
}
