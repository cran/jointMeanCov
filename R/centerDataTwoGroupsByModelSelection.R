#' Center Each Column Based on Model Selection
#'
#' This function takes a data matrix and returns a centered data
#' matrix.  For columns with indices in \code{within.group.indices},
#' centering is performed by subtracting the corresponding
#' group mean from each entry (i.e. for entries in group one,
#' the group one mean is subtracted, and for entries in group two,
#' the group two mean is subtracted).  For other columns, global
#' centering is performed (i.e. subtracting the column mean from
#' each entry).
#'
#' \bold{Example}
#' \preformatted{X <- matrix(1:12, nrow=4, ncol=3)
#' # Group center the first two columns, globally center
#' # the third column.
#' X.cen <- centerDataTwoGroupsByModelSelection(
#'   X, group.one.indices=1:2, group.two.indices=3:4,
#'   within.group.indices=1:2)
#' }
#'
#' @param X a data matrix.
#' @param group.one.indices indices of observations in group one.
#' @param group.two.indices indices of observations in group two.
#' @param within.group.indices indices of columns on which to
#' perform group centering.
#' @return Returns a centered data matrix of the same dimensions
#' as the original data matrix.
centerDataTwoGroupsByModelSelection <- function(
  X, group.one.indices, group.two.indices, within.group.indices) {

  n1 <- length(group.one.indices)
  n2 <- length(group.two.indices)
  if ((n1 == 1) | (n2 == 1)) {
    stop("current implementation does not work for n1 = 1 or n2 = 1")
  }
  n <- n1 + n2
  stopifnot(n == nrow(X))
  m <- ncol(X)

  if (length(within.group.indices) != 0) {
    stopifnot(max(within.group.indices) <= m)
  }

  colCen.indices <- setdiff(1:m, within.group.indices)

  X.cen <- matrix(NA, nrow=n, ncol=m)
  X.cen[, colCen.indices] <- scale(X[, colCen.indices], center=TRUE,
                                   scale=FALSE)
  X.cen[group.one.indices, within.group.indices] <- scale(
    X[group.one.indices, within.group.indices], scale=FALSE)
  X.cen[group.two.indices, within.group.indices] <- scale(
    X[group.two.indices, within.group.indices], scale=FALSE)

  return(X.cen)
}
