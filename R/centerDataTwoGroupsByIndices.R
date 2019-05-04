#' Center Each Column by Subtracting Group Means
#'
#' This function takes a data matrix and returns a centered
#' data matrix.  For each column, centering is performed by
#' subtracting the corresponding group mean from each entry
#' (i.e. for entries in group one, the group one mean is
#' subtracted, and for entries in group two, the group two
#' mean is subtracted).
#'
#' \bold{Example}
#' \preformatted{X <- matrix(1:12, nrow=4, ncol=3)
#' X.cen <- centerDataTwoGroupsByIndices(
#'   X, group.one.indices=1:2, group.two.indices=3:4)
#' }
#'
#' @param X a data matrix.
#' @param group.one.indices indices of observations in group one.
#' @param group.two.indices indices of observations in group two.
#' @return Returns a centered data matrix of the same dimensions
#' as the original data matrix.
centerDataTwoGroupsByIndices <- function(
  X, group.one.indices, group.two.indices) {

  n1 <- length(group.one.indices)
  n2 <- length(group.two.indices)
  if ((n1 == 1) | (n2 == 1)) {
    warning("current implementation does not work for n1 = 1 or n2 = 1")
  }
  n <- n1 + n2
  stopifnot(n == nrow(X))
  m <- ncol(X)

  X.cenColsTwoGroups <- matrix(NA, nrow=n, ncol=m)
  X.cenColsTwoGroups[group.one.indices, ] <- scale(X[group.one.indices, ],
                                                   scale=FALSE)
  X.cenColsTwoGroups[group.two.indices, ] <- scale(X[group.two.indices, ],
                                                   scale=FALSE)
  return(X.cenColsTwoGroups)
}
