#' Design Matrix for Two-Group Mean Estimation
#'
#' This function returns the design matrix for two-group mean
#' estimation.  The first column contains indicators for
#' membership in the first group, and the second column contains
#' indicators for memebership in the second group.
#'
#' \bold{Example}
#' \preformatted{D <- twoGroupDesignMatrix(1:2, 3:5)
#' # print(D) displays the following:
#'      [,1] [,2]
#' [1,]    1    0
#' [2,]    1    0
#' [3,]    0    1
#' [4,]    0    1
#' [5,]    0    1
#' }
#'
#' @param group.one.indices indices of observations in group one.
#' @param group.two.indices indices of observations in group two.
#' @return Returns a design matrix of size n by 2, where n is
#' the sample size.
twoGroupDesignMatrix <- function(group.one.indices,
                                 group.two.indices) {

  n <- length(group.one.indices) + length(group.two.indices)
  D <- matrix(0, nrow=n, ncol=2)
  D[group.one.indices, 1] <- 1
  D[group.two.indices, 2] <- 1
  return(D)
}
