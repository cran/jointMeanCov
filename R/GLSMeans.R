#' Generalized Least Squares
#'
#' This function applies generalized least squares to estimate
#' the  unknown parameters of a linear model X = D beta + E,
#' where X has dimension n by m, D has dimension n by k, and
#' beta has dimension k by m.
#'
#' \bold{Example}
#' \preformatted{X <- matrix(1:12, nrow=4, ncol=3)
#' D <- twoGroupDesignMatrix(1:2, 3:4)
#' B.inv <- diag(4)
#' beta.hat <- GLSMeans(X, D, B.inv)
#' }
#'
#' @param X data matrix.
#' @param D design matrix.
#' @param B.inv inverse covariance matrix.
#' @return Returns the estimated parameters of the linear model,
#' a matrix of dimensions k by m, where k is the number of
#' columns of D, and m is the number of columns of X.
GLSMeans <- function(X, D, B.inv) {
  LHS <- t(D) %*% B.inv %*% D
  RHS <- t(D) %*% B.inv %*% X
  beta.hat <- solve(LHS, RHS)
  return(beta.hat)
}
