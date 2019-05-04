#' Quantile Plot of Test Statistics
#'
#' This function displays a quantile plot of test statistics,
#' based on the output of the functions
#' \code{\link{jointMeanCovGroupCen}} or \code{\link{jointMeanCovModSelCen}}.
#'
#' @param x output of \code{\link{jointMeanCovGroupCen}} or \code{\link{jointMeanCovModSelCen}}.
#' @param ... other plotting arguments passed to
#' \code{\link[stats]{qqnorm}}.
#' @examples
#' # Define sample sizes
#' n1 <- 5
#' n2 <- 5
#' n <- n1 + n2
#' m <- 200
#'
#' # Generate data with row and column covariance
#' # matrices each autorogressive of order 1 with
#' # parameter 0.2.  The mean is defined so the first
#' # three columns have true differences in group means
#' # equal to four.
#' Z <- matrix(rnorm(m * n), nrow=n, ncol=m)
#' A <- outer(1:m, 1:m, function(i, j) 0.2^abs(i - j))
#' B <- outer(1:n, 1:n, function(i, j) 0.2^abs(i - j))
#' M <- matrix(0, nrow=nrow(Z), ncol=ncol(Z))
#' group.one.indices <- 1:5
#' group.two.indices <- 6:10
#' M[group.one.indices, 1:3] <- 2
#' M[group.two.indices, 1:3] <- -2
#' X <- t(chol(B)) %*% Z %*% chol(A) + M
#'
#' # Apply Algorithm 2 (jointMeanCovModSelCen) and plot the
#' # test statistics.
#' rowpen <- sqrt(log(m) / n)
#' out <- jointMeanCovModSelCen(X, group.one.indices, rowpen)
#' plot(out)
#' @export
plot.jointMeanCov <- function(x, ...) {
  stats::qqnorm(x$test.stats,
                main = "Quantile plot of test statistics", ...)
  graphics::abline(a=0, b=1)
}
