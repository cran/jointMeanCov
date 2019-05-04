#' Estimate Mean and Row-Row Correlation Matrix Using Group Centering
#'
#' This function implements Algorithm 1 from Hornstein, Fan,
#' Shedden, and Zhou (2018), doi: 10.1080/01621459.2018.1429275.
#' Given an n by m data matrix, with a vector of indices denoting
#' group membership,this function estimates the row-row inverse
#' covariance matrix after a preliminary group centering step,
#' then uses the estimated inverse covariance estimate to perform
#' GLS mean estimation. The function also returns test statistics
#' comparing the group means for each column, with standard
#' errors accounting for row-row correlation.
#'
#' @param X Data matrix.
#' @param group.one.indices Vector of indices denoting rows in group one.
#' @param rowpen Glasso penalty for estimating B, the row correlation
#' matrix.
#' @param B.inv Optional row-row covariance matrix to be used in GLS.
#' If this argument is passed, then it is used instead of estimating
#' the inverse row-row covariance.
#' @return
#' \item{B.hat.inv}{Estimated row-row inverse covariance matrix.  For identifiability,
#' A and B are scaled so that A has trace m, where m is the
#' number of columns of X.}
#' \item{corr.B.hat.inv}{Estimated row-row inverse correlation matrix.}
#' \item{gls.group.means}{Matrix with two rows and m columns, where m is
#' the number of columns of X.  Entry (i, j) contains the estimated mean
#' of the jth column for an individual in group i, with i = 1,2, and
#' j = 1, ..., m.}
#' \item{gamma.hat}{Estimated group mean differences.}
#' \item{test.stats}{Vector of test statistics of length m.}
#' \item{p.vals}{Vector of two-sided p-values, calculated using the
#' standard normal distribution.}
#' \item{p.vals.adjusted}{Vector of p-values, adjusted using the
#' Benjamini-Hochberg fdr adjustment.}
#' @examples
#' # Define sample sizes
#' n1 <- 5
#' n2 <- 5
#' n <- n1 + n2
#' m <- 200
#'
#' # Generate data with row and column covariance
#' # matrices each autorogressive of order one with
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
#' # Apply Algorithm 1 (jointMeanCovGroupCen) and
#' # plot the test statistics.
#' rowpen <- sqrt(log(m) / n)
#' out <- jointMeanCovGroupCen(X, group.one.indices, rowpen)
#' plot(out)
#' summary(out)
#' @export
jointMeanCovGroupCen <- function(X, group.one.indices, rowpen, B.inv=NULL) {
  n <- nrow(X)
  m <- ncol(X)
  group.two.indices <- setdiff(1:n, group.one.indices)

  if (is.null(B.inv)) {
    X.groupCen <- centerDataTwoGroupsByIndices(
      X, group.one.indices, group.two.indices)
    gem.B <- GeminiB(
      X.groupCen, rowpen)
    B.inv <- gem.B$B.hat.inv
    corr.B.inv <- gem.B$corr.B.hat.inv
  } else {
    corr.B.inv <- stats::cov2cor(B.inv)
  }

  D <- twoGroupDesignMatrix(group.one.indices, group.two.indices)

  gls.group.means.groupCen <- GLSMeans(
    X, D, B.inv)

  gammaHat.BiHat.groupCen <- (
    gls.group.means.groupCen[1, ] -
      gls.group.means.groupCen[2, ])

  delta <- c(1, -1)
  est.design.effects.groupCen <- as.numeric(sqrt(
    t(delta) %*% solve(t(D) %*% B.inv %*% D) %*% delta))

  test.stats <- gammaHat.BiHat.groupCen / est.design.effects.groupCen
  p.vals <- 2 * stats::pnorm(q=abs(test.stats), lower.tail=FALSE)

  out <- list()
  out$B.hat.inv <- B.inv
  out$corr.B.hat.inv <- corr.B.inv
  out$gls.group.means <- gls.group.means.groupCen
  out$gamma.hat <- gammaHat.BiHat.groupCen

  out$test.stats <- test.stats
  out$p.vals <- p.vals
  out$p.vals.adjusted <- stats::p.adjust(p.vals, method = "BH")

  class(out) <- "jointMeanCov"

  return(out)
}
