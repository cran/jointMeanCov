#' Estimate Mean and Row-Row Correlation Matrix Using Model Selection
#'
#' This function implements Algorithm 2 from Hornstein, Fan,
#' Shedden, and Zhou (2018), doi: 10.1080/01621459.2018.1429275.
#' Given an n by m data matrix, with a vector of indices denoting
#' group membership, this function estimates mean and covariance
#' structure as follows.  1. Run Algorithm 1
#' (\code{jointMeanCovGroupCen}).  2. Use a threshold to select
#' genes with the largest mean differences.  3. Group center
#' the genes with mean differences above the threshold, and
#' globally center the remaining genes.  Use the centered data
#' matrix to calculate a Gram matrix as input to Gemini.
#' 4. Use Gemini to estimate the inverse row covariance matrix,
#' and use the inverse row covariance matrix with GLS to
#' estimate group means.  5. Calculate test statistics comparing
#' group means for each column.
#'
#' @param X Data matrix.
#' @param group.one.indices Vector of indices denoting rows in group one.
#' @param rowpen Glasso penalty for estimating B, the row-row correlation
#' matrix.
#' @param B.inv Optional row-row covariance matrix to be used in GLS in
#' Algorithm 1 prior to model selection centering.
#' If this argument is passed, then it is used instead of estimating
#' the inverse row-row covariance.
#' @param rowpen.ModSel Optional Glasso penalty for estimating B
#' in the second step.
#' @param thresh Threshold for model selection centering.  If
#' group means for a column differ by less than the threshold,
#' the column is globally centered rather than group centered.  If
#' \code{thresh} is \code{NULL}, then the theoretically guided
#' threshold is used.
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
#' # Apply Algorithm 2 (jointMeanCovModSelCen) and
#' # plot the test statistics.
#' rowpen <- sqrt(log(m) / n)
#' out <- jointMeanCovModSelCen(X, group.one.indices, rowpen)
#' plot(out)
#' summary(out)
#' @export
jointMeanCovModSelCen <- function(
  X, group.one.indices, rowpen, B.inv=NULL,
  rowpen.ModSel=NULL,
  thresh=NULL) {

  n <- nrow(X)
  n1 <- length(group.one.indices)
  n2 <- n - n1
  m <- ncol(X)

  group.two.indices <- setdiff(1:n, group.one.indices)

  if (is.null(rowpen.ModSel)) {
    rowpen.ModSel <- rowpen
  }

  group.cen.out <- jointMeanCovGroupCen(
    X, group.one.indices, rowpen, B.inv)

  if (is.null(thresh)) {
    thresh <- theoryRowpenUpperBoundDiagA(
      B=group.cen.out$B.hat.inv, n1, n2, m)
  }

  groupCen.indices <- which(abs(group.cen.out$gamma.hat) >= thresh)

  X.modSelCen <- centerDataTwoGroupsByModelSelection(
    X, group.one.indices, group.two.indices,
    groupCen.indices)

  gem.B.modSel <- GeminiB(
    X.modSelCen, rowpen.ModSel)
  BiHat.modSelCen.useAllGenes  <- gem.B.modSel$B.hat.inv

  D <- twoGroupDesignMatrix(group.one.indices, group.two.indices)
  gls.group.means.modSelCen <- GLSMeans(
    X, D, BiHat.modSelCen.useAllGenes)
  gammaHat.BiHat.modSelCen <- (
    gls.group.means.modSelCen[1, ] -
      gls.group.means.modSelCen[2, ])

  delta <- c(1, -1)
  est.design.effects.modSelCen <- as.numeric(sqrt(
    t(delta) %*% solve(t(D) %*% BiHat.modSelCen.useAllGenes %*% D) %*% delta))

  test.stats<- gammaHat.BiHat.modSelCen / est.design.effects.modSelCen
  p.vals <- 2 * stats::pnorm(q=abs(test.stats), lower.tail=FALSE)

  out <- list()
  out$B.hat.inv <- BiHat.modSelCen.useAllGenes
  out$corr.B.hat.inv <- gem.B.modSel$corr.B.hat.inv
  out$gls.group.means <- gls.group.means.modSelCen
  out$gamma.hat <- gammaHat.BiHat.modSelCen

  out$test.stats <- test.stats
  out$p.vals <- p.vals
  out$p.vals.adjusted <- stats::p.adjust(p.vals, method = "BH")

  class(out) <- "jointMeanCov"

  return(out)
}
