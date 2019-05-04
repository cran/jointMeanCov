#' Estimate Mean and Correlation Structure Using Stability Selection
#'
#' Given a data matrix, this function performs stability
#' selection as described in the section "Stability of Gene
#' Sets" in the paper Joint mean and covariance estimation with
#' unreplicated matrix-variate data (2018),
#' M. Hornstein, R. Fan, K. Shedden, and S. Zhou;
#' Journal of the American Statistical Association.
#'
#' Let \code{m[i]} denote the number of group-centered genes on
#' the ith iteration of stability selection (where \code{m[i]}
#' is a decreasing sequence).
#' Estimated group means are initialized using unweighted
#' sample means.  Then, for each iteration of stability
#' selection:
#' 1. The top \code{m[i]} genes are selected for group centering
#' by ranking the estimated group differences from the previous
#' iteration.
#' 2. Group means and global means are estimated using
#' GLS, using the inverse row covariance matrix from the
#' previous iteration.  The centered data matrix is then
#' used as input to Gemini to estimate the inverse row covariance
#' matrix B.hat.inv.
#' 3. Group means are estimated using GLS with B.hat.inv.
#'
#' @param X Data matrix of size n by m.
#' @param group.one.indices Vector of indices denoting rows in group one.
#' @param rowpen Glasso penalty for estimating B, the row-row correlation
#' matrix.
#' @param n.genes.to.group.center Vector specifying the number of
#' genes to group center on each iteration of the stability
#' selection algorithm.  The length of this vector is equal to
#' the number of iterations of stability selection.  If this
#' argument is not provided, the default is a decreasing
#' sequence starting with m, followed
#' @return
#' \item{n.genes.to.group.center}{Number of group centered genes
#' on each iteration of stability selection.}
#' \item{betaHat.init}{Matrix of size 2 by m, containing
#' sample means for each group.  Row 1 contains sample means for
#' group one, and row 2 contains sample means for group two.}
#' \item{gammaHat.init}{Vector of length m containing differences
#' in sample means.}
#' \item{B.inv.list}{List of estimated row-row inverse covariance
#' matrices, where \code{B.inv.list[[i]]} corresponds
#' to the estimate from the ith iteration of the algorithm, in
#' which the number of group-centered genes is
#' \code{n.genes.to.group.center[i]}.  For identifiability,
#' A and B are scaled so that A has trace m.}
#' \item{corr.B.inv.list}{List of inverse correlation matrices
#' corresponding to the inverse covariance matrices
#' \code{B.inv.list}.}
#' \item{betaHat}{List of matrices of size 2 by m, where m is
#' the number of columns of X.  For each matrix, entry (i, j) contains the
#' estimated mean of the jth column for an individual in
#' group i, with i = 1,2, and j = 1, ..., m.  The matrix
#' \code{betaHat[[i]]} contains the estimates for the ith
#' iteration of stability selection.}
#' \item{gamma.hat}{List of vectors of estimated group mean
#' differences.  The vector \code{gammaHat[[i]]} contains
#' estimates for the ith iteration of stability selection.}
#' \item{design.effecs}{Vector containing the estimated design
#' effect for each iteration of stability selection.}
#' \item{gls.test.stats}{List of vectors of test statistics for
#' each iteration of stability selection.}
#' \item{p.vals}{List of vectors of two-sided p-values, calculated using the
#' standard normal distribution.}
#' \item{p.vals.adjusted}{List of vectors of p-values, adjusted using the
#' Benjamini-Hochberg fdr adjustment.}
#' @examples
#' # Generate matrix-variate data.
#' n1 <- 5
#' n2 <- 5
#' n <- n1 + n2
#' group.one.indices <- 1:5
#' group.two.indices <- 6:10
#' m <- 20
#' M <- matrix(0, nrow=n, ncol=m)
#' # In this example, the first three variables have nonzero
#' # mean differences.
#' M[1:n1, 1:3] <- 3
#' M[(n1 + 1):n2, 1:3] <- -3
#' X <- matrix(rnorm(n * m), nrow=n, ncol=m) + M
#'
#' # Apply the stability algorithm.
#' rowpen <- sqrt(log(m) / n)
#' n.genes.to.group.center <- c(10, 5, 2)
#' out <- jointMeanCovStability(
#'  X, group.one.indices, rowpen, c(2e3, n.genes.to.group.center))
#'
#' # Make quantile plots of the test statistics for each
#' # iteration of the stability algorithm.
#' opar <- par(mfrow=c(2, 2), pty="s")
#' qqnorm(out$gammaHat.init,
#'   main=sprintf("%d genes group centered", m))
#' abline(a=0, b=1)
#' qqnorm(out$gammaHat[[1]],
#'   main=sprintf("%d genes group centered",
#'    n.genes.to.group.center[1]))
#' abline(a=0, b=1)
#' qqnorm(out$gammaHat[[2]],
#'   main=sprintf("%d genes group centered",
#'    n.genes.to.group.center[2]))
#' abline(a=0, b=1)
#' qqnorm(out$gammaHat[[3]],
#'   main=sprintf("%d genes group centered",
#'    n.genes.to.group.center[3]))
#' abline(a=0, b=1)
#' par(opar)
#' @export
jointMeanCovStability <- function(
  X, group.one.indices, rowpen, n.genes.to.group.center=NULL) {

  # Sample sizes and group indices
  n <- nrow(X)
  m <- ncol(X)
  group.two.indices <- setdiff(1:n, group.one.indices)

  # Design matrix for two group means
  D.group <- twoGroupDesignMatrix(group.one.indices, group.two.indices)
  # Design matrix for global mean
  D.global <- matrix(1, n, 1)

  B.inv.list <- list()
  corr.B.inv.list <- list()
  betaHat <- list()
  gammaHat <- list()
  design.effects <- rep(NA, length(n.genes.to.group.center))
  gls.test.stats <- list()

  group.cen.indices <- list()
  global.cen.indices <- list()
  gls.means.forCentering <- list()

  delta <- c(1, -1)

  betaHat.init <- GLSMeans(X, D.group, diag(n))
  gammaHat.init <- betaHat.init[1, ] - betaHat.init[2, ]

  if (is.null(n.genes.to.group.center)) {
    pow <- floor(log2(m))
    n.genes.to.group.center <- c(m, round(m / 2^(1:pow)))
  } else {
    if (any(diff(n.genes.to.group.center) > 0)) {
      stop("n.genes.to.group.center must be a decreasing sequence")
    }
  }

  for (i in 1:length(n.genes.to.group.center)) {
    # Rank the columns by group mean differences from previous
    # iteration; select a subset of genes to group center.
    if (i == 1) {
      gammaHat.prev <- gammaHat.init
    } else {
      gammaHat.prev <- gammaHat[[i - 1]]
    }
    ranks <- order(abs(gammaHat.prev), decreasing=TRUE)
    group.cen.indices[[i]] <- sort(
      ranks[1:n.genes.to.group.center[i]])
    global.cen.indices[[i]] <- setdiff(
      1:m, group.cen.indices[[i]])

    # Center the data by subtracting GLS group means or global
    # means from the group-centered and globally centered
    # columns, respectively.
    if (i == 1) {
      B.inv.prev <- diag(n)
    } else {
      B.inv.prev <- B.inv.list[[i - 1]]
    }
    cen <- centerDataGLSModelSelection(
      X, B.inv.prev, group.one.indices, group.two.indices,
      group.cen.indices[[i]])
    X.modSelCen <- cen$X.cen

    # Apply Gemini to the centered data matrix.
    gem.B <- GeminiB(X.modSelCen, rowpen)
    B.inv.list[[i]] <- gem.B$B.hat.inv
    corr.B.inv.list[[i]] <- gem.B$corr.B.hat.inv

    # Estimate mean differences with GLS
    betaHat[[i]] <- GLSMeans(X, D.group, gem.B$B.hat.inv)
    gammaHat[[i]] <- betaHat[[i]][1, ] - betaHat[[i]][2, ]

    # Calculate test statistics
    design.effects[i] <- (
      t(delta) %*% solve(t(D.group) %*% gem.B$B.hat.inv
                         %*% D.group) %*% delta)
    gls.test.stats[[i]] <- gammaHat[[i]] / sqrt(design.effects[i])
  }

  p.vals.gls <- lapply(
    gls.test.stats,
    function(x) 2 * stats::pnorm(q=abs(x), lower.tail=FALSE))

  p.vals.gls.adjusted <- lapply(
    X=p.vals.gls, FUN=stats::p.adjust, method="BH")

  out <- list()
  out$n.genes.to.group.center <- n.genes.to.group.center
  out$betaHat.init <- betaHat.init
  out$gammaHat.init <- gammaHat.init
  out$B.inv.list <- B.inv.list
  out$corr.B.inv.list <- corr.B.inv.list
  out$betaHat <- betaHat
  out$gammaHat <- gammaHat
  out$design.effects <- design.effects
  out$gls.test.stats <- gls.test.stats
  out$p.vals.gls <- p.vals.gls
  out$p.vals.gls.adjusted <- p.vals.gls.adjusted

  return(out)
}




