#' Estimate Row-Row Covariance Using Gemini for a Sequence of Penalties
#'
#' GeminiBPath estimates the row-row covariance, inverse
#' covariance, correlation, and inverse correlation matrices
#' using Gemini with a sequence of penalty parameters.
#' For identifiability, the covariance factors A and B are
#' scaled so that A has trace m, where m is the number of
#' columns of X, A is the column-column covariance matrix,
#' and B is the row-row covariance matrix.
#'
#' @param X Data matrix, of dimensions n by m.
#' @param rowpen.list Vector of penalty parameters, should be
#' increasing (analogous to the \code{\link[glasso]{glassopath}}
#' function of the \code{glasso} package).
#' @param penalize.diagonal Logical indicating whether to penalize the
#' off-diagonal entries of the correlation matrix.  Default is FALSE.
#' @return
#' \item{corr.B.hat}{array of estimated correlation matrices, of
#' dimension (nrow(X), nrow(X), length(rowpen.list)).}
#' \item{corr.B.hat.inv}{array of estimated inverse correlation
#' matrices, of dimension (nrow(X), nrow(X), length(rowpen.list)).}
#' \item{B.hat}{array of estimated covariance matrices, of
#' dimension (nrow(X), nrow(X), length(rowpen.list)).}
#' \item{B.hat.inv}{array of estimated inverse covariance
#' matrices, of dimension (nrow(X), nrow(X), length(rowpen.list)).}
#' @examples
#' # Generate a data matrix.
#' n1 <- 5
#' n2 <- 5
#' n <- n1 + n2
#' m <- 20
#' X <- matrix(rnorm(n * m), nrow=n, ncol=m)
#'
#' # Apply GeminiBPath for a sequence of penalty parameters.
#' rowpen.list <- sqrt(log(m) / n) * c(1, 0.5, 0.1)
#' out <- GeminiBPath(X, rowpen.list, penalize.diagonal=FALSE)
#'
#' # Display the estimated correlation matrix corresponding
#' # to penalty 0.1, rounded to two decimal places.
#' print(round(out$corr.B.hat[, , 3], 2))
#' @export
GeminiBPath <- function (X, rowpen.list,
                         penalize.diagonal=FALSE) {
  n <- nrow(X)
  m <- ncol(X)

  Gamma.hat.B <- stats::cov2cor(X %*% t(X))

  glasso.B <- glasso::glassopath(Gamma.hat.B, rowpen.list,
                         penalize.diagonal=penalize.diagonal)

  sd.row <- sqrt(rowSums(X^2))
  inv.sd.row <- 1 / sd.row

  B.hat <- array(NA, dim=c(n, n, length(rowpen.list)))
  B.hat.inv <- array(NA, dim=c(n, n, length(rowpen.list)))
  for (i in 1:length(rowpen.list)) {
    w <- glasso.B$w[, , i]
    B.hat[, , i] <- t(t(w * sd.row) * sd.row) / m
    B.hat.inv[, , i] <-  m * t(t(w * inv.sd.row) * inv.sd.row)
  }

  output <- list(corr.B.hat=glasso.B$w,
                 corr.B.hat.inv=glasso.B$wi,
                 B.hat=B.hat,
                 B.hat.inv=B.hat.inv)
}


