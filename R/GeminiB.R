#' Estimate Row-Row Covariance Structure Using Gemini
#'
#' GeminiB estimates the row-row covariance, inverse covariance,
#' correlation, and inverse correlation matrices using Gemini.
#' For identifiability, the covariance factors A and B are scaled so
#' that A has trace m, where m is the number of columns of X,
#' A is the column-column covariance matrix, and B is the row-row
#' covariance matrix.
#'
#' @param X Data matrix, of dimensions n by m.
#' @param rowpen Glasso penalty parameter.
#' @param penalize.diagonal Logical value indicating whether to penalize the
#' off-diagonal entries of the correlation matrix.  Default is FALSE.
#' @return
#' \item{corr.B.hat}{estimated correlation matrix.}
#' \item{corr.B.hat.inv}{estimated inverse correlation matrix.}
#' \item{B.hat}{estimated covariance matrix.}
#' \item{B.hat.inv}{estimated inverse covariance matrix.}
#' @examples
#' n1 <- 5
#' n2 <- 5
#' n <- n1 + n2
#' m <- 20
#' X <- matrix(rnorm(n * m), nrow=n, ncol=m)
#' rowpen <- sqrt(log(m) / n)
#' out <- GeminiB(X, rowpen, penalize.diagonal=FALSE)
#' # Display the estimated correlation matrix rounded to two
#' # decimal places.
#' print(round(out$corr.B.hat, 2))
#' @export
GeminiB <- function (X, rowpen, penalize.diagonal=FALSE) {
  n <- nrow(X)
  m <- ncol(X)

  Gamma.hat.B <- stats::cov2cor(X %*% t(X))
  glasso.B <- glasso::glasso(Gamma.hat.B, rowpen,
                             penalize.diagonal=penalize.diagonal)

  sd.row <- sqrt(rowSums(X^2))
  inv.sd.row <- 1 / sd.row

  B.hat <- t(t(glasso.B$w * sd.row) * sd.row) / m
  B.hat.inv <- m * t(t(glasso.B$wi * inv.sd.row) * inv.sd.row)

  output <- list(corr.B.hat=glasso.B$w,
                 corr.B.hat.inv=glasso.B$wi,
                 B.hat=B.hat,
                 B.hat.inv=B.hat.inv)
}
