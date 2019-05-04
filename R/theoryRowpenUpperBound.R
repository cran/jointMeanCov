#' Penalty Parameter for Covariance Estimation Based on Theory
#'
#' This function returns a theoretically-guided choice of
#' the glasso penalty parameter, based on both the row and
#' column covariance matrices.
#'
#' @param A column covariance matrix.
#' @param B row covariance matrix.
#' @param n1 sample size of group one.
#' @param n2 sample size of group two.
#' @return Returns a theoretically guided choice of the
#' glasso penalty parameter.
#' @references Joint mean and covariance estimation with unreplicated matrix-variate data
#' Michael Hornstein, Roger Fan, Kerby Shedden, Shuheng Zhou
#' (2018).  Joint mean and covariance estimation with
#' unreplicated matrix-variate data.  Journal of the American
#' Statistical Association
#' @examples
#' # Define sample sizes
#' n1 <- 10
#' n2 <- 10
#' n <- n1 + n2
#' m <- 2e3
#' # Column covariance matrix (autoregressive of order 1)
#' A <- outer(1:n, 1:n, function(x, y) 0.2^abs(x - y))
#' # Row covariance matrix (autoregressive of order 1)
#' B <- outer(1:n, 1:n, function(x, y) 0.8^abs(x - y))
#' # Calculate theoretically guided Gemini penalty.
#' rowpen <- theoryRowpenUpperBound(A, B, n1, n2)
#' print(rowpen)
#' @export
theoryRowpenUpperBound <- function(A, B, n1, n2) {
  n <- n1 + n2
  stopifnot(n == nrow(B))
  m <- nrow(A)
  bias.term <- 3 * max(rowSums((abs(B)))) / n
  var.term <- (sqrt(sum(A^2)) / sum(diag(A))) * sqrt(log(m))
  theoretical.penalty <- var.term + bias.term
}


