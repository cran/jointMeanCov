#' Penalty Parameter for Covariance Estimation Based on Theory
#'
#' This function returns a theoretically-guided choice of
#' the glasso penalty parameter, treating the column
#' correlation matrix as the identity.
#'
#' @param B row covariance matrix.
#' @param n1 sample size of group one.
#' @param n2 sample size of group two.
#' @param m number of columns of the data matrix (where the
#' data matrix is of size n by m, with n = n1 + n2).
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
#' # Row covariance matrix (autoregressive of order 1)
#' B <- outer(1:n, 1:n, function(x, y) 0.8^abs(x - y))
#' # Calculate theoretically guided Gemini penalty.
#' rowpen <- theoryRowpenUpperBoundDiagA(B, n1, n2, m)
#' print(rowpen)
#' @export
theoryRowpenUpperBoundDiagA <- function(B, n1, n2, m) {
  n <- n1 + n2
  stopifnot(n == nrow(B))
  bias.term <- 3 * max(rowSums((abs(B)))) / n
  var.term <- sqrt(log(m) / m)
  theoretical.penalty <- var.term + bias.term
}
