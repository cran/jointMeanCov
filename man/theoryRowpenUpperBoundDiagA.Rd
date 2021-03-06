% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/theoryRowpenUpperBoundDiagA.R
\name{theoryRowpenUpperBoundDiagA}
\alias{theoryRowpenUpperBoundDiagA}
\title{Penalty Parameter for Covariance Estimation Based on Theory}
\usage{
theoryRowpenUpperBoundDiagA(B, n1, n2, m)
}
\arguments{
\item{B}{row covariance matrix.}

\item{n1}{sample size of group one.}

\item{n2}{sample size of group two.}

\item{m}{number of columns of the data matrix (where the
data matrix is of size n by m, with n = n1 + n2).}
}
\value{
Returns a theoretically guided choice of the
glasso penalty parameter.
}
\description{
This function returns a theoretically-guided choice of
the glasso penalty parameter, treating the column
correlation matrix as the identity.
}
\examples{
# Define sample sizes
n1 <- 10
n2 <- 10
n <- n1 + n2
m <- 2e3
# Row covariance matrix (autoregressive of order 1)
B <- outer(1:n, 1:n, function(x, y) 0.8^abs(x - y))
# Calculate theoretically guided Gemini penalty.
rowpen <- theoryRowpenUpperBoundDiagA(B, n1, n2, m)
print(rowpen)
}
\references{
Joint mean and covariance estimation with unreplicated matrix-variate data
Michael Hornstein, Roger Fan, Kerby Shedden, Shuheng Zhou
(2018).  Joint mean and covariance estimation with
unreplicated matrix-variate data.  Journal of the American
Statistical Association
}
