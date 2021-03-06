% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GLSMeans.R
\name{GLSMeans}
\alias{GLSMeans}
\title{Generalized Least Squares}
\usage{
GLSMeans(X, D, B.inv)
}
\arguments{
\item{X}{data matrix.}

\item{D}{design matrix.}

\item{B.inv}{inverse covariance matrix.}
}
\value{
Returns the estimated parameters of the linear model,
a matrix of dimensions k by m, where k is the number of
columns of D, and m is the number of columns of X.
}
\description{
This function applies generalized least squares to estimate
the  unknown parameters of a linear model X = D beta + E,
where X has dimension n by m, D has dimension n by k, and
beta has dimension k by m.
}
\details{
\bold{Example}
\preformatted{X <- matrix(1:12, nrow=4, ncol=3)
D <- twoGroupDesignMatrix(1:2, 3:4)
B.inv <- diag(4)
beta.hat <- GLSMeans(X, D, B.inv)
}
}
