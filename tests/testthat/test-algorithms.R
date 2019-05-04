context("Testing Algorithms 1, 2, and Stability Selection")
library(jointMeanCov)

test_that("Algorithm 1 (group centering)", {
  n1 <- 5
  n2 <- 5
  n <- n1 + n2
  group.one.indices <- 1:5
  group.two.indices <- 6:10
  m <- 20
  X <- matrix(rnorm(n * m), nrow=n, ncol=m)
  rowpen <- sqrt(log(m) / n)
  out <- jointMeanCovGroupCen(X, group.one.indices, rowpen)
  expect_equal(length(out$gamma.hat), ncol(X))
  expect_equal(nrow(out$B.hat.inv), nrow(X))
  expect_equal(
    abs(out$gamma.hat),
    abs(out$gls.group.means[2, ] - out$gls.group.means[1, ]))

  # Test when the identity matrix is passed as B.inv.
  out.diagB <- jointMeanCovGroupCen(X, group.one.indices, rowpen,
                                    B.inv=diag(n))
  expect_equal(out.diagB$B.hat.inv, diag(n))
  expect_equal(out.diagB$corr.B.hat.inv, diag(n))
  expect_equal(out.diagB$gls.group.means[1, ],
               apply(X[group.one.indices, ], 2, mean))
  expect_equal(out.diagB$gls.group.means[2, ],
               apply(X[group.two.indices, ], 2, mean))
})

test_that("Algorithm 2 (model selection centering)", {
  n1 <- 5
  n2 <- 5
  n <- n1 + n2
  group.one.indices <- 1:5
  group.two.indices <- 6:10
  m <- 20
  X <- matrix(rnorm(n * m), nrow=n, ncol=m)
  rowpen <- sqrt(log(m) / n)
  out <- jointMeanCovModSelCen(X, group.one.indices, rowpen)
  expect_equal(length(out$gamma.hat), ncol(X))
  expect_equal(nrow(out$B.hat.inv), nrow(X))
  expect_equal(
    abs(out$gamma.hat),
    abs(out$gls.group.means[2, ] - out$gls.group.means[1, ]))

  # Test a small threshold.
  n1 <- 5
  n2 <- 5
  n <- n1 + n2
  group.one.indices <- 1:5
  group.two.indices <- 6:10
  m <- 20
  X <- matrix(rnorm(n * m), nrow=n, ncol=m)
  rowpen <- sqrt(log(m) / n)
  # When the threshold is zero, model selection centering
  # should have the same output as group centering.
  out.smallThresh <- jointMeanCovModSelCen(
    X, group.one.indices, rowpen, thresh=1e-7)
  out.groupCen <- jointMeanCovGroupCen(X, group.one.indices,
                                       rowpen)
  expect_equal(out.smallThresh$gamma.hat,
               out.groupCen$gamma.hat)

})

test_that("Stability selection", {
  n1 <- 5
  n2 <- 5
  n <- n1 + n2
  group.one.indices <- 1:5
  group.two.indices <- 6:10
  m <- 20
  M <- matrix(0, nrow=n, ncol=m)
  M[1:n1, 1:3] <- 2
  M[(n1 + 1):n2, 1:3] <- -2
  X <- matrix(rnorm(n * m), nrow=n, ncol=m) + M
  rowpen <- sqrt(log(m) / n)
  n.genes.to.group.center <- c(10, 5, 2)
  out <- jointMeanCovStability(
    X, group.one.indices, rowpen, c(2e3, n.genes.to.group.center))
  expect_equal(length(out$n.genes.to.group.center),
               length(out$betaHat))
  # check that out$n.genes.to.group.center is a decreasing sequence
  expect_lte(max(diff(out$n.genes.to.group.center)), 0)

  out <- jointMeanCovStability(X, group.one.indices, rowpen)
})
