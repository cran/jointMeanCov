context("Testing GeminiB and GeminiBPath")
library(jointMeanCov)

test_that("GeminiB estimates valid (inverse) covariance matrices", {
  n1 <- 5
  n2 <- 5
  n <- n1 + n2
  m <- 20
  X <- matrix(rnorm(n * m), nrow=n, ncol=m)
  rowpen <- sqrt(log(m) / n)
  out <- GeminiB(X, rowpen, penalize.diagonal=FALSE)
  expect_gt(min(eigen(out$B.hat.inv, only.values=TRUE)$values),
            0)
  expect_equal(out$corr.B.hat, cov2cor(out$B.hat))
})


test_that("GeminiBPath estimates valid (inverse) covariance matrices", {
  n1 <- 5
  n2 <- 5
  n <- n1 + n2
  m <- 20
  X <- matrix(rnorm(n * m), nrow=n, ncol=m)
  rowpen.list <- sqrt(log(m) / n) * c(1, 0.5, 0.1)
  out <- GeminiBPath(X, rowpen.list, penalize.diagonal=FALSE)

  expect_gt(
    min(eigen(out$B.hat.inv[, , 1], only.values=TRUE)$values), 0)
  expect_gt(
    min(eigen(out$B.hat.inv[, , 2], only.values=TRUE)$values), 0)
  expect_gt(
    min(eigen(out$B.hat.inv[, , 3], only.values=TRUE)$values), 0)

  expect_equal(out$corr.B.hat[, , 1], cov2cor(out$B.hat[, , 1]))
  expect_equal(out$corr.B.hat[, , 2], cov2cor(out$B.hat[, , 2]))
  expect_equal(out$corr.B.hat[, , 3], cov2cor(out$B.hat[, , 3]))
})
