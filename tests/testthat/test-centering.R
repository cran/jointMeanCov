context("Centering helper functions")
library(jointMeanCov)

test_that("centerDataTwoGroupsByIndices works correctly", {
  X <- matrix(1:12, nrow=4, ncol=3)
  X.cen <- matrix(c(-0.5, -0.5, -0.5,
                    0.5, 0.5, 0.5,
                    -0.5, -0.5, -0.5,
                    0.5, 0.5, 0.5),
                  nrow=4, ncol=3,
                  byrow = TRUE)
  expect_equal(
    centerDataTwoGroupsByIndices(
      X, group.one.indices=1:2, group.two.indices=3:4),
    X.cen)
})

test_that("centerDataTwoGroupsByModelSelection works correctly", {
  # In the following test, we group center the first two
  # columns and globally center the third column.
  X <- matrix(1:12, nrow=4, ncol=3)
  X.cen <- matrix(c(-0.5, -0.5, -1.5,
                    0.5, 0.5, -0.5,
                    -0.5, -0.5, 0.5,
                    0.5, 0.5, 1.5),
                  nrow=4, ncol=3,
                  byrow = TRUE)
  expect_equal(
    centerDataTwoGroupsByModelSelection(
      X, group.one.indices=1:2, group.two.indices=3:4,
      within.group.indices=1:2),
    X.cen)
})

test_that("centerDataGLSModelSelection works", {
  n <- 4
  m <- 3
  X <- matrix(1:12, nrow=n, ncol=m)
  out <- centerDataGLSModelSelection(
    X, B.inv=diag(n), group.one.indices=1:2,
    group.two.indices=3:4,
    group.cen.indices=1:2)
  expect_equal(out$X.cen,
               matrix(c(-0.5, -0.5, -1.5,
                        0.5, 0.5, -0.5,
                        -0.5, -0.5, 0.5,
                        0.5, 0.5, 1.5),
                      byrow=TRUE, 4, 3))
  expect_equal(out$global.means.gls, matrix(10.5, 1, 1))
  expect_equal(out$group.means.gls,
               matrix(c(1.5, 5.5,
                        3.5, 7.5),
                      byrow=TRUE, 2, 2))
})
