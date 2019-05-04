context("Test GLS functions")
library(jointMeanCov)

test_that("twoGroupDesignMatrix works", {
  D <- twoGroupDesignMatrix(1:2, 3:5)
  D.correct <- matrix(c(1, 0,
                        1, 0,
                        0, 1,
                        0, 1,
                        0, 1),
                      nrow=5, ncol=2, byrow=TRUE)
  expect_equal(D, D.correct)
})

test_that("GLSMeans works", {
  X <- matrix(1:12, nrow=4, ncol=3)
  D <- twoGroupDesignMatrix(1:2, 3:4)
  B.inv <- diag(4)
  out <- GLSMeans(X, D, B.inv)
  expect_equal(out, matrix(c(1.5, 5.5, 9.5,
                             3.5, 7.5, 11.5),
                           byrow=TRUE, nrow=2, ncol=3))
})
