test_that("cdmm() runs successfullly", {

  y <- rnorm(50)
  x <- matrix(rnorm(50*10), nrow=50, ncol=10)
  row_sums = rowSums(x)
  x = t(t(x)/row_sums)

  result <- cdmm(y = y, x = x)

  expect_true(all(c("bet", "lam", "int") %in% names(result)))
})

test_that("cv.cdmm() runs successfullly", {

  y <- rnorm(50)
  x <- matrix(rnorm(50*10), nrow=50, ncol=10)
  row_sums = rowSums(x)
  x = t(t(x)/row_sums)

  result <- cv.cdmm(y = y, x = x)

  expect_true(all(c("bet", "lam", "int", "foldid") %in% names(result)))
})

test_that("cv.cdmm(refit=TRUE) runs successfullly", {

  y <- rnorm(50)
  x <- matrix(rnorm(50*10), nrow=50, ncol=10)
  row_sums = rowSums(x)
  x = t(t(x)/row_sums)

  result <- cv.cdmm(y = y, x = x, refit = TRUE)

  expect_true(all(c("bet", "lam", "int", "foldid") %in% names(result)))
})

test_that("cv.cdmm(std=FALSE) runs successfullly", {

  y <- rnorm(50)
  x <- matrix(rnorm(50*10), nrow=50, ncol=10)
  row_sums = rowSums(x)
  x = t(t(x)/row_sums)

  result <- cv.cdmm(y = y, x = x, std=FALSE)

  expect_true(all(c("bet", "lam", "int", "foldid") %in% names(result)))
})
