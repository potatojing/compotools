test_that("cv.coat() runs successfullly", {

  x <- matrix(rpois(10*50, lambda = 1), nrow=10, ncol=50)
  x[x==0] = 0.5
  row_sums = rowSums(x)
  x = t(t(x)/row_sums)

  result = cv.coat(x)

  expect_true(all(c("sigma", "corr", "time", "lambda") %in% names(result)))
})


test_that("coat() runs successfullly", {

  x <- matrix(rpois(10*50, lambda = 1), nrow=10, ncol=50)
  x[x==0] = 0.5
  row_sums = rowSums(x)
  x = t(t(x)/row_sums)

  cv.result = cv.coat(x)
  lambda = cv.result$lambda
  result = coat(x, lambda)

  expect_true(all(c("sigma", "corr", "lambda") %in% names(result)))
})
