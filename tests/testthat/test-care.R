test_that("cv.care() runs successfullly", {

  x <- matrix(rpois(10*50, lambda = 1), nrow=10, ncol=50)
  x[x==0] = 0.5
  row_sums = rowSums(x)
  x = t(t(x)/row_sums)

  result = cv.care(x, nlambda = 50, lambda_min = 0.01, nfold = 5)

  expect_true(all(c("Omega_hat", "lambda_op") %in% names(result)))
})

test_that("cv.care() runs successfullly", {

  x <- matrix(rpois(10*50, lambda = 1), nrow=10, ncol=50)
  x[x==0] = 0.5
  row_sums = rowSums(x)
  x = t(t(x)/row_sums)

  result = care(x, lambda_vec=rep(1,50))

  expect_true(all(c("Omega_hat", "lambda_vec") %in% names(result)))
})
