test_that("cdmm() runs successfullly", {

  y <- rnorm(50)
  x <- matrix(rnorm(50*10), nrow=50, ncol=10)
  row_sums = rowSums(x)
  x = t(t(x)/row_sums)

  result <- cdmm(y = y, x = x)

  expect_true(all(c("sol", "lam", "int") %in% names(result)))
})
