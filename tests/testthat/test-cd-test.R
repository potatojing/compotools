test_that("paired cd.test() runs successfully", {

  x <- matrix(rpois(50*10, lambda=10), nrow = 50, ncol = 10)
  y <- matrix(rpois(50*10, lambda=10), nrow = 50, ncol = 10)

  x = t(t(x)/colSums(x))
  y = t(t(y)/colSums(y))

  p <- cd.test(x,y,paired = TRUE)
  expect_true(p >= 0 && p <= 1)
})
test_that("unpaired cd.test() runs successfully", {

  x <- matrix(rpois(50*20, lambda=10), nrow = 50, ncol = 20)
  y <- matrix(rpois(50*10, lambda=10), nrow = 50, ncol = 10)

  x = t(t(x)/colSums(x))
  y = t(t(y)/colSums(y))

  p <- cd.test(x,y,paired = FALSE)
  expect_true(p >= 0 && p <= 1)
})
