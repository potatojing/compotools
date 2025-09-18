test_that("cdmm() 能正常处理模拟数据", {
  # 生成模拟数据
  y <- rnorm(50)
  x <- matrix(rnorm(50*10), nrow=50, ncol=10)
  row_sums = rowSums(x)
  x = t(t(x)/row_sums)
  # 调用函数，验证不报错
  result <- cdmm(y = y, x = x)
  # 验证输出包含预期的字段（sol、lam、int）
  expect_true(all(c("sol", "lam", "int") %in% names(result)))
})
