test_that("clamp clamps", {
  expect_equal(clamp(c(-2, -4, 0, 2, 3, 1)), c(-2, -4, 0, 2, 3, 1))
  expect_equal(clamp(c(-2, -4, 0, 2, 3, 1), lower = -1, upper = 1), c(-1, -1, 0, 1, 1, 1))
  expect_equal(clamp(c(-2, NA, 0, NA, 3, 1), lower = -1, upper = 1), c(-1, NA, 0, NA, 1, 1))
  expect_equal(clamp(c(-2, -Inf, 0, Inf, 3, 1), lower = -1, upper = 1), c(-1, -1, 0, 1, 1, 1))
  expect_equal(clamp(c(-2, -4, 0, 2, 3, 1), lower = -1), c(-1, -1, 0, 2, 3, 1))
  expect_equal(clamp(c(-2, -4, 0, 2, 3, 1), upper = 1), c(-2, -4, 0, 1, 1, 1))
  expect_equal(clamp(c(-2.2, -4.4, 0.5, 2.3, 3.1, 1.2), lower = -1.5, upper = 1.5), c(-1.5, -1.5, 0.5, 1.5, 1.5, 1.2))
})
