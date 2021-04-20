test_that("normalisation works", {
  expect_equal(normalise_er(c(-4, -4, -2, -2, -1, 0, 1, 2, 2, 4, 4)),
               c(-1, -1, -0.5, -0.5, -0.25, 0, 0.25, 0.5, 0.5, 1, 1))
  expect_equal(normalise_er(c(-4, -4, -2, -2, -1, 0, 1, 2, 2, 4, 4), q = 0.4),
               c(-2, -2, -1, -1, -0.5, 0, 0.5, 1, 1, 2, 2))
})
