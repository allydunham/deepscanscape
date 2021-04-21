test_that("normalisation works", {
  expect_equal(normalise_er(c(-4, -4, -2, -2, -1, 0, 1, 2, 2, 4, 4)),
               c(-1, -1, -0.5, -0.5, -0.25, 0, 0.25, 0.5, 0.5, 1, 1))
  expect_equal(normalise_er(c(-4, -4, -2, -2, -1, 0, 1, 2, 2, 4, 4), q = 0.4),
               c(-2, -2, -1, -1, -0.5, 0, 0.5, 1, 1, 2, 2))
})

test_that("transform_vamp works", {
  expect_equal(signif(transform_vamp(c(-0.05, 0, 0.05, 0.1, 0.5, 0.75, NA, 1, 1.5, 2)), 3),
               c(-4.39, -3.39, -2.81, -2.39, -0.807, -0.305, NA, 0.0671, 0.608, 1))
})

test_that("transformation works", {

  x <- c(0.1, 0.3, 0.5, NA, 1, 1.1)
  expect_equal(transform_er(x, trans = "log2"), log2(x))

  x <- c(-1, 0, 1, NA)
  expect_equal(transform_er(x, trans = "invert"), -x)

  x <- c(-0.05, 0, 0.05, 0.1, 0.5, 0.75, NA, 1, 1.5, 2)
  expect_equal(transform_er(x, trans = "vamp"), transform_vamp(x))

  expect_error(transform_er(1, trans = "Unrecognised_String"),
               "'arg' should be one of (\".*\",? ?)*")
  expect_error(transform_er(1, trans = NA),
               "Unrecognised transformation\\. Pass a supported string or a function")
})
