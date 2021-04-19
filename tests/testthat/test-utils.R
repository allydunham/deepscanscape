test_that("clamp clamps", {
  expect_equal(clamp(c(-2, -4, 0, 2, 3, 1)), c(-2, -4, 0, 2, 3, 1))
  expect_equal(clamp(c(-2, -4, 0, 2, 3, 1), lower = -1, upper = 1), c(-1, -1, 0, 1, 1, 1))
  expect_equal(clamp(c(-2, NA, 0, NA, 3, 1), lower = -1, upper = 1), c(-1, NA, 0, NA, 1, 1))
  expect_equal(clamp(c(-2, -Inf, 0, Inf, 3, 1), lower = -1, upper = 1), c(-1, -1, 0, 1, 1, 1))
  expect_equal(clamp(c(-2, -4, 0, 2, 3, 1), lower = -1), c(-1, -1, 0, 2, 3, 1))
  expect_equal(clamp(c(-2, -4, 0, 2, 3, 1), upper = 1), c(-2, -4, 0, 1, 1, 1))
  expect_equal(clamp(c(-2.2, -4.4, 0.5, 2.3, 3.1, 1.2), lower = -1.5, upper = 1.5), c(-1.5, -1.5, 0.5, 1.5, 1.5, 1.2))
})

test_that("cosine_distance_matrix is correct", {
  x <- matrix(c(1, 2, 4,
                2, 3, 5,
                1, 6, 2),
              nrow = 3, byrow = TRUE)

  x_na <- matrix(c(1, 2, 4,
                   2, NA, 5,
                   1, 6, 2),
                 nrow = 3, byrow = TRUE)

  y <- matrix(c(0.6, 0.12, -0.5,
                0.59, -1.2, -3.7,
                2.14, 0.31, 0.28,
                -0.82, 1.06, -0.60),
              nrow = 4, byrow = TRUE)

  x_res <- matrix(c(0.0000, 0.0423, 0.246,
                    0.0423, 0.0000, 0.225,
                    0.2460, 0.225, 0.000),
                  nrow = 3, byrow = TRUE)

  x_res_na <- matrix(c(0.000, NA, 0.246,
                          NA, NA, NA,
                       0.246, NA, 0.000),
                     nrow = 3, byrow = TRUE)

  xy_res <- matrix(c(0.604, 0.873, 0.373, 0.552,
                     0.562, 0.831, 0.336, 0.552,
                     0.480, 0.688, 0.394, 0.347),
                   nrow = 3, byrow = TRUE)

  expect_equal(signif(cosine_distance_matrix(x), 3), x_res)
  expect_equal(signif(cosine_distance_matrix(x, x), 3), x_res)
  expect_equal(signif(cosine_distance_matrix(x_na), 3), x_res_na)
  expect_equal(signif(cosine_distance_matrix(x, y), 3), xy_res)
})

test_that("permissive_positions works", {
  pos <- matrix(
    c(-0.8, 0.4, 0.9, -1.2, 1.4, 0.9, 0.1, 2.0, 1.1, 1.4, 1.1, 0.3, -0.3, 0.3, -0.6, 1.4, 0.1, -0.9, -0.1, 0.8,
      -0.3, 0.2, 0.1, -0.2, 0.1, 0.2, 0.1, 0.3, 0.2, 0.2, 0.1, 0.3, -0.3, 0.3, -0.2, 0.1, 0.1, -0.1, -0.1, 0.2),
    ncol = 20, byrow = TRUE
  )
  expect <- c(FALSE, TRUE)
  expect_equal(permissive_positions(pos), expect)
})
