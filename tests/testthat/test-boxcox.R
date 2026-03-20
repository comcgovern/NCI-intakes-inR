test_that("boxcox_transform and boxcox_inverse are inverses", {
  y <- c(100, 500, 1000, 2500, 5000)

  for (lam in c(0, 0.1, 0.25, 0.5, 1.0)) {
    z     <- boxcox_transform(y, lam)
    y_hat <- boxcox_inverse(z, lam)
    expect_equal(y_hat, y, tolerance = 1e-8,
                 label = sprintf("lambda=%.2f", lam))
  }
})

test_that("boxcox_transform with lambda=0 is log transform", {
  y <- c(1, exp(1), exp(2))
  z <- boxcox_transform(y, 0)
  expect_equal(z, c(0, 1, 2), tolerance = 1e-10)
})

test_that("boxcox_transform rejects non-positive values", {
  expect_error(boxcox_transform(c(1, 0, 2), lambda = 0.5))
  expect_error(boxcox_transform(c(1, -1, 2), lambda = 0.5))
})

test_that("boxcox_inverse handles edge cases without crashing", {
  # Very negative z can give near-zero back-transformed values
  z <- c(-1000, -10, 0, 10)
  y <- boxcox_inverse(z, lambda = 0.5)
  expect_true(all(is.finite(y)))
  expect_true(all(y >= 0))
})

test_that("find_optimal_lambda recovers lambda near 0.5 for chi-sq-like data", {
  set.seed(77)
  y <- (stats::rnorm(3000, mean = 5, sd = 2))^2
  y <- y[y > 0]
  result <- find_optimal_lambda(y, lambda_grid = seq(0.01, 1.0, by = 0.05))

  expect_type(result$lambda, "double")
  expect_type(result$lambda_raw, "double")
  expect_lte(abs(result$lambda_raw - 0.5), 0.2)
})

test_that("find_optimal_lambda uses log transform when lambda_raw < 0.15", {
  # For highly right-skewed data the optimal lambda will be near 0
  set.seed(42)
  y <- stats::rexp(2000, rate = 0.001)^3  # very skewed
  result <- find_optimal_lambda(y, lambda_grid = seq(0.01, 1.0, by = 0.01))

  if (result$used_log) {
    expect_equal(result$lambda, 0)
  } else {
    expect_gte(result$lambda_raw, 0.15)
  }
})

test_that("find_optimal_lambda errors on empty input", {
  expect_error(find_optimal_lambda(numeric(0)), "No positive values")
  expect_error(find_optimal_lambda(c(0, 0, NA)), "No positive values")
})
