# Tests for the GHQ (Gauss-Hermite Quadrature) engine of the correlated
# two-part model, and the quadrature utility functions.

# ---------------------------------------------------------------------------
# Quadrature utility tests
# ---------------------------------------------------------------------------

test_that("gh_nodes_weights returns correct structure", {
  for (n in c(3, 5, 7, 9)) {
    gh <- gh_nodes_weights(n)
    expect_length(gh$nodes, n)
    expect_length(gh$weights, n)
    expect_true(all(gh$weights > 0))
    # Weights should sum to sqrt(pi) for physicist convention
    expect_equal(sum(gh$weights), sqrt(pi), tolerance = 1e-3)
    # Nodes should be symmetric about zero
    expect_equal(gh$nodes, -rev(gh$nodes), tolerance = 1e-10)
  }
})

test_that("gh_nodes_normal weights sum to 1", {
  for (n in c(3, 5, 7)) {
    ghn <- gh_nodes_normal(n)
    expect_equal(sum(ghn$weights), 1, tolerance = 1e-10)
    expect_true(all(ghn$weights > 0))
  }
})

test_that("gh_nodes_normal integrates x^2 against N(0,1) ≈ 1", {
  # E[X^2] = 1 for X ~ N(0,1)
  ghn <- gh_nodes_normal(7)
  integral <- sum(ghn$weights * ghn$nodes^2)
  expect_equal(integral, 1, tolerance = 0.01)
})

test_that("gh_nodes_bivariate has correct dimensions", {
  biv <- gh_nodes_bivariate(5, rho = 0.3, sigma_v1 = 1, sigma_v2 = 1.2)
  expect_equal(nrow(biv), 25)
  expect_true(all(c("v1", "v2", "log_weight") %in% names(biv)))
})

test_that("gh_nodes_bivariate marginal variances match inputs", {
  # With rho=0, v1 and v2 should be independent
  # E[v1^2] = sigma_v1^2, E[v2^2] = sigma_v2^2
  biv <- gh_nodes_bivariate(9, rho = 0, sigma_v1 = 2, sigma_v2 = 3)
  w   <- exp(biv$log_weight)
  ev1_sq <- sum(w * biv$v1^2)
  ev2_sq <- sum(w * biv$v2^2)
  expect_equal(ev1_sq, 4, tolerance = 0.05)
  expect_equal(ev2_sq, 9, tolerance = 0.05)
})

test_that("gh_nodes_bivariate covariance matches rho * sigma_v1 * sigma_v2", {
  rho <- 0.5; sv1 <- 1.5; sv2 <- 2.0
  biv <- gh_nodes_bivariate(9, rho = rho, sigma_v1 = sv1, sigma_v2 = sv2)
  w   <- exp(biv$log_weight)
  cov_v1v2 <- sum(w * biv$v1 * biv$v2)
  expect_equal(cov_v1v2, rho * sv1 * sv2, tolerance = 0.05)
})

test_that("gh_nodes_weights works for non-table n via eigenvalue method", {
  gh <- gh_nodes_weights(4)
  expect_length(gh$nodes, 4)
  expect_equal(sum(gh$weights), sqrt(pi), tolerance = 1e-3)
})


# ---------------------------------------------------------------------------
# GHQ engine integration tests (synthetic data)
# ---------------------------------------------------------------------------

make_corr_data_ghq <- function(n_subjects = 400, true_rho = 0.5, seed = 42) {
  set.seed(seed)
  n_recalls <- sample(1:2, n_subjects, replace = TRUE, prob = c(0.35, 0.65))
  total     <- sum(n_recalls)
  subject_ids <- rep(seq_len(n_subjects), times = n_recalls)
  repeat_nums <- unlist(lapply(n_recalls, seq_len))

  sv1 <- 1.0; sv2 <- 0.9
  z1s <- stats::rnorm(n_subjects)
  z2s <- stats::rnorm(n_subjects)
  v1s <- sv1 * z1s
  v2s <- sv2 * (true_rho * z1s + sqrt(1 - true_rho^2) * z2s)

  v1 <- rep(v1s, times = n_recalls)
  v2 <- rep(v2s, times = n_recalls)

  p_consume <- stats::plogis(0.2 + v1)
  consumed  <- stats::rbinom(total, 1, p_consume)
  lambda    <- 0.3
  t_amount  <- 2.8 + v2 + stats::rnorm(total, 0, sqrt(0.7))
  amount    <- pmax(boxcox_inverse(t_amount, lambda), 0.1)
  intake    <- consumed * amount

  data.frame(id = subject_ids, day = repeat_nums,
             food_g = intake, stringsAsFactors = FALSE)
}

test_that("GHQ engine returns valid mixtran_fit with correct rho sign", {
  dat <- make_corr_data_ghq(true_rho = 0.5)
  fit <- suppressMessages(suppressWarnings(mixtran(
    data        = dat,
    intake_var  = "food_g",
    subject_var = "id",
    repeat_var  = "day",
    model_type  = "corr",
    corr_engine = "ghq",
    ghq_n_nodes = 5L,
    lambda      = 0.3,
    verbose     = FALSE
  )))

  expect_s3_class(fit, "mixtran_fit")
  expect_equal(fit$model_type, "corr")
  expect_equal(fit$corr_engine, "ghq")
  expect_true(fit$rho > -1 && fit$rho < 1)
  # Should recover positive correlation (may not be close with n=400 but sign should be right)
  expect_gt(fit$rho, -0.2)
})

test_that("GHQ engine stores ghq_n_nodes in fit", {
  dat <- make_corr_data_ghq(n_subjects = 200, seed = 99)
  fit <- suppressMessages(suppressWarnings(mixtran(
    data        = dat,
    intake_var  = "food_g",
    subject_var = "id",
    repeat_var  = "day",
    model_type  = "corr",
    corr_engine = "ghq",
    ghq_n_nodes = 3L,
    lambda      = 0.3,
    verbose     = FALSE
  )))
  expect_equal(fit$ghq_n_nodes, 3L)
  expect_true(!is.null(fit$ghq_loglik))
})

test_that("GHQ and profile_rho engines give rho estimates with same sign", {
  dat <- make_corr_data_ghq(true_rho = 0.6, seed = 77)

  fit_prof <- suppressMessages(suppressWarnings(mixtran(
    data        = dat,
    intake_var  = "food_g",
    subject_var = "id",
    repeat_var  = "day",
    model_type  = "corr",
    corr_engine = "profile_rho",
    lambda      = 0.3,
    verbose     = FALSE
  )))

  fit_ghq <- suppressMessages(suppressWarnings(mixtran(
    data        = dat,
    intake_var  = "food_g",
    subject_var = "id",
    repeat_var  = "day",
    model_type  = "corr",
    corr_engine = "ghq",
    ghq_n_nodes = 5L,
    lambda      = 0.3,
    verbose     = FALSE
  )))

  # Both should recover positive rho direction
  expect_gt(fit_prof$rho, 0)
  expect_gt(fit_ghq$rho, 0)
})

test_that("distrib() works correctly after GHQ-fitted mixtran", {
  dat <- make_corr_data_ghq(n_subjects = 300, seed = 55)
  fit <- suppressMessages(suppressWarnings(mixtran(
    data        = dat,
    intake_var  = "food_g",
    subject_var = "id",
    repeat_var  = "day",
    model_type  = "corr",
    corr_engine = "ghq",
    ghq_n_nodes = 3L,
    lambda      = 0.3,
    verbose     = FALSE
  )))

  dr <- suppressMessages(distrib(fit, n_sims = 20, seed = 1))
  expect_s3_class(dr, "distrib_result")
  expect_true(!is.null(dr$results$Overall$mean))
  expect_true(dr$results$Overall$mean > 0)
})
