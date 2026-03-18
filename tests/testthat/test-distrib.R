# Helper: fast amount model fit
fit_amount_for_distrib <- function(n_subjects = 600, seed = 42) {
  set.seed(seed)
  n_recalls <- sample(1:2, n_subjects, replace = TRUE, prob = c(0.4, 0.6))
  total     <- sum(n_recalls)
  sids      <- rep(seq_len(n_subjects), times = n_recalls)
  reps      <- unlist(lapply(n_recalls, seq_len))
  u_i       <- rep(stats::rnorm(n_subjects, 0, 2), times = n_recalls)
  t_intake  <- 15 + u_i + stats::rnorm(total, 0, sqrt(2.5))
  intake    <- pmax(boxcox_inverse(t_intake, 0.30), 1)
  wt        <- rep(stats::runif(n_subjects, 0.5, 2.0), times = n_recalls)
  dat       <- data.frame(SEQN = sids, day = reps, sodium_mg = intake,
                           wt = wt, stringsAsFactors = FALSE)

  suppressMessages(mixtran(
    data       = dat,
    intake_var = "sodium_mg",
    subject_var = "SEQN",
    repeat_var  = "day",
    model_type  = "amount",
    weight_var  = "wt",
    lambda      = 0.30,
    verbose     = FALSE
  ))
}

test_that("distrib returns distrib_result with correct structure", {
  fit <- fit_amount_for_distrib()
  res <- suppressMessages(distrib(fit, n_sims = 50, seed = 1))

  expect_s3_class(res, "distrib_result")
  expect_true(!is.null(res$results))
  expect_true("Overall" %in% names(res$results))
  expect_true(!is.null(res$results$Overall$mean))
  expect_true(!is.null(res$results$Overall$sd))
  expect_true(!is.null(res$results$Overall$percentiles))
})

test_that("usual intake distribution is narrower than observed", {
  fit <- fit_amount_for_distrib()
  res <- suppressMessages(distrib(fit, n_sims = 100, seed = 7))

  usual_sd <- res$results$Overall$sd
  obs_sd   <- stats::sd(fit$predicted$intake)

  expect_lt(usual_sd, obs_sd,
            label = "Usual intake SD < observed SD (within-person variance removed)")
})

test_that("distrib mean is within 10% of theoretical", {
  fit <- fit_amount_for_distrib()
  res <- suppressMessages(distrib(fit, n_sims = 200, seed = 99))

  # Theoretical mean: E[T^{-1}(beta0 + u)] where u ~ N(0, sigma2_b)
  theo_samp  <- boxcox_inverse(
    stats::rnorm(100000, fit$beta[["(Intercept)"]], sqrt(fit$sigma2_b)),
    fit$lambda
  )
  theo_mean  <- mean(theo_samp)
  est_mean   <- res$results$Overall$mean
  pct_err    <- abs(est_mean - theo_mean) / theo_mean

  expect_lte(pct_err, 0.10, label = "Mean within 10% of theoretical")
})

test_that("distrib cutpoints work correctly", {
  fit <- fit_amount_for_distrib()
  res <- suppressMessages(distrib(
    fit, n_sims = 100, seed = 5,
    cutpoints = c(1500, 2300, 3500)
  ))

  r <- res$results$Overall
  expect_true(!is.null(r$cutpoint_below))
  expect_true(!is.null(r$cutpoint_above))
  expect_equal(length(r$cutpoint_below), 3L)
  # proportions must be in [0,1]
  expect_true(all(r$cutpoint_below >= 0 & r$cutpoint_below <= 1))
  # above + below = 1
  expect_equal(unname(r$cutpoint_below + r$cutpoint_above),
               rep(1, 3), tolerance = 1e-10)
})

test_that("distrib_to_df returns a data frame with one row per subgroup", {
  fit <- fit_amount_for_distrib()
  res <- suppressMessages(distrib(fit, n_sims = 50, seed = 3))
  df  <- distrib_to_df(res)

  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 1L)            # "Overall" only
  expect_true("mean" %in% names(df))
  expect_true("sd"   %in% names(df))
})

test_that("distrib is reproducible with the same seed", {
  fit <- fit_amount_for_distrib()
  r1  <- suppressMessages(distrib(fit, n_sims = 50, seed = 42))
  r2  <- suppressMessages(distrib(fit, n_sims = 50, seed = 42))

  expect_equal(r1$results$Overall$mean, r2$results$Overall$mean)
  expect_equal(r1$results$Overall$percentiles, r2$results$Overall$percentiles)
})

test_that("weighted_quantiles returns correct quantiles", {
  # With equal weights, should match base R quantile type 7
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  w <- rep(1, 10)
  q <- nciusual:::weighted_quantiles(x, w, c(0.1, 0.5, 0.9))

  expect_equal(q[2], 5.5, tolerance = 0.1)   # median
  expect_true(q[1] < q[2] && q[2] < q[3])    # monotone
})
