# Tests for bivariate usual intake and ratio estimation (distrib_bivariate)

make_bivariate_data <- function(n_subjects = 300, seed = 123) {
  set.seed(seed)
  n_recalls   <- sample(1:2, n_subjects, replace = TRUE, prob = c(0.4, 0.6))
  total       <- sum(n_recalls)
  subject_ids <- rep(seq_len(n_subjects), times = n_recalls)
  repeat_nums <- unlist(lapply(n_recalls, seq_len))

  # Correlated between-person random effects (rho = 0.4 between sodium and energy)
  rho    <- 0.4
  u1s    <- stats::rnorm(n_subjects, 0, sqrt(1.2))    # sodium RE
  u2s    <- rho * u1s / sqrt(1.2) * sqrt(0.8) +
            sqrt(1 - rho^2) * stats::rnorm(n_subjects, 0, sqrt(0.8))

  u1 <- rep(u1s, times = n_recalls)
  u2 <- rep(u2s, times = n_recalls)

  # Sodium (log-normal, lambda ~ 0): log-scale mean ~7.2 = ~1330 mg
  sodium_log  <- 7.2 + u1 + stats::rnorm(total, 0, sqrt(0.8))
  sodium      <- pmax(exp(sodium_log), 1)

  # Energy (log-normal): log-scale mean ~7.6 = ~2000 kcal
  energy_log  <- 7.6 + u2 + stats::rnorm(total, 0, sqrt(0.5))
  energy      <- pmax(exp(energy_log), 1)

  data.frame(
    SEQN    = subject_ids,
    day     = repeat_nums,
    sodium  = sodium,
    energy  = energy,
    stringsAsFactors = FALSE
  )
}

test_that("distrib_bivariate returns correct class and components", {
  dat <- make_bivariate_data()
  biv <- suppressMessages(suppressWarnings(
    distrib_bivariate(
      data        = dat,
      intake_var1 = "sodium",
      intake_var2 = "energy",
      subject_var = "SEQN",
      repeat_var  = "day",
      lambda1     = 0,
      lambda2     = 0,
      verbose     = FALSE,
      n_sims      = 20
    )
  ))

  expect_s3_class(biv, "bivariate_distrib")
  expect_true(!is.null(biv$component1))
  expect_true(!is.null(biv$component2))
  expect_true(!is.null(biv$ratio))
  expect_true(!is.null(biv$cross_rho))
  expect_true(biv$cross_rho > -1 && biv$cross_rho < 1)
})

test_that("distrib_bivariate component means are in plausible range", {
  dat <- make_bivariate_data()
  biv <- suppressMessages(suppressWarnings(
    distrib_bivariate(
      data        = dat,
      intake_var1 = "sodium",
      intake_var2 = "energy",
      subject_var = "SEQN",
      repeat_var  = "day",
      lambda1     = 0,
      lambda2     = 0,
      verbose     = FALSE,
      n_sims      = 50
    )
  ))

  mean1 <- biv$component1$Overall$mean
  mean2 <- biv$component2$Overall$mean

  # Sodium usual mean should be in 500-5000 range
  expect_true(mean1 > 500 && mean1 < 5000,
              label = sprintf("Sodium mean %.0f out of expected range", mean1))

  # Energy usual mean should be in 500-8000 range
  expect_true(mean2 > 500 && mean2 < 8000,
              label = sprintf("Energy mean %.0f out of expected range", mean2))
})

test_that("distrib_bivariate ratio is positive and finite", {
  dat <- make_bivariate_data()
  biv <- suppressMessages(suppressWarnings(
    distrib_bivariate(
      data        = dat,
      intake_var1 = "sodium",
      intake_var2 = "energy",
      subject_var = "SEQN",
      repeat_var  = "day",
      lambda1     = 0,
      lambda2     = 0,
      verbose     = FALSE,
      n_sims      = 30
    )
  ))

  ratio_mean <- biv$ratio$Overall$mean
  expect_true(is.finite(ratio_mean))
  expect_true(ratio_mean > 0)
  # Sodium/energy ratio (mg/kcal) should be in rough 0.3-5 range
  expect_true(ratio_mean > 0.1 && ratio_mean < 20,
              label = sprintf("Ratio mean %.3f out of expected range", ratio_mean))
})

test_that("bivariate_to_df returns 3-row data frame", {
  dat <- make_bivariate_data(n_subjects = 100, seed = 7)
  biv <- suppressMessages(suppressWarnings(
    distrib_bivariate(
      data        = dat,
      intake_var1 = "sodium",
      intake_var2 = "energy",
      subject_var = "SEQN",
      repeat_var  = "day",
      lambda1     = 0,
      lambda2     = 0,
      verbose     = FALSE,
      n_sims      = 10
    )
  ))

  df <- bivariate_to_df(biv)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 3)
  expect_true("component" %in% names(df))
  expect_true("mean" %in% names(df))
})

test_that("print.bivariate_distrib runs without error", {
  dat <- make_bivariate_data(n_subjects = 100, seed = 8)
  biv <- suppressMessages(suppressWarnings(
    distrib_bivariate(
      data        = dat,
      intake_var1 = "sodium",
      intake_var2 = "energy",
      subject_var = "SEQN",
      repeat_var  = "day",
      lambda1     = 0,
      lambda2     = 0,
      verbose     = FALSE,
      n_sims      = 10
    )
  ))
  expect_output(print(biv), "Bivariate Usual Intake")
})

test_that("distrib_bivariate custom fn_ratio works (e.g. sodium density = Na/1000*kcal)", {
  dat <- make_bivariate_data(n_subjects = 200, seed = 99)
  biv <- suppressMessages(suppressWarnings(
    distrib_bivariate(
      data        = dat,
      intake_var1 = "sodium",
      intake_var2 = "energy",
      subject_var = "SEQN",
      repeat_var  = "day",
      lambda1     = 0,
      lambda2     = 0,
      verbose     = FALSE,
      n_sims      = 20,
      ratio_label = "sodium_density",
      fn_ratio    = function(u1, u2) (u1 / u2) * 1000  # mg/1000 kcal
    )
  ))
  expect_equal(biv$ratio_label, "sodium_density")
  ratio_mean <- biv$ratio$Overall$mean
  expect_true(is.finite(ratio_mean) && ratio_mean > 0)
})
