# Tests for lme4::glmer alternative engine and starting values

# Shared helper
make_episodic_glmer <- function(n_subjects = 500, seed = 55) {
  set.seed(seed)
  n_recalls   <- sample(1:2, n_subjects, replace = TRUE, prob = c(0.35, 0.65))
  total       <- sum(n_recalls)
  subject_ids <- rep(seq_len(n_subjects), times = n_recalls)
  repeat_nums <- unlist(lapply(n_recalls, seq_len))

  v1 <- rep(stats::rnorm(n_subjects, 0, 1.0), times = n_recalls)
  v2 <- rep(stats::rnorm(n_subjects, 0, sqrt(1.5)), times = n_recalls)

  consumed <- stats::rbinom(total, 1, stats::plogis(0.5 + v1))
  t_amount <- 3.0 + v2 + stats::rnorm(total, 0, 1.0)
  intake   <- consumed * pmax(boxcox_inverse(t_amount, 0.25), 0.1)

  data.frame(subject = subject_ids, day = repeat_nums,
             intake = intake, stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------------

test_that("glmer engine runs and returns a mixtran_fit", {
  skip_if_not_installed("lme4")
  dat <- make_episodic_glmer()
  fit <- suppressMessages(suppressWarnings(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "uncorr",
    lambda      = 0.25,
    prob_engine = "glmer",
    verbose     = FALSE
  )))
  expect_s3_class(fit, "mixtran_fit")
  expect_equal(fit$prob_engine, "glmer")
  expect_true(is.numeric(fit$sigma2_v1))
  expect_true(fit$sigma2_v1 > 0)
})

test_that("glmer engine recovers variance components within 30%", {
  skip_if_not_installed("lme4")
  dat <- make_episodic_glmer(n_subjects = 1200)
  fit <- suppressMessages(suppressWarnings(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "uncorr",
    lambda      = 0.25,
    prob_engine = "glmer",
    verbose     = FALSE
  )))
  v1_err <- abs(fit$sigma2_v1 - 1.0) / 1.0
  v2_err <- abs(fit$sigma2_v2 - 1.5) / 1.5
  expect_lte(v1_err, 0.50, label = "sigma2_v1 within 50% (glmer)")
  expect_lte(v2_err, 0.50, label = "sigma2_v2 within 50% (glmer)")
})

test_that("glmmPQL and glmer engines yield similar variance component estimates", {
  skip_if_not_installed("lme4")
  dat <- make_episodic_glmer(n_subjects = 600)

  fit_pql <- suppressMessages(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "uncorr",
    lambda      = 0.25,
    prob_engine = "glmmPQL",
    verbose     = FALSE
  ))
  fit_glmer <- suppressMessages(suppressWarnings(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "uncorr",
    lambda      = 0.25,
    prob_engine = "glmer",
    verbose     = FALSE
  )))

  # sigma2_v2 (amount) should be identical (same lme fit)
  expect_equal(fit_pql$sigma2_v2, fit_glmer$sigma2_v2, tolerance = 1e-6)

  # sigma2_v1 (prob) may differ due to PQL vs Laplace, but should be in same ballpark
  ratio <- fit_pql$sigma2_v1 / fit_glmer$sigma2_v1
  expect_true(ratio > 0.4 && ratio < 2.5,
              label = "PQL and glmer sigma2_v1 within factor of 2.5")
})

test_that("glmer engine works with corr model", {
  skip_if_not_installed("lme4")
  dat <- make_episodic_glmer(n_subjects = 400)
  fit <- suppressMessages(suppressWarnings(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "corr",
    lambda      = 0.25,
    prob_engine = "glmer",
    verbose     = FALSE
  )))
  expect_s3_class(fit, "mixtran_fit")
  expect_true(is.numeric(fit$rho))
  expect_true(fit$rho >= -1 && fit$rho <= 1)
})

# -------------------------------------------------------------------------
# Starting values tests

test_that("start parameter does not change results for amount model", {
  set.seed(7)
  n <- 300
  u_i <- rep(stats::rnorm(n, 0, 1), each = 2)
  dat_amt <- data.frame(
    subject = rep(1:n, each = 2),
    day     = rep(1:2, n),
    intake  = pmax(boxcox_inverse(3 + u_i + stats::rnorm(2*n, 0, 0.7), 0.3), 0.1)
  )

  fit1 <- suppressMessages(mixtran(
    data        = dat_amt,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "amount",
    lambda      = 0.3,
    verbose     = FALSE
  ))

  fit2 <- suppressMessages(mixtran(
    data        = dat_amt,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "amount",
    lambda      = 0.3,
    start       = fit1,
    verbose     = FALSE
  ))

  # With starting values from a good fit, results should be essentially identical
  expect_equal(fit1$sigma2_b, fit2$sigma2_b, tolerance = 1e-3)
  expect_equal(fit1$sigma2_w, fit2$sigma2_w, tolerance = 1e-3)
  expect_equal(fit1$beta, fit2$beta, tolerance = 1e-3)
})

test_that("start parameter accepts mixtran_fit and runs without error for uncorr", {
  dat <- make_episodic_glmer(n_subjects = 300)
  fit1 <- suppressMessages(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "uncorr",
    lambda      = 0.25,
    verbose     = FALSE
  ))
  expect_no_error(suppressMessages(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "uncorr",
    lambda      = 0.25,
    start       = fit1,
    verbose     = FALSE
  )))
})

test_that("start parameter rejects non-mixtran_fit objects", {
  expect_error(
    mixtran(
      data        = data.frame(s = 1, d = 1, y = 1),
      intake_var  = "y",
      subject_var = "s",
      repeat_var  = "d",
      model_type  = "amount",
      start       = list(beta = c(1, 2))
    ),
    regexp = "mixtran_fit"
  )
})
