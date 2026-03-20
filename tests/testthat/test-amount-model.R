# Helper: generate synthetic sodium-like data with known parameters
make_amount_data <- function(n_subjects    = 800,
                              true_lambda   = 0.30,
                              true_beta0    = 15.0,
                              true_sigma2_b = 4.0,
                              true_sigma2_w = 2.5,
                              seed          = 42) {
  set.seed(seed)

  n_recalls <- sample(1:2, n_subjects, replace = TRUE, prob = c(0.4, 0.6))
  total     <- sum(n_recalls)

  subject_ids  <- rep(seq_len(n_subjects), times = n_recalls)
  repeat_nums  <- unlist(lapply(n_recalls, seq_len))
  weekend      <- stats::rbinom(total, 1, 0.28)
  u_i          <- rep(stats::rnorm(n_subjects, 0, sqrt(true_sigma2_b)), times = n_recalls)
  t_intake     <- true_beta0 + 0.5 * weekend +
    (-0.2) * as.numeric(repeat_nums > 1) +
    u_i +
    stats::rnorm(total, 0, sqrt(true_sigma2_w))
  intake       <- pmax(boxcox_inverse(t_intake, true_lambda), 1)

  data.frame(
    SEQN       = subject_ids,
    day        = repeat_nums,
    sodium_mg  = intake,
    weekend    = weekend,
    stringsAsFactors = FALSE
  )
}

test_that("amount model recovers variance components within 15%", {
  dat <- make_amount_data()
  fit <- suppressMessages(mixtran(
    data       = dat,
    intake_var = "sodium_mg",
    subject_var = "SEQN",
    repeat_var  = "day",
    model_type  = "amount",
    weekend_var = "weekend",
    lambda      = 0.30,
    verbose     = FALSE
  ))

  expect_s3_class(fit, "mixtran_fit")
  expect_equal(fit$model_type, "amount")
  expect_true(fit$converged)

  b_err <- abs(fit$sigma2_b - 4.0) / 4.0
  w_err <- abs(fit$sigma2_w - 2.5) / 2.5
  expect_lte(b_err, 0.15, label = "sigma2_b within 15% of true")
  expect_lte(w_err, 0.15, label = "sigma2_w within 15% of true")
})

test_that("amount model returns required fields", {
  dat <- make_amount_data(n_subjects = 300)
  fit <- suppressMessages(mixtran(
    data       = dat,
    intake_var = "sodium_mg",
    subject_var = "SEQN",
    repeat_var  = "day",
    model_type  = "amount",
    lambda      = 0.30,
    verbose     = FALSE
  ))

  expect_true(!is.null(fit$beta))
  expect_true(!is.null(fit$sigma2_b))
  expect_true(!is.null(fit$sigma2_w))
  expect_true(!is.null(fit$var_ratio))
  expect_true(!is.null(fit$predicted))
  expect_true(!is.null(fit$model_fit))
  expect_true(!is.null(fit$n_subjects))
  expect_equal(fit$n_subjects, 300L)
})

test_that("mixtran stops on missing variables", {
  dat <- make_amount_data(n_subjects = 100)
  expect_error(
    suppressMessages(mixtran(
      data       = dat,
      intake_var = "sodium_mg",
      subject_var = "SEQN",
      repeat_var  = "day",
      covariates  = "nonexistent_var",
      lambda      = 0.30,
      verbose     = FALSE
    )),
    "not found in data"
  )
})

test_that("amount model auto-estimates lambda", {
  dat <- make_amount_data(n_subjects = 400)
  fit <- suppressMessages(mixtran(
    data       = dat,
    intake_var = "sodium_mg",
    subject_var = "SEQN",
    repeat_var  = "day",
    model_type  = "amount",
    lambda      = NULL,          # auto-estimate
    verbose     = FALSE
  ))
  expect_true(is.numeric(fit$lambda))
  expect_true(fit$lambda >= 0 && fit$lambda <= 1)
})

test_that("zero replacement works for amount model", {
  dat <- make_amount_data(n_subjects = 200)
  # Introduce some zeros
  dat$sodium_mg[sample(nrow(dat), 20)] <- 0
  # Should not error; should replace zeros silently
  fit <- suppressMessages(mixtran(
    data       = dat,
    intake_var = "sodium_mg",
    subject_var = "SEQN",
    repeat_var  = "day",
    model_type  = "amount",
    lambda      = 0.30,
    verbose     = FALSE
  ))
  expect_s3_class(fit, "mixtran_fit")
})

test_that("all-zero intake errors by default, returns NULL with skip_if_empty = TRUE", {
  dat <- make_amount_data(n_subjects = 100)
  dat$sodium_mg <- 0  # all zeros

  # Default: should error
  expect_error(
    suppressMessages(mixtran(
      data        = dat,
      intake_var  = "sodium_mg",
      subject_var = "SEQN",
      repeat_var  = "day",
      model_type  = "amount",
      lambda      = 0.30,
      verbose     = FALSE
    )),
    "No positive intake values"
  )

  # skip_if_empty = TRUE: should warn and return NULL
  result <- expect_warning(
    mixtran(
      data          = dat,
      intake_var    = "sodium_mg",
      subject_var   = "SEQN",
      repeat_var    = "day",
      model_type    = "amount",
      lambda        = 0.30,
      skip_if_empty = TRUE,
      verbose       = FALSE
    ),
    "no positive intake values"
  )
  expect_null(result)
})
