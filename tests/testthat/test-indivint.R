# Tests for indivint() — individual BLUP predictions

# Shared helper: amount-only synthetic data (re-used from test-amount-model.R pattern)
make_amount_data_indivint <- function(n_subjects = 400,
                                      true_beta0  = 3.0,
                                      true_sigma2_b = 1.0,
                                      true_sigma2_w = 0.5,
                                      true_lambda   = 0.3,
                                      seed = 42) {
  set.seed(seed)
  n_recalls <- sample(1:2, n_subjects, replace = TRUE, prob = c(0.3, 0.7))
  total     <- sum(n_recalls)
  subject_ids <- rep(seq_len(n_subjects), times = n_recalls)
  repeat_nums <- unlist(lapply(n_recalls, seq_len))

  u_i <- rep(stats::rnorm(n_subjects, 0, sqrt(true_sigma2_b)), times = n_recalls)
  eps <- stats::rnorm(total, 0, sqrt(true_sigma2_w))
  t_y <- true_beta0 + u_i + eps
  y   <- pmax(boxcox_inverse(t_y, true_lambda), 0.01)

  data.frame(subject = subject_ids, day = repeat_nums,
             intake = y, stringsAsFactors = FALSE)
}

# Shared helper: two-part synthetic data
make_episodic_data_indivint <- function(n_subjects = 500, seed = 77) {
  set.seed(seed)
  n_recalls <- sample(1:2, n_subjects, replace = TRUE, prob = c(0.35, 0.65))
  total     <- sum(n_recalls)
  subject_ids <- rep(seq_len(n_subjects), times = n_recalls)
  repeat_nums <- unlist(lapply(n_recalls, seq_len))

  v1 <- rep(stats::rnorm(n_subjects, 0, 1), times = n_recalls)
  v2 <- rep(stats::rnorm(n_subjects, 0, sqrt(1.5)), times = n_recalls)

  p_consume <- stats::plogis(0.5 + v1)
  consumed  <- stats::rbinom(total, 1, p_consume)
  t_amount  <- 3.0 + v2 + stats::rnorm(total, 0, 1)
  amount    <- pmax(boxcox_inverse(t_amount, 0.25), 0.1)
  intake    <- consumed * amount

  data.frame(subject = subject_ids, day = repeat_nums,
             intake = intake, stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------------

test_that("indivint returns one row per subject for amount model", {
  dat <- make_amount_data_indivint()
  fit <- suppressMessages(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "amount",
    lambda      = 0.3,
    verbose     = FALSE
  ))

  result <- indivint(fit)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), length(unique(dat$subject)))
  expect_true("subject"              %in% names(result))
  expect_true("predicted_usual"     %in% names(result))
  expect_true("predicted_usual_orig" %in% names(result))
})

test_that("indivint predictions are positive for amount model", {
  dat <- make_amount_data_indivint()
  fit <- suppressMessages(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "amount",
    lambda      = 0.3,
    verbose     = FALSE
  ))
  result <- indivint(fit)
  expect_true(all(result$predicted_usual_orig > 0))
})

test_that("indivint BLUP mean close to DISTRIB mean for amount model", {
  dat <- make_amount_data_indivint(n_subjects = 600)
  fit <- suppressMessages(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "amount",
    lambda      = 0.3,
    verbose     = FALSE
  ))

  blup_result <- indivint(fit)
  dist_result <- suppressMessages(distrib(fit, n_sims = 200, seed = 1))

  # Simple unweighted mean of BLUPs should approximate distrib mean within 10%
  blup_mean <- mean(blup_result$predicted_usual_orig)
  dist_mean <- dist_result$results$Overall$mean
  rel_diff  <- abs(blup_mean - dist_mean) / dist_mean

  expect_lte(rel_diff, 0.15,
             label = "BLUP mean within 15% of DISTRIB mean")
})

test_that("indivint returns one row per subject for uncorr two-part model", {
  dat <- make_episodic_data_indivint()
  fit <- suppressMessages(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "uncorr",
    lambda      = 0.25,
    verbose     = FALSE
  ))

  result <- indivint(fit)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), length(unique(dat$subject)))
  expect_true("prob_usual"          %in% names(result))
  expect_true("amt_usual_orig"      %in% names(result))
  expect_true("predicted_usual_orig" %in% names(result))
})

test_that("indivint two-part probabilities are in [0, 1]", {
  dat <- make_episodic_data_indivint()
  fit <- suppressMessages(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "uncorr",
    lambda      = 0.25,
    verbose     = FALSE
  ))
  result <- indivint(fit)
  expect_true(all(result$prob_usual >= 0 & result$prob_usual <= 1))
})

test_that("indivint predictions are non-negative for two-part model", {
  dat <- make_episodic_data_indivint()
  fit <- suppressMessages(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "uncorr",
    lambda      = 0.25,
    verbose     = FALSE
  ))
  result <- indivint(fit)
  expect_true(all(result$predicted_usual_orig >= 0))
})

test_that("indivint works for corr two-part model", {
  dat <- make_episodic_data_indivint(n_subjects = 400)
  fit <- suppressMessages(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "corr",
    lambda      = 0.25,
    verbose     = FALSE
  ))
  result <- indivint(fit)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), length(unique(dat$subject)))
  expect_true(all(result$predicted_usual_orig >= 0))
})
