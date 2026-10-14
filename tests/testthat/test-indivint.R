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

  # E[T | recalls] averages to E[T] (law of total expectation), so the
  # unweighted mean of INDIVINT predictions approximates the DISTRIB mean
  blup_mean <- mean(blup_result$predicted_usual_orig)
  dist_mean <- dist_result$results$Overall$mean
  rel_diff  <- abs(blup_mean - dist_mean) / dist_mean

  expect_lte(rel_diff, 0.03,
             label = "INDIVINT mean within 3% of DISTRIB mean")
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

# ---------------------------------------------------------------------------
# Regression: two-part BLUP predictions must vary across subjects
#
# Guards against the random-effect extraction bug where the per-subject
# random effects collapsed to a single unnamed value, so every subject
# lookup missed and all subjects received the identical population-level
# prediction — silently defeating the purpose of INDIVINT (BLUP).
# ---------------------------------------------------------------------------

test_that("indivint two-part predictions vary across subjects (BLUP is used)", {
  dat <- make_episodic_data_indivint(n_subjects = 500, seed = 77)
  fit <- suppressWarnings(suppressMessages(mixtran(
    data        = dat,
    intake_var  = "intake",
    subject_var = "subject",
    repeat_var  = "day",
    model_type  = "uncorr",
    lambda      = 0.25,
    verbose     = FALSE
  )))
  result <- indivint(fit)

  # Predictions must reflect individual random effects, not a single constant.
  expect_gt(length(unique(round(result$amt_usual_orig, 6))), 50)
  expect_gt(stats::sd(result$predicted_usual_orig), 0)

  # The amount random effects should correlate with the true simulated v2.
  # (Subjects who consume more should have higher predicted amounts.)
  expect_gt(stats::sd(result$amt_usual_orig), 0)
})


test_that("indivint (amount model) equals the conditional expectation by integrate()", {
  dat <- make_amount_data_indivint(n_subjects = 200)
  fit <- suppressMessages(mixtran(dat, "intake", "subject", "day",
                                  model_type = "amount", verbose = FALSE))
  res <- indivint(fit)
  w <- fit$data; lam <- fit$lambda
  X <- cbind(1, as.matrix(w[, fit$amt_covs_used]))
  for (s in unique(w$subject)[c(1, 5, 9)]) {
    i <- w$subject == s
    z <- boxcox_transform(w$intake[i], lam); mu <- as.numeric(X[i, , drop = FALSE] %*% fit$beta)
    lik <- function(u) vapply(u, function(uu) prod(dnorm(z, mu + uu, sqrt(fit$sigma2_w))) *
                                dnorm(uu, 0, sqrt(fit$sigma2_b)), 1)
    x0 <- fit$predicted$x2b2[fit$predicted$subject == s]
    g  <- function(u) nciusual:::.backtransform_indivint(x0 + u, lam, fit$sigma2_w,
                                                          0.5 * fit$min_pos_raw)
    num <- integrate(function(u) lik(u) * g(u), -Inf, Inf, rel.tol = 1e-10)$value
    den <- integrate(lik, -Inf, Inf, rel.tol = 1e-10)$value
    expect_equal(res$predicted_usual_orig[res$subject == s], num / den, tolerance = 1e-6)
  }
})

test_that("indivint two-part mean approximates the DISTRIB mean", {
  dat <- make_episodic_data_indivint(n_subjects = 800)
  fit <- suppressMessages(mixtran(dat, "intake", "subject", "day",
                                  model_type = "corr", lambda = 0.25, verbose = FALSE))
  res <- indivint(fit)
  dm  <- distrib(fit, n_sims = 200, seed = 2)$results$Overall$mean
  expect_lt(abs(mean(res$predicted_usual_orig) - dm) / dm, 0.03)
})

test_that("indivint lamt option returns E[BoxCox(T, lamt)]", {
  dat <- make_amount_data_indivint(n_subjects = 150)
  fit <- suppressMessages(mixtran(dat, "intake", "subject", "day",
                                  model_type = "amount", verbose = FALSE))
  r1 <- indivint(fit); r2 <- indivint(fit, lamt = 0.5)
  # Jensen: E[BC(T)] <= BC(E[T]) for the concave transform
  expect_true(all(r2$indusint <= (r1$predicted_usual_orig^0.5 - 1) / 0.5 + 1e-10))
  expect_error(indivint(fit, lamt = 0), "non-zero")
})
