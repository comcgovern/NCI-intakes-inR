# =====================================================================
# Parity tests against ANALYTIC GROUND-TRUTH reference outputs.
#
# The fixtures in tests/testthat/fixtures/ contain fixed synthetic recall
# datasets generated from KNOWN parameters and the true usual-intake
# distribution they imply under the NCI definition of usual intake (the
# expected daily intake: the back-transformation integrated over
# within-person error, with weekday/weekend averaging for the amount model),
# computed by large-sample Monte Carlo. See data-raw/generate_reference.R.
#
# Direct parity with SAS output (NCI's published Example 1 and 2 listings)
# is checked in test-sas-examples.R and data-raw/sas_parity_example1.R.
# =====================================================================

ref_outputs <- read.csv(test_path("fixtures", "reference_outputs.csv"),
                        stringsAsFactors = FALSE)
ref_params  <- read.csv(test_path("fixtures", "reference_params.csv"),
                        stringsAsFactors = FALSE)

ref_value <- function(model, statistic) {
  v <- ref_outputs$value[ref_outputs$model == model &
                           ref_outputs$statistic == statistic]
  stopifnot(length(v) == 1L)
  v
}
true_param <- function(model, parameter) {
  v <- ref_params$value[ref_params$model == model &
                          ref_params$parameter == parameter]
  stopifnot(length(v) == 1L)
  v
}

# Compare a fitted distrib() result to the reference for a set of statistics.
expect_parity <- function(dist, model, statistics, tol) {
  d <- dist$results$Overall
  est_lookup <- c(list(mean = d$mean), as.list(d$percentiles))
  for (s in statistics) {
    est <- as.numeric(est_lookup[[s]])
    tru <- ref_value(model, s)
    rel <- abs(est - tru) / abs(tru)
    expect_lt(rel, tol,
      label = sprintf("%s %s: est=%.4f true=%.4f rel_err=%.1f%% (tol=%.0f%%)",
                      model, s, est, tru, rel * 100, tol * 100))
  }
}

# ---------------------------------------------------------------------
# Amount-only model
# ---------------------------------------------------------------------
test_that("amount model matches analytic ground-truth distribution", {
  dat <- read.csv(test_path("fixtures", "amount_recall_data.csv"))
  fit <- suppressMessages(mixtran(
    dat, intake_var = "intake", subject_var = "SEQN", repeat_var = "day",
    model_type = "amount", weekend_var = "weekend", verbose = FALSE))

  # lambda is estimated by maximum likelihood (true 0.30; the MLE for this
  # fixed dataset is about 0.25 -- the likelihood is flat in lambda)
  expect_true(fit$converged)
  expect_lt(abs(fit$lambda - true_param("amount", "lambda")), 0.08)

  # With lambda fixed at the truth the variance components are on the
  # generating scale and are recovered
  fit_l <- suppressMessages(mixtran(
    dat, intake_var = "intake", subject_var = "SEQN", repeat_var = "day",
    model_type = "amount", weekend_var = "weekend",
    lambda = true_param("amount", "lambda"), verbose = FALSE))
  expect_lt(abs(fit_l$sigma2_b - true_param("amount", "sigma2_b")) /
              true_param("amount", "sigma2_b"), 0.10)
  expect_lt(abs(fit_l$sigma2_w - true_param("amount", "sigma2_w")) /
              true_param("amount", "sigma2_w"), 0.10)

  dist <- distrib(fit, n_sims = 500, seed = 999)
  expect_parity(dist, "amount",
                c("mean", "p5", "p10", "p25", "p50", "p75", "p90", "p95"),
                tol = 0.04)
})

# ---------------------------------------------------------------------
# Two-part uncorrelated model
# ---------------------------------------------------------------------
test_that("uncorrelated two-part model matches analytic ground-truth distribution", {
  dat <- read.csv(test_path("fixtures", "uncorr_recall_data.csv"))
  fit <- suppressMessages(mixtran(
    dat, intake_var = "intake", subject_var = "SEQN", repeat_var = "day",
    model_type = "uncorr", lambda = true_param("uncorr", "lambda"),
    verbose = FALSE))

  for (vc in c("sigma2_v1", "sigma2_v2", "sigma2_e")) {
    expect_lt(abs(fit[[vc]] - true_param("uncorr", vc)) /
                true_param("uncorr", vc), 0.15)
  }

  dist <- distrib(fit, n_sims = 500, seed = 999)
  expect_parity(dist, "uncorr",
                c("mean", "p5", "p10", "p25", "p50", "p75", "p90", "p95"),
                tol = 0.06)
})

# ---------------------------------------------------------------------
# Two-part correlated model: joint maximum likelihood recovers rho and the
# whole distribution, including the lower tail.
# ---------------------------------------------------------------------
test_that("correlated two-part model matches the full ground-truth distribution", {
  dat <- read.csv(test_path("fixtures", "corr_recall_data.csv"))
  fit <- suppressMessages(mixtran(
    dat, intake_var = "intake", subject_var = "SEQN", repeat_var = "day",
    model_type = "corr", lambda = true_param("corr", "lambda"),
    verbose = FALSE))

  expect_lt(abs(fit$rho - true_param("corr", "rho")), 0.15)
  expect_lt(abs(as.numeric(fit$beta["(Intercept)"]) - true_param("corr", "beta0")), 0.1)
  for (vc in c("sigma2_v1", "sigma2_v2", "sigma2_e")) {
    expect_lt(abs(fit[[vc]] - true_param("corr", vc)) /
                true_param("corr", vc), 0.2)
  }

  dist <- distrib(fit, n_sims = 500, seed = 999)
  expect_parity(dist, "corr",
                c("mean", "p5", "p10", "p25", "p50", "p75", "p90", "p95"),
                tol = 0.07)
})
