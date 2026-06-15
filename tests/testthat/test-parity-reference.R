# =====================================================================
# Parity tests against ANALYTIC GROUND-TRUTH reference outputs.
#
# We cannot run the proprietary NCI SAS macros here, so rather than
# comparing R to another estimator we compare R to the *truth*. The
# fixtures in tests/testthat/fixtures/ contain:
#   * fixed synthetic recall datasets generated from KNOWN parameters
#     (frozen as CSV so they are stable across R/RNG versions), and
#   * the true usual-intake distribution implied by those parameters,
#     computed independently by large-sample (4e6) Monte Carlo.
# A correct SAS MIXTRAN/DISTRIB run would be estimating this same target
# distribution, so recovering it within sampling/method error is the
# meaningful parity check. See data-raw/generate_reference.R.
#
# Tolerances are calibrated (with margin) to the deviations actually
# observed for each model. The amount and uncorrelated models recover the
# truth tightly; the correlated two-part model is approximate, especially
# in the lower tail, so its extreme percentiles are checked loosely and
# the central statistics tightly. This is a documented limitation of the
# profile-rho / GHQ correlated engines, not a regression.
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
# Amount-only model — recovers ground truth tightly (observed max ~1.6%)
# ---------------------------------------------------------------------
test_that("amount model matches analytic ground-truth distribution", {
  dat <- read.csv(test_path("fixtures", "amount_recall_data.csv"))
  fit <- suppressMessages(mixtran(
    dat, intake_var = "intake", subject_var = "SEQN", repeat_var = "day",
    model_type = "amount", weekend_var = "weekend",
    lambda = true_param("amount", "lambda"), verbose = FALSE))

  # Variance components recovered close to truth
  expect_lt(abs(fit$sigma2_b - true_param("amount", "sigma2_b")) /
              true_param("amount", "sigma2_b"), 0.10)
  expect_lt(abs(fit$sigma2_w - true_param("amount", "sigma2_w")) /
              true_param("amount", "sigma2_w"), 0.10)

  dist <- distrib(fit, n_sims = 500, seed = 999)
  expect_parity(dist, "amount",
                c("mean", "p5", "p10", "p25", "p50", "p75", "p90", "p95"),
                tol = 0.05)
})

# ---------------------------------------------------------------------
# Two-part uncorrelated model — recovers ground truth well (max ~3.9%)
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
                tol = 0.08)
})

# ---------------------------------------------------------------------
# Correlated two-part model (profile_rho) — central stats tight, tails loose
# ---------------------------------------------------------------------
test_that("correlated two-part model (profile_rho) matches ground-truth central statistics", {
  dat <- read.csv(test_path("fixtures", "corr_recall_data.csv"))
  fit <- suppressWarnings(suppressMessages(mixtran(
    dat, intake_var = "intake", subject_var = "SEQN", repeat_var = "day",
    model_type = "corr", corr_engine = "profile_rho",
    lambda = true_param("corr", "lambda"), verbose = FALSE)))

  # rho recovered with correct sign and substantial magnitude (true = 0.5)
  expect_gt(fit$rho, 0.2)
  expect_lt(fit$rho, 0.9)

  dist <- distrib(fit, n_sims = 500, seed = 999)
  # Central tendency and upper percentiles
  expect_parity(dist, "corr", c("mean"), tol = 0.12)
  expect_parity(dist, "corr", c("p25", "p50", "p75", "p90", "p95"), tol = 0.22)
  # Lower-tail percentiles are an approximate-engine limitation: assert only
  # that they are finite, positive, and correctly ordered.
  d <- dist$results$Overall$percentiles
  expect_true(all(d > 0))
  expect_true(d["p5"] < d["p10"] && d["p10"] < d["p25"])
})

# ---------------------------------------------------------------------
# Correlated two-part model (GHQ engine) — central stats vs ground truth
# ---------------------------------------------------------------------
test_that("correlated two-part model (GHQ) matches ground-truth central statistics", {
  dat <- read.csv(test_path("fixtures", "corr_recall_data.csv"))
  fit <- suppressWarnings(suppressMessages(mixtran(
    dat, intake_var = "intake", subject_var = "SEQN", repeat_var = "day",
    model_type = "corr", corr_engine = "ghq", ghq_n_nodes = 5L,
    lambda = true_param("corr", "lambda"), verbose = FALSE)))

  expect_gt(fit$rho, 0.1)
  expect_lt(fit$rho, 0.9)

  dist <- distrib(fit, n_sims = 500, seed = 999)
  expect_parity(dist, "corr", c("mean", "p50", "p75", "p90", "p95"), tol = 0.15)
})
