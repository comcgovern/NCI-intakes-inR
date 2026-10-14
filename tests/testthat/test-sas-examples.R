# Parity with NCI's published SAS output for MIXTRAN/DISTRIB v2.1 Examples 1
# (amount model) and 2 (correlated two-part model), run on NCI's own
# analytic dataset.
#
# These tests download NCI's ~14 MB example dataset, so they only run when
# NCIUSUAL_SAS_PARITY=true (optionally set NCIUSUAL_NCI_DATA to a local copy
# of indiv_ims.sas7bdat). Run with:
#   NCIUSUAL_SAS_PARITY=true Rscript -e 'testthat::test_file("tests/testthat/test-sas-examples.R")'
#
# Tolerances. MIXTRAN: SAS's NLMIXED stops at a relative-gradient criterion,
# so its estimates sit slightly short of the maximum (our -2 log likelihood is
# lower than SAS's at SAS's own estimates); coefficients agree to ~0.02.
# DISTRIB: SAS's numbers are one Monte Carlo run (100 sims/subject), so
# percentiles (3-4%) and cutpoint proportions (0.02; a single run in a
# subgroup of ~1000 subjects carries about +-0.015 of simulation error) are
# compared with tolerances covering that error.

skip_sas <- function() {
  skip_if_not(nci_examples_enabled(), "set NCIUSUAL_SAS_PARITY=true to run")
  skip_if_not_installed("haven")
}

sas_val <- function(ref, ex, kind, stat, group = "") {
  v <- ref$value[ref$example == ex & ref$kind == kind & ref$statistic == stat &
                   ref$group == group]
  stopifnot(length(v) == 1L)
  v
}

check_distrib <- function(res, ref, ex, ncut, pct_tol, cut_tol) {
  for (g in names(res$results)) {
    r <- res$results[[g]]
    expect_lt(abs(r$mean / sas_val(ref, ex, "distrib", "mean", g) - 1), pct_tol,
              label = sprintf("Example %d %s mean", ex, g))
    for (p in c(5, 10, 25, 50, 75, 90, 95)) {
      est <- r$percentiles[[paste0("p", p)]]
      expect_lt(abs(est / sas_val(ref, ex, "distrib", paste0("p", p), g) - 1), pct_tol,
                label = sprintf("Example %d %s p%d", ex, g, p))
    }
    for (k in seq_len(ncut)) {
      expect_lt(abs(r$cutpoint_below[[k]] - sas_val(ref, ex, "distrib", paste0("cutprob", k), g)),
                cut_tol, label = sprintf("Example %d %s cutprob%d", ex, g, k))
    }
  }
}

test_that("Example 1 (amount model): MIXTRAN estimates match SAS", {
  skip_sas()
  ref <- sas_reference()
  fit <- suppressMessages(nci_example_fit(1))
  expect_true(fit$converged)
  for (nm in c("(Intercept)", "age4", "age5", "age6", "age7", "age8", "race1",
               "race3", "weekend", "seq_num")) {
    expect_lt(abs(fit$beta[[nm]] - sas_val(ref, 1, "mixtran_amount", paste0("A.", nm))), 0.025,
              label = paste("beta", nm))
  }
  expect_lt(abs(fit$lambda - sas_val(ref, 1, "mixtran_amount", "lambda")), 0.002)
  expect_lt(abs(fit$sigma2_b / sas_val(ref, 1, "mixtran_amount", "var_u2") - 1), 0.005)
  expect_lt(abs(fit$sigma2_w / sas_val(ref, 1, "mixtran_amount", "var_e") - 1), 0.005)
  # SAS prints -2LL to 5 significant digits
  expect_lt(abs(fit$m2loglik / sas_val(ref, 1, "mixtran_amount", "m2loglik") - 1), 1e-4)
})

test_that("Example 1 (amount model): DISTRIB matches SAS", {
  skip_sas()
  ref <- sas_reference()
  fit <- suppressMessages(nci_example_fit(1))
  res <- distrib(fit, subgroup_var = "agegrp", cutpoints = nci_example_cutpoints[[1]],
                 percentiles = c(5, 10, 25, 50, 75, 90, 95), n_sims = 100, seed = 5454768)
  check_distrib(res, ref, 1, 10, pct_tol = 0.03, cut_tol = 0.02)
})

test_that("Example 2 (correlated model): MIXTRAN estimates match SAS", {
  skip_sas()
  ref <- sas_reference()
  fit <- suppressMessages(nci_example_fit(2))
  expect_true(fit$converged)
  # Uncorrelated stage (MIXTRAN nocorr NLMIXED fits)
  th <- fit$theta_uncorr
  for (nm in setdiff(ref$statistic[ref$kind == "mixtran_nocorr"], "m2loglik")) {
    expect_lt(abs(th[[nm]] - sas_val(ref, 2, "mixtran_nocorr", nm)), 0.003,
              label = paste("nocorr", nm))
  }
  expect_lt(abs(fit$m2loglik_uncorr / sas_val(ref, 2, "mixtran_nocorr", "m2loglik") - 1), 1e-5)
  # Correlated model: the likelihood is flat along (var_u1, rho), so SAS's
  # stopping point differs a little more; our -2LL is lower than SAS's
  # (checked in data-raw/sas_parity_examples.R).
  for (nm in c("(Intercept)", "sex1", "age2", "race1", "race3", "weekend", "seq_num")) {
    expect_lt(abs(fit$alpha[[nm]] - sas_val(ref, 2, "mixtran_corr", paste0("P.", nm))), 0.03)
    expect_lt(abs(fit$beta[[nm]]  - sas_val(ref, 2, "mixtran_corr", paste0("A.", nm))), 0.01)
  }
  expect_lt(abs(fit$lambda - sas_val(ref, 2, "mixtran_corr", "A_LAMBDA")), 0.002)
  expect_lt(abs(fit$rho - sas_val(ref, 2, "mixtran_corr", "rho")), 0.03)
  expect_lt(abs(fit$sigma2_v2 / sas_val(ref, 2, "mixtran_corr", "var_u2") - 1), 0.01)
  expect_lt(abs(fit$sigma2_e / sas_val(ref, 2, "mixtran_corr", "var_e") - 1), 0.005)
  expect_lt(abs(fit$sigma2_v1 / exp(2 * sas_val(ref, 2, "mixtran_corr", "P_LOGSDU1")) - 1), 0.06)
  expect_lt(abs(fit$m2loglik / sas_val(ref, 2, "mixtran_corr", "m2loglik") - 1), 1e-3)
})

test_that("Example 2 (correlated model): DISTRIB matches SAS", {
  skip_sas()
  ref <- sas_reference()
  fit <- suppressMessages(nci_example_fit(2))
  res <- distrib(fit, subgroup_var = "agegrp", cutpoints = nci_example_cutpoints[[2]],
                 percentiles = c(5, 10, 25, 50, 75, 90, 95), n_sims = 100, seed = 5454768)
  check_distrib(res, ref, 2, 3, pct_tol = 0.04, cut_tol = 0.02)
})
