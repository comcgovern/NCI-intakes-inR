# =============================================================================
# SAS parity report: nciusual vs NCI's published MIXTRAN/DISTRIB v2.1 output
# =============================================================================
#
# Runs nciusual on NCI's own analytic dataset for the macro example programs
#   Example 1: added sugars, males 9+, amount-only model
#   Example 2: total fruit, children 1-8, correlated two-part model
# (covariates, seq = day2, weekend, replicate_var = rndw1, DISTRIB with
# nsim_mc = 100, seed = 5454768, wkend_prop = 3/7, subgroup = agegrp) and
# prints the estimates side by side with SAS's listing output, which is frozen
# in tests/testthat/fixtures/sas_examples_reference.csv. For each model it
# also evaluates nciusual's likelihood at SAS's estimates, showing whether
# remaining parameter differences come from where SAS's optimiser stopped.
#
# Sources:
#   https://epi.grants.cancer.gov/diet/usualintakes/macros_single.html
#   https://epi.grants.cancer.gov/diet/usualintakes/dataset.html
#
# Run from the package root (downloads ~14 MB unless NCIUSUAL_NCI_DATA points
# to a local indiv_ims.sas7bdat):
#   Rscript data-raw/sas_parity_examples.R
# =============================================================================

suppressMessages(library(nciusual))
test_path <- function(...) file.path("tests", "testthat", ...)
source(test_path("helper-nci-examples.R"))
ref <- sas_reference()
sv <- function(ex, kind, stat, group = "") {
  ref$value[ref$example == ex & ref$kind == kind & ref$statistic == stat & ref$group == group]
}

m2ll_at <- function(fit, pars, part) {
  d <- nciusual:::.nci_design(fit$data, fit$prob_covs_used, fit$amt_covs_used,
                              NA_real_, fit$n_quad, has_prob = fit$model_type != "amount")
  2 * d$sum_w * nciusual:::.neg_wll(pars, d, part)
}

distrib_table <- function(fit, ex) {
  res <- distrib(fit, subgroup_var = "agegrp", cutpoints = nci_example_cutpoints[[ex]],
                 percentiles = c(1, 5, 10, 25, 50, 75, 90, 95, 99),
                 n_sims = 100, seed = 5454768)
  out <- lapply(names(res$results), function(g) {
    r <- res$results[[g]]
    stats <- c("mean", paste0("p", c(1, 5, 10, 25, 50, 75, 90, 95, 99)))
    est <- c(r$mean, r$percentiles)
    sas <- vapply(stats, function(s) sv(ex, "distrib", s, g), 1)
    cp <- r$cutpoint_below
    cps <- vapply(seq_along(cp), function(k) sv(ex, "distrib", paste0("cutprob", k), g), 1)
    data.frame(group = g, statistic = c(stats, names(cp)),
               nciusual = round(c(est, cp), 5), SAS = round(c(sas, cps), 5),
               diff = c(sprintf("%+.1f%%", 100 * (est / sas - 1)),
                        sprintf("%+.4f", cp - cps)))
  })
  do.call(rbind, out)
}

# ---- Example 1 ---------------------------------------------------------------
cat("\n=== Example 1: amount model (added sugars, males 9+) ===\n")
f1 <- suppressMessages(nci_example_fit(1))
k  <- ref$kind == "mixtran_amount" & ref$example == 1
nm <- setdiff(ref$statistic[k], "m2loglik")
ours <- c(f1$beta[sub("^A\\.", "", grep("^A\\.", nm, value = TRUE))],
          lambda = f1$lambda, var_u2 = f1$sigma2_b, var_e = f1$sigma2_w)
print(data.frame(parameter = nm, nciusual = round(unname(ours), 4),
                 SAS = vapply(nm, function(s) sv(1, "mixtran_amount", s), 1), row.names = NULL))
sas_th <- f1$theta
for (s in grep("^A\\.", nm, value = TRUE)) sas_th[s] <- sv(1, "mixtran_amount", s)
sas_th["A_LAMBDA"] <- sv(1, "mixtran_amount", "lambda")
sas_th["A_LOGSDU2"] <- 0.5 * log(sv(1, "mixtran_amount", "var_u2"))
sas_th["A_LOGSDE"] <- 0.5 * log(sv(1, "mixtran_amount", "var_e"))
cat(sprintf("-2 log L: nciusual %.1f | at SAS estimates %.1f | SAS printed %.5g\n",
            f1$m2loglik, m2ll_at(f1, sas_th, "amount"), sv(1, "mixtran_amount", "m2loglik")))
print(distrib_table(f1, 1), row.names = FALSE)

# ---- Example 2 ---------------------------------------------------------------
cat("\n=== Example 2: correlated two-part model (total fruit, children 1-8) ===\n")
f2 <- suppressMessages(nci_example_fit(2))
for (stage in c("nocorr", "corr")) {
  kind <- paste0("mixtran_", stage)
  nm <- intersect(setdiff(ref$statistic[ref$kind == kind], "m2loglik"),
                  names(f2$theta))
  th <- if (stage == "nocorr") f2$theta_uncorr else f2$theta
  cat(sprintf("\n-- %s stage --\n", stage))
  print(data.frame(parameter = nm, nciusual = round(unname(th[nm]), 4),
                   SAS = vapply(nm, function(s) sv(2, kind, s), 1), row.names = NULL))
}
cat(sprintf("rho: nciusual %.4f | SAS %.4f\n", f2$rho, sv(2, "mixtran_corr", "rho")))
sas_th <- f2$theta
for (s in intersect(names(sas_th), ref$statistic[ref$kind == "mixtran_corr"])) {
  sas_th[s] <- sv(2, "mixtran_corr", s)
}
sas_th["Z_U"] <- atanh(sv(2, "mixtran_corr", "rho"))
cat(sprintf("-2 log L (corr): nciusual %.1f | at SAS estimates %.1f | SAS printed %.4g\n",
            f2$m2loglik, m2ll_at(f2, sas_th, "corr"), sv(2, "mixtran_corr", "m2loglik")))
print(distrib_table(f2, 2), row.names = FALSE)
