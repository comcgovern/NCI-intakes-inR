#' Bivariate Usual Intake and Ratio Estimation
#'
#' Functions for estimating the joint distribution of two dietary components
#' and derived statistics (e.g., the ratio of usual intake of one nutrient to
#' another, or to total energy).
#'
#' The NCI multivariate SAS macros fit a joint mixed model with correlated
#' random effects across dietary components. This function is NOT a port of
#' those macros: it fits each component with [mixtran()], approximates the
#' between-component correlation from the empirical Bayes estimates, and
#' propagates it through a DISTRIB-style joint Monte Carlo simulation.
#'
#' @name bivariate
NULL


#' Estimate the bivariate usual intake distribution for two dietary components
#'
#' Fits separate `mixtran()` models for two dietary components, estimates the
#' cross-component correlation between their random effects, then jointly
#' simulates usual intakes for both components using correlated draws. This
#' allows computing distribution statistics for linear combinations or ratios
#' (e.g., percent energy from fat = fat_kcal / total_energy × 100).
#'
#' @param data Data frame in long format (one row per recall per person).
#' @param intake_var1 Name of the first intake variable (character).
#' @param intake_var2 Name of the second intake variable (character).
#' @param subject_var Name of the subject ID variable.
#' @param repeat_var Name of the recall sequence variable.
#' @param model_type1,model_type2 Model type for each component: `"amount"`,
#'   `"uncorr"`, or `"corr"`. Defaults to `"amount"`.
#' @param covariates Optional character vector of covariate names shared by
#'   both models. Per-component covariates are not currently supported; pass
#'   `NULL` and handle complex cases via separate `mixtran()` calls.
#' @param weekend_var,weight_var,lambda1,lambda2,verbose Passed through to
#'   `mixtran()` for each component.
#' @param n_sims Number of Monte Carlo replications (default 100).
#' @param seed Random seed.
#' @param percentiles Percentiles to estimate for each component and the ratio.
#' @param cutpoints Optional numeric vector of cutpoints for proportion
#'   below/above.
#' @param ratio_label Label used in output for the derived ratio
#'   (e.g., `"pct_energy_fat"`). Default `"ratio"`.
#' @param wkend_prop Weekend weight passed to the DISTRIB-style simulation.
#' @param fn_ratio Function of (usual_1, usual_2) that defines the derived
#'   quantity.  Default `function(u1, u2) u1 / u2` (simple ratio). Could also
#'   be `function(u1, u2) u1 / u2 * 100` for percentage, or
#'   `function(u1, u2) u1 + u2` for a sum, etc.
#' @return A `bivariate_distrib` object (a list) with:
#'   \describe{
#'     \item{component1}{`distrib_result`-style summary for intake_var1}
#'     \item{component2}{`distrib_result`-style summary for intake_var2}
#'     \item{ratio}{`distrib_result`-style summary for the derived quantity}
#'     \item{cross_rho}{Estimated between-component random-effect correlation}
#'     \item{fit1,fit2}{`mixtran_fit` objects for each component}
#'   }
#' @export
distrib_bivariate <- function(data,
                               intake_var1,
                               intake_var2,
                               subject_var,
                               repeat_var,
                               model_type1 = "amount",
                               model_type2 = "amount",
                               covariates  = NULL,
                               weekend_var = NULL,
                               weight_var  = NULL,
                               lambda1     = NULL,
                               lambda2     = NULL,
                               verbose     = TRUE,
                               n_sims      = 100L,
                               seed        = 12345L,
                               percentiles = c(5, 10, 25, 50, 75, 90, 95),
                               cutpoints   = NULL,
                               ratio_label = "ratio",
                               fn_ratio    = function(u1, u2) u1 / u2,
                               wkend_prop  = 3 / 7) {

  data <- as.data.frame(data)

  # --- Fit each component separately ---
  if (verbose) message("Fitting model for component 1: ", intake_var1)
  fit1 <- mixtran(
    data        = data,
    intake_var  = intake_var1,
    subject_var = subject_var,
    repeat_var  = repeat_var,
    model_type  = model_type1,
    covariates  = covariates,
    weekend_var = weekend_var,
    weight_var  = weight_var,
    lambda      = lambda1,
    verbose     = verbose
  )

  if (verbose) message("Fitting model for component 2: ", intake_var2)
  fit2 <- mixtran(
    data        = data,
    intake_var  = intake_var2,
    subject_var = subject_var,
    repeat_var  = repeat_var,
    model_type  = model_type2,
    covariates  = covariates,
    weekend_var = weekend_var,
    weight_var  = weight_var,
    lambda      = lambda2,
    verbose     = verbose
  )

  # --- Estimate cross-component correlation from random effects ---
  cross_rho <- .estimate_cross_rho(fit1, fit2, verbose)
  if (verbose) message(sprintf("  Cross-component RE correlation: %.4f", cross_rho))

  # --- Align subjects (intersection, in the same order) ---
  pred1 <- fit1$predicted
  pred2 <- fit2$predicted
  common_subj <- intersect(pred1$subject, pred2$subject)
  pred1 <- pred1[match(common_subj, pred1$subject), , drop = FALSE]
  pred2 <- pred2[match(common_subj, pred2$subject), , drop = FALSE]
  n_subj <- length(common_subj)
  weights_use <- pred1$weight  # use component-1 weights (should be identical)

  # --- Joint Monte Carlo simulation ---
  if (verbose) message("Running joint Monte Carlo simulation...")

  sim_result <- .with_seed(seed, .simulate_bivariate(
    subj1      = pred1,
    subj2      = pred2,
    fit1       = fit1,
    fit2       = fit2,
    cross_rho  = cross_rho,
    n_sims     = n_sims,
    fn_ratio   = fn_ratio,
    wkend_prop = wkend_prop
  ))

  # --- Compute distribution summaries ---
  wts_expanded <- rep(weights_use, each = n_sims)

  sum1  <- compute_distribution_summary(sim_result$ui1,   wts_expanded, percentiles, cutpoints)
  sum2  <- compute_distribution_summary(sim_result$ui2,   wts_expanded, percentiles, cutpoints)
  sum_r <- compute_distribution_summary(sim_result$ratio, wts_expanded, percentiles,
                                         cutpoints = NULL)  # cutpoints on ratio rarely used

  out <- list(
    component1  = list(Overall = sum1),
    component2  = list(Overall = sum2),
    ratio       = list(Overall = sum_r),
    cross_rho   = cross_rho,
    fit1        = fit1,
    fit2        = fit2,
    intake_var1 = intake_var1,
    intake_var2 = intake_var2,
    ratio_label = ratio_label,
    n_sims      = n_sims,
    seed        = seed
  )
  class(out) <- "bivariate_distrib"
  out
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Estimate cross-component correlation between random effects
#'
#' Uses the empirical Bayes estimates of the amount random effects from both
#' fitted models and computes their Pearson correlation. This is an
#' approximation: EBEs are shrunken, and the NCI multivariate macros instead
#' estimate the cross-component correlation within a joint model.
#'
#' @keywords internal
.estimate_cross_rho <- function(fit1, fit2, verbose = FALSE) {

  re1 <- .extract_re(fit1)
  re2 <- .extract_re(fit2)

  # Match on subject IDs
  common <- intersect(names(re1), names(re2))
  if (length(common) < 10) {
    if (verbose) warning("Fewer than 10 subjects for cross-rho estimation; ",
                         "using rho = 0.")
    return(0)
  }

  r1 <- re1[common]
  r2 <- re2[common]
  stats::cor(r1, r2, use = "complete.obs")
}


#' Extract the amount random-effect EBEs (named by subject ID as character)
#' @keywords internal
.extract_re <- function(fit) {
  stats::setNames(fit$predicted$u2, as.character(fit$predicted$subject))
}


#' Joint simulation for two dietary components with correlated random effects
#'
#' @keywords internal
.simulate_bivariate <- function(subj1, subj2, fit1, fit2,
                                 cross_rho, n_sims, fn_ratio,
                                 wkend_prop = 3 / 7) {

  n_subj   <- nrow(subj1)
  rho_safe <- min(max(cross_rho, -0.9999), 0.9999)

  # Standardised amount random effects of the two components, correlated
  # through cross_rho; each component's probability random effect is then
  # correlated with its own amount random effect (DISTRIB-style draws).
  w1 <- matrix(stats::rnorm(n_subj * n_sims), n_subj, n_sims)
  w2 <- rho_safe * w1 +
    sqrt(1 - rho_safe^2) * matrix(stats::rnorm(n_subj * n_sims), n_subj, n_sims)
  e1 <- matrix(stats::rnorm(n_subj * n_sims), n_subj, n_sims)
  e2 <- matrix(stats::rnorm(n_subj * n_sims), n_subj, n_sims)

  ui1 <- as.vector(t(.simulate_usual(subj1, fit1, e1, w1, wkend_prop)$t))
  ui2 <- as.vector(t(.simulate_usual(subj2, fit2, e2, w2, wkend_prop)$t))

  ratio <- fn_ratio(ui1, ui2)
  ratio[!is.finite(ratio)] <- NA_real_

  list(ui1 = ui1, ui2 = ui2, ratio = ratio)
}


# ---------------------------------------------------------------------------
# Print / summary / plot for bivariate_distrib
# ---------------------------------------------------------------------------

#' Print a bivariate_distrib object
#' @param x A `bivariate_distrib` object from distrib_bivariate()
#' @param ... Unused
#' @export
print.bivariate_distrib <- function(x, ...) {
  cat("NCI Bivariate Usual Intake Distribution\n")
  cat(sprintf("  Component 1 : %s\n", x$intake_var1))
  cat(sprintf("  Component 2 : %s\n", x$intake_var2))
  cat(sprintf("  Derived     : %s  (fn_ratio applied)\n", x$ratio_label))
  cat(sprintf("  Cross-rho   : %.4f\n", x$cross_rho))
  cat(sprintf("  N sims      : %d\n", x$n_sims))
  cat("\n")

  .print_biv_section <- function(label, res_list) {
    r <- res_list$Overall
    cat(sprintf("  %s — Mean (SD): %.2f (%.2f)\n", label, r$mean, r$sd))
    cat(sprintf("    Percentiles: %s\n",
                paste(sprintf("%s=%.2f", names(r$percentiles), r$percentiles),
                      collapse = "  ")))
  }

  .print_biv_section(x$intake_var1, x$component1)
  .print_biv_section(x$intake_var2, x$component2)
  .print_biv_section(x$ratio_label, x$ratio)
  invisible(x)
}


#' Summarise a bivariate_distrib object
#' @param object A `bivariate_distrib` object from distrib_bivariate()
#' @param ... Unused
#' @export
summary.bivariate_distrib <- function(object, ...) {
  print(object, ...)
  invisible(object)
}


#' Extract bivariate distribution results as a tidy data frame
#'
#' @param biv_obj A `bivariate_distrib` object from distrib_bivariate()
#' @return Data frame with one row per (component, statistic).
#' @export
bivariate_to_df <- function(biv_obj) {
  comps <- list(
    biv_obj$component1$Overall,
    biv_obj$component2$Overall,
    biv_obj$ratio$Overall
  )
  labels <- c(biv_obj$intake_var1, biv_obj$intake_var2, biv_obj$ratio_label)

  rows <- mapply(function(r, lbl) {
    row <- data.frame(component = lbl, mean = r$mean, sd = r$sd,
                      stringsAsFactors = FALSE)
    for (nm in names(r$percentiles)) row[[nm]] <- r$percentiles[[nm]]
    if (!is.null(r$cutpoint_below)) {
      for (nm in names(r$cutpoint_below)) row[[nm]] <- r$cutpoint_below[[nm]]
    }
    row
  }, comps, labels, SIMPLIFY = FALSE)

  do.call(rbind, rows)
}
