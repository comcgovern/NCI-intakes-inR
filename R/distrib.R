#' DISTRIB: Estimate Usual Intake Distributions
#'
#' R port of the NCI DISTRIB SAS macro (v2.1). Uses the parameter estimates
#' and subject-level predictions from [mixtran()] and Monte Carlo simulation of
#' the random effects to estimate the distribution of usual intake for a
#' population or subpopulations.
#'
#' @name distrib
NULL

#' Estimate usual intake distribution (DISTRIB equivalent)
#'
#' For each subject and each of `n_sims` repetitions, DISTRIB draws the random
#' effects `(u1, u2)` from their fitted (bivariate) normal distribution and
#' computes usual intake on the original scale as in the DISTRIB macro:
#'
#' * probability of consumption `p = logistic(x1b1 + u1)` (1 for amount-only
#'   models);
#' * consumption-day amount `a = E_e[T^-1(x2b2 + u2 + e)]`, the expectation
#'   over the within-person error `e ~ N(0, var_e)`, computed with DISTRIB's
#'   9-point approximation (exact `exp(x2b2 + u2 + var_e/2)` when lambda = 0)
#'   and bounded below by half of the minimum amount;
#' * usual intake `t = p * a`.
#'
#' The sequence effect is removed (first-recall reference). When the model
#' includes a weekend indicator, `t`, `p` and `a` are computed for a weekday
#' and a weekend day and averaged with weights `1 - wkend_prop` and
#' `wkend_prop` (DISTRIB default 4/7 and 3/7). Covariates are kept at each
#' subject's observed (first-record) values.
#'
#' Percentiles and cutpoint proportions use the DISTRIB PC-macro algorithm:
#' weights are normalised within each group, tied values are pooled,
#' cumulative proportions are taken at the midpoint of each value's weight,
#' and values are interpolated (and extrapolated at the ends) linearly.
#'
#' @param mixtran_obj Output from [mixtran()].
#' @param subgroup_var Optional name of a subgroup variable. It must be in
#'   `mixtran_obj$predicted` (pass `subgroup_var` to `mixtran()`) or in
#'   `subgroup_data`.
#' @param cutpoints Numeric vector of cutpoints; the weighted proportion of
#'   usual intake below (and above) each cutpoint is reported.
#' @param percentiles Percentiles to estimate (0-100).
#' @param n_sims Number of Monte Carlo repetitions per subject (`nsim_mc`).
#' @param seed Random seed for reproducibility. The caller's RNG state is
#'   restored on exit.
#' @param wkend_prop Weight given to weekend days when the model includes a
#'   weekend indicator (DISTRIB `wkend_prop`; default 3/7).
#' @param subgroup_data Optional data frame with the subject ID column
#'   (`mixtran_obj$subject_var`) and `subgroup_var`, like DISTRIB's `add_da`.
#' @return A `distrib_result` object.
#' @export
distrib <- function(mixtran_obj,
                    subgroup_var = NULL,
                    cutpoints = NULL,
                    percentiles = c(5, 10, 25, 50, 75, 90, 95),
                    n_sims = 100,
                    seed = 12345,
                    wkend_prop = 3 / 7,
                    subgroup_data = NULL) {

  stopifnot(inherits(mixtran_obj, "mixtran_fit"))
  if (!is.numeric(wkend_prop) || wkend_prop <= 0 || wkend_prop >= 1) {
    stop("'wkend_prop' must be strictly between 0 and 1.")
  }
  .with_seed(seed, {
    pred <- mixtran_obj$predicted
    if (!is.null(subgroup_var)) {
      pred <- .attach_subgroup(pred, mixtran_obj, subgroup_var, subgroup_data)
    }
    n_subj <- nrow(pred)
    z1 <- matrix(stats::rnorm(n_subj * n_sims), n_subj, n_sims)
    z2 <- matrix(stats::rnorm(n_subj * n_sims), n_subj, n_sims)
    sims <- .simulate_usual(pred, mixtran_obj, z1, z2, wkend_prop)
  })

  groups <- list(Overall = seq_len(n_subj))
  if (!is.null(subgroup_var)) {
    g <- pred[[subgroup_var]]
    lev <- if (is.factor(g)) levels(droplevels(g)) else sort(unique(g[!is.na(g)]))
    for (lv in lev) groups[[as.character(lv)]] <- which(!is.na(g) & g == lv)
    groups <- c(groups[-1], groups[1])   # subgroups first, then Overall
  }

  results <- lapply(groups, function(idx) {
    compute_distribution_summary(
      as.numeric(t(sims$t[idx, , drop = FALSE])),
      rep(pred$weight[idx], each = n_sims),
      percentiles, cutpoints
    )
  })

  output <- list(
    results = results,
    simulated = list(usual_intake = as.vector(t(sims$t)),
                     prob = as.vector(t(sims$p)),
                     amount = as.vector(t(sims$a)),
                     n_subjects = n_subj, n_sims = n_sims),
    mixtran_obj = mixtran_obj,
    n_sims = n_sims,
    seed = seed,
    wkend_prop = wkend_prop,
    percentiles_requested = percentiles,
    cutpoints = cutpoints,
    subgroup_var = subgroup_var
  )
  class(output) <- "distrib_result"
  output
}


# Evaluate `expr` with a given seed and restore the caller's RNG state.
.with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(".Random.seed", envir = globalenv())
    } else assign(".Random.seed", old, envir = globalenv())
  })
  set.seed(seed)
  eval.parent(substitute(expr))
}


.attach_subgroup <- function(pred, fit, subgroup_var, subgroup_data) {
  if (subgroup_var %in% names(pred)) return(pred)
  if (is.null(subgroup_data)) {
    stop("Subgroup variable '", subgroup_var, "' not found in predicted data. ",
         "Pass subgroup_var to mixtran() or supply subgroup_data.")
  }
  sid <- fit$subject_var
  if (!all(c(sid, subgroup_var) %in% names(subgroup_data))) {
    stop("subgroup_data must contain '", sid, "' and '", subgroup_var, "'.")
  }
  sd <- subgroup_data[!duplicated(subgroup_data[[sid]]), , drop = FALSE]
  pred[[subgroup_var]] <- sd[[subgroup_var]][match(pred$subject, sd[[sid]])]
  pred
}


#' DISTRIB back-transformation of the consumption-day amount
#'
#' Expected original-scale amount given the transformed-scale mean `mu`,
#' integrating over the within-person error `e ~ N(0, var_e)`: exact
#' (`exp(mu + var_e/2)`) when `lambda = 0`, otherwise DISTRIB's 9-point
#' approximation. The result is bounded below by `0.5 * min_amt`.
#'
#' @param mu Numeric vector or matrix of transformed-scale means.
#' @param lambda Box-Cox parameter.
#' @param var_e Within-person variance on the transformed scale.
#' @param min_amt Minimum amount (see `mixtran_fit$min_amt`).
#' @return Back-transformed amounts with the same shape as `mu`.
#' @keywords internal
backtransform_9pt <- function(mu, lambda, var_e, min_amt) {
  floor_amt <- 0.5 * min_amt
  if (abs(lambda) < 1e-10) {
    out <- exp(mu + var_e / 2)
  } else {
    cj <- c(-2.1, -1.3, -0.8, -0.5, 0, 0.5, 0.8, 1.3, 2.1)
    wj <- c(0.063345, 0.080255, 0.070458, 0.159698, 0.252489,
            0.159698, 0.070458, 0.080255, 0.063345)
    sd_e <- sqrt(var_e)
    out <- 0
    for (j in seq_along(cj)) {
      inner <- (mu + cj[j] * sd_e) * lambda + 1
      out <- out + wj[j] * pmax(inner, 0)^(1 / lambda)
    }
  }
  pmax(out, floor_amt)
}


#' Simulate usual intake for each subject (DISTRIB macro MC)
#'
#' @param pred Subject-level predicted data from `mixtran()`.
#' @param fit A `mixtran_fit`.
#' @param z1,z2 Independent standard-normal matrices (n_subjects x n_sims);
#'   `z2` is the standardised amount random effect.
#' @param wkend_prop Weekend weight.
#' @return List of n_subjects x n_sims matrices `t` (usual intake), `p`
#'   (probability), `a` (amount).
#' @keywords internal
.simulate_usual <- function(pred, fit, z1, z2, wkend_prop = 3 / 7) {
  lambda <- fit$lambda
  if (fit$model_type == "amount") {
    sd1 <- 0; sd2 <- sqrt(fit$sigma2_b); ve <- fit$sigma2_w; rho <- 0
  } else {
    sd1 <- sqrt(fit$sigma2_v1); sd2 <- sqrt(fit$sigma2_v2)
    ve  <- fit$sigma2_e; rho <- fit$rho
  }
  # z2 drives the amount random effect; u1 is correlated with it via rho
  u2 <- sd2 * z2
  u1 <- sd1 * (rho * z2 + sqrt(1 - rho^2) * z1)

  day <- function(x1, x2) {
    p <- if (fit$model_type == "amount") 1 else stats::plogis(x1 + u1)
    a <- backtransform_9pt(x2 + u2, lambda, ve, fit$min_amt)
    if (length(p) == 1L) p <- matrix(1, nrow(a), ncol(a))
    list(p = p, a = a, t = p * a)
  }
  d0 <- day(pred$x1b1, pred$x2b2)
  if (!"x2b2_1" %in% names(pred)) return(d0)
  d1 <- day(pred$x1b1_1, pred$x2b2_1)
  wd <- 1 - wkend_prop
  list(t = wd * d0$t + wkend_prop * d1$t,
       p = wd * d0$p + wkend_prop * d1$p,
       a = wd * d0$a + wkend_prop * d1$a)
}


#' Compute distribution summary statistics
#'
#' Weighted mean, SD, percentiles and cutpoint proportions of simulated usual
#' intakes, using the DISTRIB PC-macro interpolation.
#'
#' @keywords internal
compute_distribution_summary <- function(sim_vals, weights, percentiles, cutpoints) {

  valid <- !is.na(sim_vals) & !is.na(weights)
  sim_vals <- sim_vals[valid]
  weights  <- weights[valid]

  na_cut <- function(prefix) if (is.null(cutpoints)) NULL else
    stats::setNames(rep(NA_real_, length(cutpoints)), paste0(prefix, cutpoints))
  if (length(sim_vals) == 0L || sum(weights) <= 0) {
    warning("No valid simulated values (or zero total weight); ",
            "returning NA for all summary statistics.", call. = FALSE)
    return(list(
      mean = NA_real_, sd = NA_real_,
      percentiles = stats::setNames(rep(NA_real_, length(percentiles)),
                                    paste0("p", percentiles)),
      cutpoint_below = na_cut("pct_below_"),
      cutpoint_above = na_cut("pct_above_"),
      n_simulated = 0L
    ))
  }

  wmean <- sum(sim_vals * weights) / sum(weights)
  wsd   <- sqrt(sum(weights * (sim_vals - wmean)^2) / sum(weights))

  ec <- .sas_ecdf(sim_vals, weights)
  wpercentiles <- .sas_percentiles(ec, percentiles / 100)
  names(wpercentiles) <- paste0("p", percentiles)

  cutpoint_below <- cutpoint_above <- NULL
  if (!is.null(cutpoints)) {
    cutpoint_below <- .sas_cutprob(ec, cutpoints)
    cutpoint_above <- 1 - cutpoint_below
    names(cutpoint_below) <- paste0("pct_below_", cutpoints)
    names(cutpoint_above) <- paste0("pct_above_", cutpoints)
  }

  list(
    mean = wmean,
    sd = wsd,
    percentiles = wpercentiles,
    cutpoint_below = cutpoint_below,
    cutpoint_above = cutpoint_above,
    n_simulated = length(sim_vals)
  )
}


# Unique sorted values with their mid-weight cumulative proportions
# (DISTRIB PC macro: wsumplush = wsum + adjwt_ties / 2).
.sas_ecdf <- function(x, w) {
  o  <- order(x)
  x  <- x[o]; w <- w[o] / sum(w)
  ux <- unique(x)
  wt <- rowsum(w, match(x, ux), reorder = FALSE)[, 1]
  list(x = ux, pos = cumsum(wt) - wt / 2)
}

# Linear interpolation of x at cumulative proportion h, extrapolating along
# the first/last segment, bounded below by 0 (tpercentile = max(0, ...)).
.sas_percentiles <- function(ec, probs) {
  x <- ec$x; pos <- ec$pos; n <- length(x)
  if (n == 1L) return(rep(max(0, x), length(probs)))
  vapply(probs, function(h) {
    k <- which(pos >= h)[1]
    if (is.na(k)) k <- n
    k <- max(k, 2L)
    slope <- (pos[k] - pos[k - 1]) / (x[k] - x[k - 1])
    max(0, x[k - 1] + (h - pos[k - 1]) / slope)
  }, numeric(1))
}

# Interpolated cumulative proportion at each cutpoint, clipped to [0, 1].
.sas_cutprob <- function(ec, cuts) {
  x <- ec$x; pos <- ec$pos; n <- length(x)
  if (n == 1L) return(as.numeric(cuts > x))
  vapply(cuts, function(cp) {
    k <- which(x >= cp)[1]
    if (is.na(k)) k <- n
    k <- max(k, 2L)
    slope <- (pos[k] - pos[k - 1]) / (x[k] - x[k - 1])
    min(1, max(0, pos[k - 1] + slope * (cp - x[k - 1])))
  }, numeric(1))
}


#' Weighted quantiles
#'
#' Weighted quantiles using the DISTRIB interpolation (mid-weight cumulative
#' proportions, linear interpolation). For equal weights and distinct values
#' this equals `quantile(x, probs, type = 5)` inside the data range.
#'
#' @keywords internal
weighted_quantiles <- function(x, w, probs) {
  if (!is.finite(sum(w)) || sum(w) <= 0) {
    warning("All survey weights are zero or non-finite; falling back to unweighted quantiles.",
            call. = FALSE)
    w <- rep(1, length(w))
  }
  .sas_percentiles(.sas_ecdf(x, w), probs)
}


#' Print method for distrib_result
#' @export
print.distrib_result <- function(x, ...) {
  cat("NCI Usual Intake Distribution (DISTRIB)\n")
  cat(sprintf("  Model type: %s\n", x$mixtran_obj$model_type))
  cat(sprintf("  Intake var: %s\n", x$mixtran_obj$intake_var))
  cat(sprintf("  N sims:     %d\n", x$n_sims))
  cat(sprintf("  Seed:       %d\n", x$seed))
  cat("\n")

  for (grp in names(x$results)) {
    r <- x$results[[grp]]
    cat(sprintf("--- %s ---\n", grp))
    cat(sprintf("  Mean: %.4f  (SD: %.4f)\n", r$mean, r$sd))
    cat("  Percentiles:\n")
    for (nm in names(r$percentiles)) {
      cat(sprintf("    %s: %.4f\n", nm, r$percentiles[nm]))
    }
    if (!is.null(r$cutpoint_below)) {
      cat("  Cutpoints (% below):\n")
      for (nm in names(r$cutpoint_below)) {
        cat(sprintf("    %s: %.2f%%\n", nm, r$cutpoint_below[nm] * 100))
      }
    }
    cat("\n")
  }
  invisible(x)
}


#' Extract results as a data frame
#'
#' @param distrib_obj Output from distrib()
#' @return Data frame with one row per subgroup
#' @export
distrib_to_df <- function(distrib_obj) {
  rows <- lapply(names(distrib_obj$results), function(grp) {
    r <- distrib_obj$results[[grp]]
    row <- data.frame(subgroup = grp, mean = r$mean, sd = r$sd,
                      stringsAsFactors = FALSE)
    for (nm in names(r$percentiles)) row[[nm]] <- r$percentiles[[nm]]
    if (!is.null(r$cutpoint_below)) {
      for (nm in names(r$cutpoint_below)) row[[nm]] <- r$cutpoint_below[[nm]]
    }
    row
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
