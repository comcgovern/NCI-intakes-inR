#' DISTRIB: Estimate Usual Intake Distributions
#'
#' R port of the NCI DISTRIB SAS macro. Uses parameter estimates from
#' mixtran() and Monte Carlo simulation to estimate the distribution of
#' usual intake for a population or subpopulation.
#'
#' The key insight: observed 24-hour recall data conflates between-person
#' and within-person variation. DISTRIB removes within-person variation via
#' simulation, yielding the distribution of *usual* (long-run average) intake.
#'
#' @name distrib
NULL

#' Estimate usual intake distribution (DISTRIB equivalent)
#'
#' @param mixtran_obj Output from mixtran()
#' @param subgroup_var Optional: name of subgroup variable in the predicted data
#'   for which to compute separate distributions
#' @param cutpoints Numeric vector of cutpoints for computing % above/below
#' @param percentiles Percentiles to estimate (default: standard set)
#' @param n_sims Number of Monte Carlo simulations per person (default: 100)
#' @param seed Random seed for reproducibility
#' @return A `distrib_result` object with estimated usual intake distribution
#' @export
distrib <- function(mixtran_obj,
                    subgroup_var = NULL,
                    cutpoints = NULL,
                    percentiles = c(5, 10, 25, 50, 75, 90, 95),
                    n_sims = 100,
                    seed = 12345) {

  stopifnot(inherits(mixtran_obj, "mixtran_fit"))
  set.seed(seed)

  pred <- mixtran_obj$predicted
  model_type <- mixtran_obj$model_type
  lambda <- mixtran_obj$lambda

  # --- Remove nuisance covariate effects (NCI convention) ---
  # DISTRIB estimates usual intake at reference covariate values:
  # seq_num = 0 (first recall) and weekend = 0 (weekday).
  # User covariates (demographics) are retained at observed values.
  pred <- .remove_nuisance_effects_distrib(pred, mixtran_obj)

  # Get unique subjects with their linear predictors and weights
  subj_df <- get_subject_level_data(pred, model_type, mixtran_obj)

  # Determine subgroups
  if (!is.null(subgroup_var)) {
    if (!(subgroup_var %in% names(pred))) {
      stop("Subgroup variable '", subgroup_var, "' not found in predicted data.")
    }
    # Attach subgroup to subject-level data
    subj_map <- pred[!duplicated(pred$subject), c("subject", subgroup_var)]
    subj_df <- merge(subj_df, subj_map, by = "subject")
    subgroup_levels <- c(unique(subj_df[[subgroup_var]]), "Overall")
  } else {
    subgroup_levels <- "Overall"
  }

  # --- Monte Carlo simulation ---
  if (model_type == "amount") {
    simulated <- simulate_amount_model(
      subj_df = subj_df,
      beta = mixtran_obj$beta,
      sigma2_b = mixtran_obj$sigma2_b,
      sigma2_w = mixtran_obj$sigma2_w,
      lambda = lambda,
      n_sims = n_sims
    )
  } else {
    simulated <- simulate_twopart_model(
      subj_df = subj_df,
      alpha = mixtran_obj$alpha,
      beta = mixtran_obj$beta,
      sigma2_v1 = mixtran_obj$sigma2_v1,
      sigma2_v2 = mixtran_obj$sigma2_v2,
      sigma2_e = mixtran_obj$sigma2_e,
      rho = mixtran_obj$rho,
      lambda = lambda,
      n_sims = n_sims
    )
  }

  # Compute summaries by subgroup
  results <- list()
  for (grp in subgroup_levels) {
    if (grp == "Overall") {
      sim_vals <- simulated$usual_intake
      wts <- rep(subj_df$weight, each = n_sims)
    } else {
      grp_idx <- which(subj_df[[subgroup_var]] == grp)
      sim_idx <- as.vector(outer(seq_len(n_sims), (grp_idx - 1) * n_sims, "+"))
      sim_vals <- simulated$usual_intake[sim_idx]
      wts <- rep(subj_df$weight[grp_idx], each = n_sims)
    }

    results[[grp]] <- compute_distribution_summary(
      sim_vals, wts, percentiles, cutpoints
    )
  }

  output <- list(
    results = results,
    simulated = simulated,
    mixtran_obj = mixtran_obj,
    n_sims = n_sims,
    seed = seed,
    percentiles_requested = percentiles,
    cutpoints = cutpoints,
    subgroup_var = subgroup_var
  )
  class(output) <- "distrib_result"
  return(output)
}


#' Remove nuisance covariate effects from predicted linear predictors
#'
#' NCI convention: DISTRIB computes usual intake at reference covariate values
#' (seq_num = 0, weekend = 0). This function subtracts these effects from the
#' population-level linear predictors so that the Monte Carlo simulation
#' produces usual intake for a reference ("usual") day.
#'
#' @keywords internal
.remove_nuisance_effects_distrib <- function(pred, mixtran_obj) {
  model_type <- mixtran_obj$model_type

  if (model_type == "amount") {
    beta <- mixtran_obj$beta
    # Remove seq_num effect
    if ("seq_num" %in% names(beta) && "seq_num" %in% names(pred)) {
      pred$linpred <- pred$linpred - beta["seq_num"] * pred$seq_num
    }
    # Remove weekend effect
    if ("weekend" %in% names(beta) && "weekend" %in% names(pred)) {
      pred$linpred <- pred$linpred - beta["weekend"] * pred$weekend
    }
  } else {
    alpha <- mixtran_obj$alpha
    beta  <- mixtran_obj$beta

    # Probability sub-model: remove seq_num and weekend from prob_linpred
    if ("seq_num" %in% names(alpha) && "seq_num" %in% names(pred)) {
      pred$prob_linpred <- pred$prob_linpred - alpha["seq_num"] * pred$seq_num
    }
    if ("weekend" %in% names(alpha) && "weekend" %in% names(pred)) {
      pred$prob_linpred <- pred$prob_linpred - alpha["weekend"] * pred$weekend
    }

    # Amount sub-model: remove seq_num and weekend from amt_linpred
    # (only for consumed rows where amt_linpred is not NA)
    pos_mask <- !is.na(pred$amt_linpred)
    if ("seq_num" %in% names(beta) && "seq_num" %in% names(pred)) {
      pred$amt_linpred[pos_mask] <- pred$amt_linpred[pos_mask] -
        beta["seq_num"] * pred$seq_num[pos_mask]
    }
    if ("weekend" %in% names(beta) && "weekend" %in% names(pred)) {
      pred$amt_linpred[pos_mask] <- pred$amt_linpred[pos_mask] -
        beta["weekend"] * pred$weekend[pos_mask]
    }
  }

  pred
}


#' Get subject-level data from predicted dataset
#' @keywords internal
get_subject_level_data <- function(pred, model_type, mixtran_obj = NULL) {
  # For each subject, get one row with their linear predictor(s) and weight.
  # Nuisance effects (seq_num, weekend) should already be removed by
  # .remove_nuisance_effects_distrib() before calling this function.

  subjects <- unique(pred$subject)

  if (model_type == "amount") {
    subj_df <- do.call(rbind, lapply(subjects, function(s) {
      rows <- pred[pred$subject == s, ]
      data.frame(
        subject = s,
        linpred = mean(rows$linpred),
        weight = rows$weight[1],
        stringsAsFactors = FALSE
      )
    }))
  } else {
    # For subjects who never consumed, amt_linpred is all NA. Use the
    # population-level intercept (beta_0) as a fallback — this is the NCI
    # convention: never-consumers' amount component is the population mean
    # on the transformed scale.
    global_amt_mean <- mean(pred$amt_linpred, na.rm = TRUE)
    if (is.nan(global_amt_mean)) {
      warning("All amt_linpred values are NA (no consumers in sample); ",
              "amount sub-model estimates will be the population intercept only.",
              call. = FALSE)
    }

    subj_df <- do.call(rbind, lapply(subjects, function(s) {
      rows <- pred[pred$subject == s, ]
      amt_lp <- mean(rows$amt_linpred, na.rm = TRUE)
      # NaN means all amt_linpred were NA (never-consumer): use population mean
      if (is.nan(amt_lp)) amt_lp <- global_amt_mean
      data.frame(
        subject = s,
        prob_linpred = mean(rows$prob_linpred, na.rm = TRUE),
        amt_linpred = amt_lp,
        weight = rows$weight[1],
        stringsAsFactors = FALSE
      )
    }))
  }
  return(subj_df)
}


#' Simulate usual intakes from amount-only model
#'
#' For each person i, simulate n_sims draws of their usual intake:
#'   T(UI_i^(r)) = X_i'beta + u_i^(r)
#'   where u_i^(r) ~ N(0, sigma2_b)
#'   then back-transform: UI_i^(r) = T^{-1}(T(UI_i^(r)))
#'
#' Note: Within-person error (epsilon) is NOT added because we want
#' usual (habitual) intake, not single-day intake.
#'
#' Vectorized: all n_subj × n_sims random effects drawn in one call;
#' the Box-Cox inverse is applied to the full matrix — no per-subject loop.
#'
#' @keywords internal
simulate_amount_model <- function(subj_df, beta, sigma2_b, sigma2_w,
                                   lambda, n_sims) {

  n_subj <- nrow(subj_df)
  sd_b   <- sqrt(sigma2_b)

  # Draw all between-person random effects at once: n_subj × n_sims
  u_mat <- matrix(stats::rnorm(n_subj * n_sims, mean = 0, sd = sd_b),
                  nrow = n_subj, ncol = n_sims)

  # Broadcast: linpred[i] is added to row i of u_mat (R recycling is column-major,
  # so we rely on the fact that linpred[i] is added to u_mat[i, j] for all j)
  t_usual_mat <- subj_df$linpred + u_mat  # n_subj × n_sims

  # Back-transform and clamp
  ui_mat <- pmax(boxcox_inverse(t_usual_mat, lambda), 0)

  # Flatten row-major: subject i occupies positions (i-1)*n_sims+1 : i*n_sims
  usual_intake <- as.vector(t(ui_mat))

  list(
    usual_intake = usual_intake,
    n_subjects   = n_subj,
    n_sims       = n_sims
  )
}


#' Simulate usual intakes from two-part model
#'
#' For each person, simulate correlated random effects and compute:
#'   P(consume) = logistic(Z'alpha + v1_i)
#'   Amount|consume = T^{-1}(X'beta + v2_i)
#'   Usual intake = P(consume) * E[Amount|consume]
#'
#' Vectorized: all n_subj × n_sims random effects drawn as matrices;
#' probability and amount computed via matrix operations — no per-subject loop.
#'
#' @keywords internal
simulate_twopart_model <- function(subj_df, alpha, beta,
                                    sigma2_v1, sigma2_v2, sigma2_e,
                                    rho, lambda, n_sims) {

  n_subj    <- nrow(subj_df)
  sd_v1     <- sqrt(sigma2_v1)
  sd_v2     <- sqrt(sigma2_v2)
  rho_safe  <- min(max(rho, -0.9999), 0.9999)
  sqrt_1r2  <- sqrt(1 - rho_safe^2)

  # Draw independent standard normals: n_subj × n_sims
  z1_mat <- matrix(stats::rnorm(n_subj * n_sims), nrow = n_subj, ncol = n_sims)
  z2_mat <- matrix(stats::rnorm(n_subj * n_sims), nrow = n_subj, ncol = n_sims)

  # Correlated random effects via Cholesky (bivariate normal)
  v1_mat <- sd_v1 * z1_mat
  v2_mat <- sd_v2 * (rho_safe * z1_mat + sqrt_1r2 * z2_mat)

  # Probability: broadcast prob_linpred[i] across n_sims columns
  eta_prob_mat <- subj_df$prob_linpred + v1_mat  # n_subj × n_sims
  prob_mat     <- 1 / (1 + exp(-eta_prob_mat))

  # Amount: broadcast amt_linpred[i] across n_sims columns
  t_amount_mat <- subj_df$amt_linpred + v2_mat
  amount_mat   <- pmax(boxcox_inverse(t_amount_mat, lambda), 0)

  # Usual intake = P(consume) * E[amount | consume]
  ui_mat <- prob_mat * amount_mat

  # Flatten row-major: subject i occupies positions (i-1)*n_sims+1 : i*n_sims
  usual_intake <- as.vector(t(ui_mat))

  list(
    usual_intake = usual_intake,
    n_subjects   = n_subj,
    n_sims       = n_sims
  )
}


#' Compute distribution summary statistics
#'
#' Weighted percentiles, mean, and cutpoint proportions from simulated data
#'
#' @keywords internal
compute_distribution_summary <- function(sim_vals, weights, percentiles, cutpoints) {

  # Remove NAs
  valid <- !is.na(sim_vals)
  sim_vals <- sim_vals[valid]
  weights <- weights[valid]

  if (length(sim_vals) == 0L) {
    warning("No valid simulated values after removing NAs; ",
            "returning NA for all summary statistics.", call. = FALSE)
    na_pct <- setNames(rep(NA_real_, length(percentiles)), paste0("p", percentiles))
    return(list(
      mean = NA_real_, sd = NA_real_,
      percentiles = na_pct,
      cutpoint_below = if (!is.null(cutpoints)) setNames(rep(NA_real_, length(cutpoints)), paste0("pct_below_", cutpoints)) else NULL,
      cutpoint_above = if (!is.null(cutpoints)) setNames(rep(NA_real_, length(cutpoints)), paste0("pct_above_", cutpoints)) else NULL,
      n_simulated = 0L
    ))
  }

  # Weighted mean
  wmean <- stats::weighted.mean(sim_vals, weights)

  # Weighted percentiles
  wpercentiles <- weighted_quantiles(sim_vals, weights, percentiles / 100)
  names(wpercentiles) <- paste0("p", percentiles)

  # Weighted SD (for reference)
  wvar <- stats::weighted.mean((sim_vals - wmean)^2, weights)
  wsd <- sqrt(wvar)

  # Cutpoint proportions
  if (!is.null(cutpoints)) {
    cutpoint_below <- vapply(cutpoints, function(cp) {
      stats::weighted.mean(sim_vals < cp, weights)
    }, numeric(1))
    cutpoint_above <- 1 - cutpoint_below
    names(cutpoint_below) <- paste0("pct_below_", cutpoints)
    names(cutpoint_above) <- paste0("pct_above_", cutpoints)
  } else {
    cutpoint_below <- NULL
    cutpoint_above <- NULL
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


#' Weighted quantiles
#'
#' Compute weighted quantiles using linear interpolation
#' @keywords internal
weighted_quantiles <- function(x, w, probs) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]

  # Cumulative normalized weights
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    warning("All survey weights are zero or non-finite; falling back to unweighted quantiles.",
            call. = FALSE)
    w <- rep(1, length(w))
    sw <- length(w)
  }
  cw <- cumsum(w) / sw

  vapply(probs, function(p) {
    if (p <= cw[1]) return(x[1])
    if (p >= cw[length(cw)]) return(x[length(x)])
    idx <- max(which(cw <= p))
    if (idx >= length(x)) return(x[length(x)])
    # Linear interpolation
    frac <- (p - cw[idx]) / (cw[idx + 1] - cw[idx])
    x[idx] + frac * (x[idx + 1] - x[idx])
  }, numeric(1))
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
    cat(sprintf("  Mean: %.2f  (SD: %.2f)\n", r$mean, r$sd))
    cat("  Percentiles:\n")
    for (nm in names(r$percentiles)) {
      cat(sprintf("    %s: %.2f\n", nm, r$percentiles[nm]))
    }
    if (!is.null(r$cutpoint_below)) {
      cat("  Cutpoints (% below):\n")
      for (nm in names(r$cutpoint_below)) {
        cat(sprintf("    %s: %.1f%%\n", nm, r$cutpoint_below[nm] * 100))
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
    row <- data.frame(
      subgroup = grp,
      mean = r$mean,
      sd = r$sd,
      stringsAsFactors = FALSE
    )
    # Add percentiles
    for (nm in names(r$percentiles)) {
      row[[nm]] <- r$percentiles[nm]
    }
    # Add cutpoints
    if (!is.null(r$cutpoint_below)) {
      for (nm in names(r$cutpoint_below)) {
        row[[nm]] <- r$cutpoint_below[nm]
      }
    }
    row
  })
  do.call(rbind, rows)
}
