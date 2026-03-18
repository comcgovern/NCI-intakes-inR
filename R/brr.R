#' BRR Standard Errors for Usual Intake Distribution Statistics
#'
#' Implements Balanced Repeated Replication (BRR) with Fay's modification
#' to compute standard errors for percentiles, means, and cutpoint proportions
#' of the estimated usual intake distribution. This mirrors the approach used
#' in the NCI SAS macros for NHANES complex survey analyses.
#'
#' BRR works by running the full MIXTRAN -> DISTRIB pipeline K times, each
#' time using a different set of replicate weights. The variance of each
#' statistic is estimated from the spread of replicate estimates around the
#' full-sample estimate.
#'
#' @name brr
NULL


#' Compute BRR standard errors for usual intake distribution
#'
#' @param data Data frame in long format (same as passed to mixtran())
#' @param replicate_weights Matrix of BRR replicate weights. Rows correspond
#'   to subjects in the order they appear in unique(data[[subject_var]]).
#'   Columns are replicates (K columns). Alternatively, a data frame with
#'   subject IDs in the first column and replicate weights in remaining columns.
#' @param subject_var Name of the subject ID variable (character). Used to
#'   merge replicate weights to data rows.
#' @param mixtran_args Named list of arguments to pass to mixtran(). Must
#'   include at minimum intake_var, subject_var, repeat_var. Do not include
#'   weight_var here; that is handled internally.
#' @param distrib_args Named list of arguments to pass to distrib(). Optional.
#' @param fay_factor Fay's adjustment factor (default 0.3, the NHANES standard)
#' @param parallel Logical. If TRUE, use parallel::mclapply() for replicate
#'   runs. If the future.apply package is available, future_lapply() is used
#'   instead (supports more parallel backends). Default FALSE.
#' @param n_cores Number of cores to use when parallel = TRUE. Default uses
#'   half the available cores (at least 1).
#' @param verbose Print progress messages. Default TRUE.
#' @return A `brr_result` object with full-sample estimates, replicate
#'   estimates, BRR standard errors, and 95% confidence intervals.
#' @export
brr_usual_intake <- function(data,
                             replicate_weights,
                             subject_var,
                             mixtran_args,
                             distrib_args = list(),
                             fay_factor = 0.3,
                             parallel = FALSE,
                             n_cores = NULL,
                             verbose = TRUE) {

  data <- as.data.frame(data)

  # --- Parse replicate weights ---
  rw <- parse_replicate_weights(replicate_weights, data, subject_var)
  K <- ncol(rw$weights_matrix)
  subj_ids <- rw$subject_ids

  if (verbose) {
    message(sprintf("BRR: %d replicates, Fay factor = %.2f", K, fay_factor))
  }

  # --- Full-sample estimate ---
  if (verbose) message("Running full-sample estimate...")

  full_args <- c(list(data = data), mixtran_args)
  full_fit  <- do.call(mixtran, full_args)

  full_dist_args <- c(list(mixtran_obj = full_fit), distrib_args)
  full_dist <- do.call(distrib, full_dist_args)
  theta0    <- distrib_to_vector(full_dist)

  if (verbose) {
    message(sprintf("  Full-sample statistics computed (%d values)", length(theta0)))
  }

  # --- Replicate estimates ---
  if (verbose) message(sprintf("Running %d BRR replicates...", K))

  run_replicate <- function(k) {
    tryCatch({
      # Merge replicate weight k into data rows
      wt_k <- rw$weights_matrix[, k]
      names(wt_k) <- subj_ids
      data_k <- data
      data_k$.brr_wt_ <- wt_k[as.character(data_k[[subject_var]])]

      # Subjects with weight 0 should be excluded from model fitting
      # (Fay's method perturbs weights but doesn't zero them; however the
      #  user may supply exact BRR weights with some zeros from Hadamard
      #  construction — we keep zero-weight subjects in the model but give
      #  them down-weighted influence)
      args_k <- c(list(data = data_k), mixtran_args)
      args_k$weight_var <- ".brr_wt_"

      fit_k  <- do.call(mixtran, args_k)
      dist_args_k <- c(list(mixtran_obj = fit_k), distrib_args)
      dist_k <- do.call(distrib, dist_args_k)
      distrib_to_vector(dist_k)
    }, error = function(e) {
      warning(sprintf("Replicate %d failed: %s", k, e$message))
      rep(NA_real_, length(theta0))
    })
  }

  if (parallel) {
    n_cores <- n_cores %||% max(1L, parallel::detectCores() %/% 2L)
    if (requireNamespace("future.apply", quietly = TRUE)) {
      theta_k_list <- future.apply::future_lapply(
        seq_len(K), run_replicate,
        future.seed = TRUE
      )
    } else {
      theta_k_list <- parallel::mclapply(
        seq_len(K), run_replicate,
        mc.cores = n_cores
      )
    }
  } else {
    theta_k_list <- lapply(seq_len(K), function(k) {
      if (verbose && K > 5 && k %% max(1L, K %/% 5L) == 0) {
        message(sprintf("  Replicate %d / %d", k, K))
      }
      run_replicate(k)
    })
  }

  # Assemble replicate matrix (rows = statistics, cols = replicates)
  theta_matrix <- do.call(cbind, theta_k_list)

  # --- BRR variance with Fay's modification ---
  # Var_BRR(theta) = 1 / (K * (1-F)^2) * sum_k (theta_k - theta_0)^2
  fay_denom <- K * (1 - fay_factor)^2
  diff_sq <- (theta_matrix - theta0)^2
  # Handle NAs (failed replicates) — average over available replicates
  brr_var <- rowSums(diff_sq, na.rm = TRUE) /
    (fay_denom * rowMeans(!is.na(theta_matrix)))
  brr_se  <- sqrt(brr_var)

  # 95% CI (normal approximation)
  ci_lower <- theta0 - 1.96 * brr_se
  ci_upper <- theta0 + 1.96 * brr_se

  n_failed <- sum(apply(theta_matrix, 2, function(x) any(is.na(x))))
  if (n_failed > 0 && verbose) {
    message(sprintf("  Warning: %d replicate(s) failed and were excluded.", n_failed))
  }

  result <- list(
    estimates    = theta0,
    se           = brr_se,
    ci_lower     = ci_lower,
    ci_upper     = ci_upper,
    theta_matrix = theta_matrix,
    K            = K,
    K_failed     = n_failed,
    fay_factor   = fay_factor,
    full_fit     = full_fit,
    full_dist    = full_dist
  )
  class(result) <- "brr_result"
  return(result)
}


#' Parse replicate weights into a standardised matrix form
#' @keywords internal
parse_replicate_weights <- function(replicate_weights, data, subject_var) {

  if (is.data.frame(replicate_weights)) {
    # First column is subject ID, rest are replicate weights
    subject_ids   <- as.character(replicate_weights[[1]])
    weights_matrix <- as.matrix(replicate_weights[, -1, drop = FALSE])
  } else if (is.matrix(replicate_weights)) {
    # Must match unique subjects in data in order
    subject_ids <- as.character(unique(data[[subject_var]]))
    if (nrow(replicate_weights) != length(subject_ids)) {
      stop(
        "replicate_weights matrix has ", nrow(replicate_weights), " rows but ",
        "there are ", length(subject_ids), " unique subjects in data. ",
        "Pass a data frame with a subject ID column to avoid ambiguity."
      )
    }
    weights_matrix <- replicate_weights
  } else {
    stop("replicate_weights must be a matrix or data frame.")
  }

  list(subject_ids = subject_ids, weights_matrix = weights_matrix)
}


#' Flatten a distrib_result into a named numeric vector for BRR
#' @keywords internal
distrib_to_vector <- function(distrib_obj) {
  out <- c()
  for (grp in names(distrib_obj$results)) {
    r    <- distrib_obj$results[[grp]]
    pfx  <- if (grp == "Overall") "" else paste0(grp, ".")

    out[paste0(pfx, "mean")] <- r$mean
    out[paste0(pfx, "sd")]   <- r$sd

    for (nm in names(r$percentiles)) {
      out[paste0(pfx, nm)] <- r$percentiles[[nm]]
    }

    if (!is.null(r$cutpoint_below)) {
      for (nm in names(r$cutpoint_below)) {
        out[paste0(pfx, nm)] <- r$cutpoint_below[[nm]]
      }
    }
  }
  out
}


#' Null-coalescing operator (internal)
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b


#' Print method for brr_result
#' @export
print.brr_result <- function(x, ...) {
  cat("BRR Standard Errors — NCI Usual Intake\n")
  cat(sprintf("  Replicates:  %d  (failed: %d)\n", x$K, x$K_failed))
  cat(sprintf("  Fay factor:  %.2f\n", x$fay_factor))
  cat("\n")
  cat(sprintf("  %-28s %10s %8s  %s\n", "Statistic", "Estimate", "SE", "95% CI"))
  cat(sprintf("  %s\n", strrep("-", 65)))

  nms <- names(x$estimates)
  for (nm in nms) {
    cat(sprintf("  %-28s %10.3f %8.3f  [%.3f, %.3f]\n",
                nm,
                x$estimates[nm],
                x$se[nm],
                x$ci_lower[nm],
                x$ci_upper[nm]))
  }
  invisible(x)
}


#' Extract BRR results as a tidy data frame
#'
#' @param brr_obj Output from brr_usual_intake()
#' @return Data frame with one row per statistic
#' @export
brr_to_df <- function(brr_obj) {
  data.frame(
    statistic = names(brr_obj$estimates),
    estimate  = as.numeric(brr_obj$estimates),
    se        = as.numeric(brr_obj$se),
    ci_lower  = as.numeric(brr_obj$ci_lower),
    ci_upper  = as.numeric(brr_obj$ci_upper),
    stringsAsFactors = FALSE
  )
}


#' Generate BRR replicate weights from NHANES strata/PSU structure
#'
#' Creates BRR replicate weights using the Fay modification of the
#' balanced half-sample method. Requires NHANES masked variance unit
#' variables (SDMVSTRA, SDMVPSU).
#'
#' This is a thin wrapper around survey::svrepdesign() / survey::as.svrepdesign()
#' that returns the replicate weight matrix in the format expected by
#' brr_usual_intake().
#'
#' @param data Data frame containing NHANES design variables
#' @param strata_var Name of the stratum variable (default "SDMVSTRA")
#' @param psu_var Name of the PSU variable (default "SDMVPSU")
#' @param weight_var Name of the sampling weight variable
#' @param subject_var Name of the subject ID variable
#' @param fay_factor Fay modification factor (default 0.3)
#' @return Data frame with subject IDs in column 1 and replicate weights
#'   in subsequent columns, suitable for passing to brr_usual_intake().
#' @export
generate_nhanes_brr_weights <- function(data,
                                        strata_var  = "SDMVSTRA",
                                        psu_var     = "SDMVPSU",
                                        weight_var,
                                        subject_var,
                                        fay_factor  = 0.3) {

  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("The 'survey' package is required. Install with: install.packages('survey')")
  }

  # One row per subject (design variables are at person level)
  subj_data <- data[!duplicated(data[[subject_var]]), ]
  subj_data <- subj_data[order(subj_data[[subject_var]]), ]

  svy_design <- survey::svydesign(
    ids     = stats::as.formula(paste("~", psu_var)),
    strata  = stats::as.formula(paste("~", strata_var)),
    weights = stats::as.formula(paste("~", weight_var)),
    data    = subj_data,
    nest    = TRUE
  )

  brr_design <- survey::as.svrepdesign(
    svy_design,
    type       = "Fay",
    fay.rho    = fay_factor,
    compress   = FALSE
  )

  wt_matrix <- stats::weights(brr_design, type = "analysis")

  result <- as.data.frame(wt_matrix)
  colnames(result) <- paste0("brr_wt_", seq_len(ncol(result)))
  result <- cbind(
    data.frame(subject = subj_data[[subject_var]], stringsAsFactors = FALSE),
    result
  )
  names(result)[1] <- subject_var

  result
}
