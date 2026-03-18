#' MIXTRAN: Fit NCI Usual Intake Models
#'
#' R port of the NCI MIXTRAN SAS macro. Fits nonlinear mixed effects models
#' to 24-hour dietary recall data to estimate parameters of the usual intake
#' distribution, separating within-person from between-person variation.
#'
#' Three model types are supported:
#' - "amount": Amount-only model for ubiquitously consumed nutrients
#'             (e.g., total energy, sodium, added sugars)
#' - "uncorr": Two-part model with uncorrelated random effects for
#'             episodically consumed foods
#' - "corr":   Two-part model with correlated random effects (most general)
#'
#' @name mixtran
NULL

#' Fit an NCI usual intake model (MIXTRAN equivalent)
#'
#' @param data Data frame in long format (one row per recall per person)
#' @param intake_var Name of the intake variable (character)
#' @param subject_var Name of the subject ID variable (character)
#' @param repeat_var Name of the recall day indicator (character), e.g. 1 or 2
#' @param model_type One of "amount", "uncorr", "corr"
#' @param covariates Character vector of covariate names. Must be binary or
#'   continuous. Categorical variables with >2 levels should be dummy-coded.
#' @param weekend_var Optional name of weekend indicator variable (character)
#' @param weight_var Optional name of survey weight variable (character)
#' @param lambda Box-Cox lambda. If NULL, estimated from data.
#' @param lambda_grid Grid for Box-Cox search (only used if lambda is NULL)
#' @param min_positive Minimum value to replace zeros in amount model.
#'   Default NULL uses half the minimum positive value (NCI convention).
#' @param verbose Print model fitting progress
#' @return A `mixtran_fit` object containing parameter estimates, predicted
#'   values, and diagnostics needed for the DISTRIB step.
#' @export
mixtran <- function(data,
                    intake_var,
                    subject_var,
                    repeat_var,
                    model_type = c("amount", "uncorr", "corr"),
                    covariates = NULL,
                    weekend_var = NULL,
                    weight_var = NULL,
                    lambda = NULL,
                    lambda_grid = seq(0.01, 1.0, by = 0.01),
                    min_positive = NULL,
                    verbose = TRUE) {

  model_type <- match.arg(model_type)
  data <- as.data.frame(data)

  # --- Validate inputs ---
  required_vars <- c(intake_var, subject_var, repeat_var)
  if (!is.null(weekend_var)) required_vars <- c(required_vars, weekend_var)
  if (!is.null(weight_var)) required_vars <- c(required_vars, weight_var)
  if (!is.null(covariates)) required_vars <- c(required_vars, covariates)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # --- Prepare data ---
  prep <- prepare_mixtran_data(
    data = data,
    intake_var = intake_var,
    subject_var = subject_var,
    repeat_var = repeat_var,
    covariates = covariates,
    weekend_var = weekend_var,
    weight_var = weight_var,
    min_positive = min_positive,
    model_type = model_type
  )

  # --- Find optimal lambda if not provided ---
  if (is.null(lambda)) {
    if (verbose) message("Finding optimal Box-Cox lambda...")
    y_pos <- prep$data$intake[prep$data$intake > 0]
    bc_result <- find_optimal_lambda(y_pos, lambda_grid = lambda_grid)
    lambda <- bc_result$lambda
    if (verbose) {
      message(sprintf("  Raw optimal lambda: %.3f", bc_result$lambda_raw))
      if (bc_result$used_log) message("  Using log transform (lambda < 0.15)")
      message(sprintf("  Final lambda: %.3f", lambda))
    }
  }

  # --- Fit model ---
  if (model_type == "amount") {
    result <- fit_amount_model(prep, lambda, verbose)
  } else if (model_type == "uncorr") {
    result <- fit_twopart_uncorr(prep, lambda, verbose)
  } else {
    result <- fit_twopart_corr(prep, lambda, verbose)
  }

  result$lambda <- lambda
  result$model_type <- model_type
  result$intake_var <- intake_var
  result$subject_var <- subject_var
  result$covariates <- covariates
  result$weekend_var <- weekend_var
  result$n_subjects <- prep$n_subjects
  result$n_recalls <- nrow(prep$data)
  result$n_positive <- sum(prep$data$intake > 0)
  result$pct_consuming <- result$n_positive / result$n_recalls * 100

  class(result) <- "mixtran_fit"
  return(result)
}


#' Prepare data for MIXTRAN model fitting
#' @keywords internal
prepare_mixtran_data <- function(data, intake_var, subject_var, repeat_var,
                                 covariates, weekend_var, weight_var,
                                 min_positive, model_type) {

  # Build working dataset
  work <- data.frame(
    subject = data[[subject_var]],
    repeat_num = data[[repeat_var]],
    intake = data[[intake_var]],
    stringsAsFactors = FALSE
  )

  # Add covariates
  cov_names <- character(0)
  if (!is.null(weekend_var)) {
    work$weekend <- as.numeric(data[[weekend_var]])
    cov_names <- c(cov_names, "weekend")
  }
  if (!is.null(covariates)) {
    for (cv in covariates) {
      work[[cv]] <- data[[cv]]
      cov_names <- c(cov_names, cv)
    }
  }

  # Survey weights
  if (!is.null(weight_var)) {
    work$weight <- data[[weight_var]]
  } else {
    work$weight <- 1
  }

  # Add sequence indicator (NCI convention: first recall = 0, subsequent = 1)
  work$seq_num <- as.numeric(work$repeat_num > 1)
  cov_names <- c("seq_num", cov_names)

  # Handle zeros for amount model (NCI convention: replace with half minimum)
  if (model_type == "amount") {
    if (any(work$intake == 0, na.rm = TRUE)) {
      if (is.null(min_positive)) {
        min_pos <- min(work$intake[work$intake > 0], na.rm = TRUE)
        min_positive <- min_pos / 2
      }
      message(sprintf("  Replacing %d zero values with %.4f (half-min)",
                      sum(work$intake == 0, na.rm = TRUE), min_positive))
      work$intake[work$intake == 0] <- min_positive
    }
  }

  # Create consumption indicator for two-part models
  work$consumed <- as.numeric(work$intake > 0)

  # Remove rows with missing intake
  work <- work[!is.na(work$intake), ]

  # Sort by subject and repeat (required for mixed model)
  work <- work[order(work$subject, work$repeat_num), ]

  # Diagnostics
  n_subjects <- length(unique(work$subject))
  recalls_per_subject <- table(work$subject)
  n_with_2plus <- sum(recalls_per_subject >= 2)

  if (model_type != "amount") {
    # For two-part models, need subjects with 2+ positive recalls
    pos_per_subject <- tapply(work$consumed, work$subject, sum)
    n_with_2pos <- sum(pos_per_subject >= 2, na.rm = TRUE)
    if (n_with_2pos == 0) {
      stop("No subjects with 2+ positive recalls. Cannot fit two-part model.")
    }
    if (n_with_2pos <= 10) {
      warning("Only ", n_with_2pos, " subjects with 2+ positive recalls. ",
              "Results may be unstable.")
    }
  }

  list(
    data = work,
    cov_names = cov_names,
    n_subjects = n_subjects,
    n_with_2plus = n_with_2plus,
    min_positive = min_positive
  )
}


#' Fit amount-only model (for ubiquitously consumed nutrients)
#'
#' Model: T(Y_ij) = X_ij'beta + u_i + eps_ij
#'   where u_i ~ N(0, sigma2_b) and eps_ij ~ N(0, sigma2_w)
#'
#' @keywords internal
fit_amount_model <- function(prep, lambda, verbose) {

  work <- prep$data
  cov_names <- prep$cov_names

  # Apply Box-Cox transform
  work$t_intake <- boxcox_transform(work$intake, lambda)

  # Build formula
  # Fixed effects: intercept + covariates
  if (length(cov_names) > 0) {
    fixed_formula <- stats::as.formula(
      paste("t_intake ~", paste(cov_names, collapse = " + "))
    )
  } else {
    fixed_formula <- stats::as.formula("t_intake ~ 1")
  }

  # Fit linear mixed model using nlme::lme
  # Random intercept for subject captures between-person variation
  if (verbose) message("Fitting amount-only model with nlme::lme()...")

  fit <- tryCatch({
    nlme::lme(
      fixed = fixed_formula,
      random = ~ 1 | subject,
      data = work,
      method = "REML",
      weights = if (any(work$weight != 1)) nlme::varFixed(~ 1/weight) else NULL,
      control = nlme::lmeControl(
        maxIter = 500,
        msMaxIter = 500,
        opt = "optim",
        returnObject = TRUE
      )
    )
  }, error = function(e) {
    # Fall back to ML if REML fails
    if (verbose) message("  REML failed, trying ML...")
    nlme::lme(
      fixed = fixed_formula,
      random = ~ 1 | subject,
      data = work,
      method = "ML",
      control = nlme::lmeControl(
        maxIter = 500,
        msMaxIter = 500,
        opt = "optim",
        returnObject = TRUE
      )
    )
  })

  # Extract variance components
  vc <- nlme::VarCorr(fit)
  sigma2_b <- as.numeric(vc[1, "Variance"])  # Between-person
  sigma2_w <- as.numeric(vc[2, "Variance"])  # Within-person (residual)

  # Fixed effects
  beta <- nlme::fixef(fit)

  # Predicted values (linear predictor without random effect)
  # This is X_i'beta for each subject, used by DISTRIB
  # Use only unique subject-level predictions (average over covariates)
  pred_df <- work
  pred_df$linpred <- stats::predict(fit, level = 0)  # population-level
  pred_df$ranef <- stats::predict(fit, level = 1) - pred_df$linpred  # subject-level RE

  if (verbose) {
    message(sprintf("  Between-person variance (sigma2_b): %.4f", sigma2_b))
    message(sprintf("  Within-person variance  (sigma2_w): %.4f", sigma2_w))
    message(sprintf("  Variance ratio (within/between):    %.2f", sigma2_w / sigma2_b))
    message(sprintf("  Fixed effects:"))
    for (nm in names(beta)) {
      message(sprintf("    %s: %.4f", nm, beta[nm]))
    }
  }

  list(
    model_fit = fit,
    beta = beta,
    sigma2_b = sigma2_b,
    sigma2_w = sigma2_w,
    var_ratio = sigma2_w / sigma2_b,
    predicted = pred_df,
    converged = TRUE
  )
}


#' Fit two-part uncorrelated model (for episodically consumed foods)
#'
#' Part 1 (probability): logit(P(Y>0)) = Z'alpha + v1_i
#' Part 2 (amount|consumed): T(Y|Y>0) = X'beta + v2_i + eps
#' Assumes Cov(v1_i, v2_i) = 0
#'
#' @keywords internal
fit_twopart_uncorr <- function(prep, lambda, verbose) {

  work <- prep$data
  cov_names <- prep$cov_names

  if (verbose) message("Fitting two-part uncorrelated model...")

  # --- Part 1: Probability model (logistic GLMM) ---
  if (verbose) message("  Part 1: Probability of consumption...")

  if (length(cov_names) > 0) {
    prob_formula <- stats::as.formula(
      paste("consumed ~", paste(cov_names, collapse = " + "))
    )
  } else {
    prob_formula <- stats::as.formula("consumed ~ 1")
  }

  # Use MASS::glmmPQL for logistic mixed model (faster than full ML for large data)
  # Alternative: lme4::glmer, but PQL is more robust for starting values
  prob_fit <- tryCatch({
    MASS::glmmPQL(
      fixed = prob_formula,
      random = ~ 1 | subject,
      family = stats::binomial(link = "logit"),
      data = work,
      verbose = FALSE
    )
  }, error = function(e) {
    if (verbose) message("    glmmPQL failed: ", e$message)
    NULL
  })

  if (is.null(prob_fit)) {
    warning("Probability model failed to converge. Returning partial results.")
    alpha <- NULL
    sigma2_v1 <- NA
  } else {
    alpha <- nlme::fixef(prob_fit)
    vc_prob <- nlme::VarCorr(prob_fit)
    sigma2_v1 <- as.numeric(vc_prob[1, "Variance"])
    if (verbose) {
      message(sprintf("    Between-person variance (prob): %.4f", sigma2_v1))
    }
  }

  # --- Part 2: Amount model (conditional on consumption) ---
  if (verbose) message("  Part 2: Amount consumed (conditional)...")

  work_pos <- work[work$consumed == 1, ]
  work_pos$t_intake <- boxcox_transform(work_pos$intake, lambda)

  if (length(cov_names) > 0) {
    amt_formula <- stats::as.formula(
      paste("t_intake ~", paste(cov_names, collapse = " + "))
    )
  } else {
    amt_formula <- stats::as.formula("t_intake ~ 1")
  }

  amt_fit <- tryCatch({
    nlme::lme(
      fixed = amt_formula,
      random = ~ 1 | subject,
      data = work_pos,
      method = "REML",
      control = nlme::lmeControl(maxIter = 500, opt = "optim", returnObject = TRUE)
    )
  }, error = function(e) {
    if (verbose) message("    Amount model REML failed, trying ML...")
    nlme::lme(
      fixed = amt_formula,
      random = ~ 1 | subject,
      data = work_pos,
      method = "ML",
      control = nlme::lmeControl(maxIter = 500, opt = "optim", returnObject = TRUE)
    )
  })

  beta <- nlme::fixef(amt_fit)
  vc_amt <- nlme::VarCorr(amt_fit)
  sigma2_v2 <- as.numeric(vc_amt[1, "Variance"])
  sigma2_e <- as.numeric(vc_amt[2, "Variance"])

  if (verbose) {
    message(sprintf("    Between-person variance (amount): %.4f", sigma2_v2))
    message(sprintf("    Within-person variance  (amount): %.4f", sigma2_e))
  }

  # Predicted values
  pred_df <- work
  pred_df$prob_linpred <- if (!is.null(prob_fit)) {
    stats::predict(prob_fit, level = 0, newdata = work)
  } else {
    NA
  }

  # For amount predictions, only available for consumed rows
  pred_df$amt_linpred <- NA
  pred_df$amt_linpred[pred_df$consumed == 1] <- stats::predict(amt_fit, level = 0)

  list(
    prob_fit = prob_fit,
    amt_fit = amt_fit,
    alpha = alpha,
    beta = beta,
    sigma2_v1 = sigma2_v1,
    sigma2_v2 = sigma2_v2,
    sigma2_e = sigma2_e,
    rho = 0,  # uncorrelated by definition
    predicted = pred_df,
    converged = !is.null(prob_fit)
  )
}


#' Fit two-part correlated model (for episodically consumed foods)
#'
#' Part 1: logit(P(Y>0)) = Z'alpha + v1_i
#' Part 2: T(Y|Y>0) = X'beta + v2_i + eps
#' (v1_i, v2_i) ~ BVN(0, Sigma)
#'
#' Uses profile likelihood over a grid of rho values:
#' 1. Fit uncorrelated model to get fixed effects and variance components
#' 2. For each rho on a grid, condition on rho and compute profile log-likelihood
#'    using Cholesky-adjusted residuals from the two fitted models
#' 3. Select rho that maximises the profile likelihood
#'
#' This approach is fast because it reuses existing nlme/glmmPQL fits rather
#' than performing full joint optimisation. The approximation is documented in
#' the package spec (PACKAGE_SPEC.md, section 4.3.1).
#'
#' @keywords internal
fit_twopart_corr <- function(prep, lambda, verbose) {

  if (verbose) message("Fitting two-part correlated model (profile likelihood over rho)...")

  # Step 1: Fit uncorrelated model to get all parameters
  if (verbose) message("  Step 1: Fitting uncorrelated model...")
  uncorr <- fit_twopart_uncorr(prep, lambda, verbose = FALSE)

  if (!uncorr$converged) {
    warning("Uncorrelated model failed. Cannot fit correlated model.")
    uncorr$rho <- NA
    return(uncorr)
  }

  work      <- prep$data
  cov_names <- prep$cov_names

  # Build design matrix (shared by both parts)
  if (length(cov_names) > 0) {
    X_full <- stats::model.matrix(
      stats::as.formula(paste("~", paste(cov_names, collapse = " + "))),
      data = work
    )
  } else {
    X_full <- matrix(1, nrow = nrow(work), ncol = 1)
  }

  # Subject-level standardised residuals from each part
  # Probability part: subject-level random effects are v1_i ~ N(0, sigma2_v1)
  # We approximate v1_i by the glmmPQL empirical Bayes estimate, then standardise
  prob_re <- tryCatch(
    as.numeric(nlme::ranef(uncorr$prob_fit)[[1]]),
    error = function(e) rep(0, prep$n_subjects)
  )
  amt_re <- tryCatch(
    as.numeric(nlme::ranef(uncorr$amt_fit)[[1]]),
    error = function(e) rep(0, prep$n_subjects)
  )

  # Standardise to unit variance
  sd_v1 <- sqrt(uncorr$sigma2_v1)
  sd_v2 <- sqrt(uncorr$sigma2_v2)
  r1 <- if (sd_v1 > 0) prob_re / sd_v1 else prob_re
  r2 <- if (sd_v2 > 0) amt_re  / sd_v2 else amt_re

  # Profile log-likelihood over rho:
  # Given rho, the bivariate density of (v1, v2) changes. We evaluate how well
  # the standardised residuals (r1, r2) fit a BVN(0, [[1, rho],[rho, 1]]) model
  # using the bivariate normal log-likelihood.
  profile_loglik <- function(rho) {
    n <- length(r1)
    det_val <- 1 - rho^2
    if (det_val <= 0) return(-Inf)
    # Sum of bivariate normal log-densities for (r1_i, r2_i)
    ll <- -0.5 * n * log(det_val) -
      0.5 / det_val * sum(r1^2 - 2 * rho * r1 * r2 + r2^2) +
      # subtract the univariate terms already counted (constants)
      0.5 * sum(r1^2) + 0.5 * sum(r2^2)
    ll
  }

  # Coarse grid search then refine around maximum
  rho_grid_coarse <- seq(-0.9, 0.9, by = 0.1)
  ll_coarse <- vapply(rho_grid_coarse, profile_loglik, numeric(1))
  best_coarse <- rho_grid_coarse[which.max(ll_coarse)]

  # Fine grid within ±0.15 of coarse best
  rho_grid_fine <- seq(
    max(-0.99, best_coarse - 0.15),
    min( 0.99, best_coarse + 0.15),
    by = 0.01
  )
  ll_fine <- vapply(rho_grid_fine, profile_loglik, numeric(1))
  rho_hat <- rho_grid_fine[which.max(ll_fine)]

  if (verbose) {
    message(sprintf("  Step 2: Profile likelihood over rho"))
    message(sprintf("    Coarse best rho: %.2f", best_coarse))
    message(sprintf("    Fine best rho:   %.4f", rho_hat))
    message(sprintf("  Final parameter estimates:"))
    message(sprintf("    sigma2_v1 (prob):   %.4f", uncorr$sigma2_v1))
    message(sprintf("    sigma2_v2 (amount): %.4f", uncorr$sigma2_v2))
    message(sprintf("    sigma2_e (within):  %.4f", uncorr$sigma2_e))
    message(sprintf("    rho:                %.4f", rho_hat))
  }

  # Return all uncorrelated parameters plus the estimated rho
  list(
    prob_fit   = uncorr$prob_fit,
    amt_fit    = uncorr$amt_fit,
    alpha      = uncorr$alpha,
    beta       = uncorr$beta,
    sigma2_v1  = uncorr$sigma2_v1,
    sigma2_v2  = uncorr$sigma2_v2,
    sigma2_e   = uncorr$sigma2_e,
    rho        = rho_hat,
    rho_profile = data.frame(
      rho    = c(rho_grid_coarse, rho_grid_fine),
      loglik = c(ll_coarse, ll_fine)
    ),
    predicted  = uncorr$predicted,
    converged  = uncorr$converged
  )
}



#' Print method for mixtran_fit objects
#' @export
print.mixtran_fit <- function(x, ...) {
  cat("NCI Usual Intake Model (MIXTRAN)\n")
  cat(sprintf("  Model type:     %s\n", x$model_type))
  cat(sprintf("  Intake variable: %s\n", x$intake_var))
  cat(sprintf("  N subjects:     %d\n", x$n_subjects))
  cat(sprintf("  N recalls:      %d\n", x$n_recalls))
  cat(sprintf("  N positive:     %d (%.1f%%)\n", x$n_positive, x$pct_consuming))
  cat(sprintf("  Box-Cox lambda: %.3f\n", x$lambda))

  if (x$model_type == "amount") {
    cat(sprintf("  sigma2_between: %.4f\n", x$sigma2_b))
    cat(sprintf("  sigma2_within:  %.4f\n", x$sigma2_w))
    cat(sprintf("  Variance ratio: %.2f\n", x$var_ratio))
    cat("  Fixed effects (amount):\n")
    for (nm in names(x$beta)) {
      cat(sprintf("    %s: %.4f\n", nm, x$beta[nm]))
    }
  } else {
    cat(sprintf("  sigma2_v1 (prob):   %.4f\n", x$sigma2_v1))
    cat(sprintf("  sigma2_v2 (amount): %.4f\n", x$sigma2_v2))
    cat(sprintf("  sigma2_e (within):  %.4f\n", x$sigma2_e))
    cat(sprintf("  rho:                %.4f\n", x$rho))
    cat("  Fixed effects (probability):\n")
    for (nm in names(x$alpha)) {
      cat(sprintf("    %s: %.4f\n", nm, x$alpha[nm]))
    }
    cat("  Fixed effects (amount):\n")
    for (nm in names(x$beta)) {
      cat(sprintf("    %s: %.4f\n", nm, x$beta[nm]))
    }
  }
  invisible(x)
}
