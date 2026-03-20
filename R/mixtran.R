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
#' @param prob_engine Engine for the logistic GLMM probability sub-model in
#'   two-part models. `"glmmPQL"` (default) uses `MASS::glmmPQL()` —
#'   fast but uses penalized quasi-likelihood (PQL), which can attenuate
#'   variance components. `"glmer"` uses `lme4::glmer()` with Laplace
#'   approximation for full ML, giving closer parity with SAS PROC NLMIXED.
#'   Requires the **lme4** package.
#' @param corr_engine Engine for the correlated two-part model
#'   (`model_type = "corr"` only). `"profile_rho"` (default) uses profile
#'   likelihood over a grid of ρ values — fast and robust, but approximate.
#'   `"ghq"` uses Gauss-Hermite quadrature to integrate out the bivariate
#'   random effects and optimises over (σ²_v1, σ²_v2, σ²_e, ρ) jointly,
#'   providing more accurate variance component and ρ estimates at the cost
#'   of additional computation. The number of quadrature nodes is controlled
#'   by `ghq_n_nodes`.
#' @param ghq_n_nodes Number of Gauss-Hermite quadrature nodes per dimension
#'   when `corr_engine = "ghq"`. Total bivariate nodes = `ghq_n_nodes^2`.
#'   Supported values: 3, 5 (default), 7, 9. Larger values are more accurate
#'   but slower. 5 is typically sufficient for variance component estimation.
#' @param start Optional `mixtran_fit` object from a previous call. When
#'   supplied, its fixed-effect estimates are used as starting values for the
#'   optimiser, which can improve convergence speed and stability (useful in
#'   BRR replicates or when re-fitting after a lambda update).
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
                    prob_engine = c("glmmPQL", "glmer"),
                    corr_engine = c("profile_rho", "ghq"),
                    ghq_n_nodes = 5L,
                    start = NULL,
                    verbose = TRUE) {

  prob_engine <- match.arg(prob_engine)
  corr_engine <- match.arg(corr_engine)
  if (!is.null(start) && !inherits(start, "mixtran_fit")) {
    stop("'start' must be a mixtran_fit object or NULL.")
  }

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
    result <- fit_amount_model(prep, lambda, verbose, start = start)
  } else if (model_type == "uncorr") {
    result <- fit_twopart_uncorr(prep, lambda, verbose,
                                  prob_engine = prob_engine, start = start)
  } else {
    if (corr_engine == "ghq") {
      result <- fit_twopart_corr_ghq(prep, lambda, verbose,
                                      prob_engine = prob_engine,
                                      n_nodes = ghq_n_nodes,
                                      start = start)
    } else {
      result <- fit_twopart_corr(prep, lambda, verbose,
                                  prob_engine = prob_engine, start = start)
    }
  }
  result$prob_engine  <- prob_engine
  result$corr_engine  <- if (model_type == "corr") corr_engine else NULL

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


#' Drop covariates that are constant within a data slice
#'
#' A covariate with fewer than 2 distinct non-NA values produces a singular
#' design matrix. This helper silently removes such columns and emits a
#' message() naming each dropped covariate so callers are aware.
#'
#' @param cov_names Character vector of candidate covariate names.
#' @param data Data frame to check values against.
#' @return Character vector of covariate names with constant columns removed.
#' @keywords internal
drop_constant_covs <- function(cov_names, data) {
  if (length(cov_names) == 0L) return(cov_names)
  keep <- vapply(cov_names, function(v) {
    vals <- data[[v]]
    length(unique(vals[!is.na(vals)])) >= 2L
  }, logical(1L))
  dropped <- cov_names[!keep]
  if (length(dropped) > 0L) {
    message("Dropping constant covariate(s) with fewer than 2 distinct non-NA values: ",
            paste(dropped, collapse = ", "))
  }
  cov_names[keep]
}


#' Run nlme::lme with optim; on Singularity/MEEM/backsolve error retry with nlminb
#'
#' Both the REML and ML nlme::lme() calls can throw uncaught MEEM() /
#' backsolve() errors when the data are near-singular and opt = "optim" is used.
#' This wrapper catches those specific errors and retries with opt = "nlminb".
#'
#' @param fixed,random,data,method,start,weights Arguments forwarded to nlme::lme().
#' @param verbose Logical; if TRUE emit a message when the nlminb retry fires.
#' @return An nlme "lme" fit object.
#' @keywords internal
.lme_robust <- function(fixed, random, data, method,
                        start = NULL, weights = NULL, verbose = FALSE) {
  # Check whether the installed nlme version supports the `start` parameter
  # (added in nlme >= 3.1-163)
  nlme_has_start <- "start" %in% names(formals(nlme::lme))
  use_start <- !is.null(start) && nlme_has_start

  # Helper: call lme() from a function whose enclosing environment is
  # parented by nlme's namespace.
  #
  # WHY THIS IS NECESSARY:
  # nlme's lme.formula() captures the call via match.call() and later
  # re-evaluates it with eval(Call, parent.frame()).  parent.frame()
  # returns the *call frame* of the caller — and R's eval(expr, envir)
  # does NOT change what parent.frame() returns.  So the previous approach
  # of eval(cl, envir = <nlme-parented env>) failed: parent.frame()
  # still returned call_lme's frame (in nciusual's namespace), which
  # lacks nlme's unexported helpers like base_lme.
  #
  # The fix: create an actual *function* whose enclosing environment is
  # parented by asNamespace("nlme"), then call it.  Now parent.frame()
  # from within lme.formula() → lme() generic → this wrapper function,
  # and symbol lookups walk: wrapper's frame → wrapper_env → nlme NS,
  # which contains all unexported helpers.
  call_lme <- function(ctl) {
    # Collect arguments into an environment parented by nlme's namespace
    args <- list(
      fixed   = fixed,
      random  = random,
      data    = data,
      method  = method,
      control = ctl
    )
    if (use_start)          args$start   <- start
    if (!is.null(weights))  args$weights <- weights

    wrapper_env <- list2env(args, parent = asNamespace("nlme"))

    # Build a function whose body is lme(fixed=fixed, random=random, ...)
    # and whose enclosing environment is wrapper_env.  Calling this function
    # makes parent.frame() from within lme.formula() resolve through
    # wrapper_env → asNamespace("nlme"), finding all internal helpers.
    arg_names <- names(args)
    call_expr <- as.call(c(
      quote(lme),
      stats::setNames(lapply(arg_names, as.symbol), arg_names)
    ))
    wrapper <- function() {}
    body(wrapper) <- call_expr
    environment(wrapper) <- wrapper_env
    wrapper()
  }

  ctl_optim <- nlme::lmeControl(
    maxIter = 500, msMaxIter = 500, opt = "optim", returnObject = TRUE
  )
  tryCatch(
    call_lme(ctl_optim),
    error = function(e) {
      if (grepl("Singularity|MEEM|backsolve", conditionMessage(e),
                ignore.case = TRUE)) {
        if (verbose) message(sprintf(
          "  %s Singularity/MEEM detected – retrying with opt = \"nlminb\"...",
          method))
        ctl_nlminb <- nlme::lmeControl(
          maxIter = 500, msMaxIter = 500, opt = "nlminb", returnObject = TRUE
        )
        tryCatch(
          call_lme(ctl_nlminb),
          error = function(e2) {
            stop(sprintf(
              "nlme::lme() failed with both optim and nlminb (method = %s): %s",
              method, conditionMessage(e2)
            ))
          }
        )
      } else {
        stop(e)
      }
    }
  )
}


#' Fit amount-only model (for ubiquitously consumed nutrients)
#'
#' Model: T(Y_ij) = X_ij'beta + u_i + eps_ij
#'   where u_i ~ N(0, sigma2_b) and eps_ij ~ N(0, sigma2_w)
#'
#' @keywords internal
fit_amount_model <- function(prep, lambda, verbose, start = NULL) {

  work <- prep$data
  cov_names <- prep$cov_names

  # Apply Box-Cox transform
  work$t_intake <- boxcox_transform(work$intake, lambda)

  # Drop covariates that are constant within this data slice (would produce a
  # singular design matrix and an uncaught error from nlme::lme).
  cov_names <- drop_constant_covs(cov_names, work)

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

  # Starting values from a prior fit improve convergence in BRR replicates
  lme_start <- NULL
  if (!is.null(start) && !is.null(start$beta)) {
    beta_names <- names(stats::model.matrix(fixed_formula, data = work[1, , drop = FALSE]))
    shared <- intersect(names(start$beta), beta_names)
    if (length(shared) > 0) {
      sv <- rep(0, length(beta_names))
      names(sv) <- beta_names
      sv[shared] <- start$beta[shared]
      lme_start <- list(fixed = sv)
    }
  }

  wts <- if (any(work$weight != 1)) nlme::varFixed(~ 1/weight) else NULL

  fit <- tryCatch({
    # Primary: REML with optim; retries with nlminb on Singularity/MEEM/backsolve
    .lme_robust(
      fixed = fixed_formula, random = ~ 1 | subject, data = work,
      method = "REML", start = lme_start, weights = wts, verbose = verbose
    )
  }, error = function(e) {
    # Fall back to ML if REML fails for any other reason
    if (verbose) message("  REML failed, trying ML...")
    tryCatch(
      .lme_robust(
        fixed = fixed_formula, random = ~ 1 | subject, data = work,
        method = "ML", start = NULL, weights = wts, verbose = verbose
      ),
      error = function(e2) {
        stop(sprintf(
          "nlme::lme() failed for both REML and ML in fit_amount_model(): %s",
          conditionMessage(e2)
        ))
      }
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
fit_twopart_uncorr <- function(prep, lambda, verbose,
                                prob_engine = "glmmPQL", start = NULL) {

  work <- prep$data
  cov_names <- prep$cov_names

  if (verbose) {
    eng_label <- if (prob_engine == "glmer") "lme4::glmer" else "MASS::glmmPQL"
    message(sprintf("Fitting two-part uncorrelated model (prob engine: %s)...", eng_label))
  }

  # --- Part 1: Probability model (logistic GLMM) ---
  if (verbose) message("  Part 1: Probability of consumption...")

  # Drop covariates that are constant in this data slice (all rows for prob model)
  prob_cov_names <- drop_constant_covs(cov_names, work)

  if (length(prob_cov_names) > 0) {
    prob_formula_fixed <- stats::as.formula(
      paste("consumed ~", paste(prob_cov_names, collapse = " + "))
    )
  } else {
    prob_formula_fixed <- stats::as.formula("consumed ~ 1")
  }

  prob_fit <- NULL
  alpha <- NULL
  sigma2_v1 <- NA

  if (prob_engine == "glmer") {
    # --- lme4::glmer: full ML via Laplace approximation ---
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("lme4 is required for prob_engine='glmer'. Install with: install.packages('lme4')")
    }
    glmer_formula <- stats::as.formula(
      paste("consumed ~",
            paste(if (length(prob_cov_names) > 0) prob_cov_names else "1", collapse = " + "),
            "+ (1 | subject)")
    )
    prob_fit <- tryCatch({
      lme4::glmer(
        formula = glmer_formula,
        family  = stats::binomial(link = "logit"),
        data    = work,
        control = lme4::glmerControl(
          optimizer = "bobyqa",
          optCtrl   = list(maxfun = 2e5)
        )
      )
    }, warning = function(w) {
      # glmer issues convergence warnings; catch and return the fit anyway
      withCallingHandlers(
        lme4::glmer(
          formula = glmer_formula,
          family  = stats::binomial(link = "logit"),
          data    = work,
          control = lme4::glmerControl(
            optimizer = "bobyqa",
            optCtrl   = list(maxfun = 2e5)
          )
        ),
        warning = function(w) invokeRestart("muffleWarning")
      )
    }, error = function(e) {
      if (verbose) message("    glmer failed: ", e$message,
                           "\n    Falling back to glmmPQL...")
      NULL
    })

    if (!is.null(prob_fit)) {
      alpha <- lme4::fixef(prob_fit)
      sigma2_v1 <- as.numeric(lme4::VarCorr(prob_fit)$subject[1, 1])
      attr(prob_fit, ".engine") <- "glmer"
      if (verbose) {
        message(sprintf("    Between-person variance (prob, glmer): %.4f", sigma2_v1))
      }
    }
    # Fall through to glmmPQL if glmer returned NULL
    if (is.null(prob_fit)) prob_engine <- "glmmPQL"
  }

  if (prob_engine == "glmmPQL") {
    # --- MASS::glmmPQL: penalized quasi-likelihood (faster, default) ---
    prob_fit <- tryCatch({
      MASS::glmmPQL(
        fixed   = prob_formula_fixed,
        random  = ~ 1 | subject,
        family  = stats::binomial(link = "logit"),
        data    = work,
        verbose = FALSE
      )
    }, error = function(e) {
      if (verbose) message("    glmmPQL failed: ", e$message)
      NULL
    })

    if (!is.null(prob_fit)) {
      attr(prob_fit, ".engine") <- "glmmPQL"
      alpha <- nlme::fixef(prob_fit)
      vc_prob <- nlme::VarCorr(prob_fit)
      sigma2_v1 <- as.numeric(vc_prob[1, "Variance"])
      if (verbose) {
        message(sprintf("    Between-person variance (prob, PQL): %.4f", sigma2_v1))
      }
    }
  }

  if (is.null(prob_fit)) {
    warning("Probability model failed to converge. Returning partial results.")
  }

  # --- Part 2: Amount model (conditional on consumption) ---
  if (verbose) message("  Part 2: Amount consumed (conditional)...")

  work_pos <- work[work$consumed == 1, ]
  work_pos$t_intake <- boxcox_transform(work_pos$intake, lambda)

  # Drop covariates that are constant in the positive-only slice
  amt_cov_names <- drop_constant_covs(cov_names, work_pos)

  if (length(amt_cov_names) > 0) {
    amt_formula <- stats::as.formula(
      paste("t_intake ~", paste(amt_cov_names, collapse = " + "))
    )
  } else {
    amt_formula <- stats::as.formula("t_intake ~ 1")
  }

  # Starting values for amount fixed effects from prior fit
  amt_start <- NULL
  if (!is.null(start) && !is.null(start$beta)) {
    beta_names <- names(stats::model.matrix(amt_formula, data = work_pos[1, , drop = FALSE]))
    shared <- intersect(names(start$beta), beta_names)
    if (length(shared) > 0) {
      sv <- rep(0, length(beta_names))
      names(sv) <- beta_names
      sv[shared] <- start$beta[shared]
      amt_start <- list(fixed = sv)
    }
  }

  amt_wts <- if (any(work_pos$weight != 1)) nlme::varFixed(~ 1/weight) else NULL

  amt_fit <- tryCatch({
    # Primary: REML with optim; retries with nlminb on Singularity/MEEM/backsolve
    .lme_robust(
      fixed = amt_formula, random = ~ 1 | subject, data = work_pos,
      method = "REML", start = amt_start, weights = amt_wts, verbose = verbose
    )
  }, error = function(e) {
    if (verbose) message("    Amount model REML failed, trying ML...")
    tryCatch(
      .lme_robust(
        fixed = amt_formula, random = ~ 1 | subject, data = work_pos,
        method = "ML", start = NULL, weights = amt_wts, verbose = verbose
      ),
      error = function(e2) {
        stop(sprintf(
          "nlme::lme() failed for both REML and ML in fit_twopart_uncorr() amount sub-model: %s",
          conditionMessage(e2)
        ))
      }
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

  # Predicted values (population-level, i.e. without subject random effects)
  pred_df <- work
  pred_df$prob_linpred <- if (!is.null(prob_fit)) {
    prob_popn_predict(prob_fit, newdata = work)
  } else {
    NA_real_
  }

  pred_df$amt_linpred <- NA_real_
  pred_df$amt_linpred[pred_df$consumed == 1] <- stats::predict(amt_fit, level = 0)

  list(
    prob_fit  = prob_fit,
    amt_fit   = amt_fit,
    alpha     = alpha,
    beta      = beta,
    sigma2_v1 = sigma2_v1,
    sigma2_v2 = sigma2_v2,
    sigma2_e  = sigma2_e,
    rho       = 0,  # uncorrelated by definition
    predicted = pred_df,
    converged = !is.null(prob_fit)
  )
}


# ---------------------------------------------------------------------------
# Engine-agnostic helpers for probability sub-model
# ---------------------------------------------------------------------------

#' Population-level predictions from prob model (no random effects)
#' Works for both glmmPQL (nlme) and glmer (lme4) fits.
#' @keywords internal
prob_popn_predict <- function(prob_fit, newdata) {
  if (inherits(prob_fit, "glmerMod")) {
    stats::predict(prob_fit, newdata = newdata, re.form = NA, type = "link")
  } else {
    stats::predict(prob_fit, level = 0, newdata = newdata)
  }
}

#' Subject-level random effects from prob model
#' Returns a named numeric vector (name = subject ID as character).
#' @keywords internal
prob_ranef <- function(prob_fit) {
  if (inherits(prob_fit, "glmerMod")) {
    re <- lme4::ranef(prob_fit)$subject  # data frame, 1 col per RE term
    setNames(as.numeric(re[[1]]), rownames(re))
  } else {
    re <- nlme::ranef(prob_fit)[[1]]     # data frame, 1 col ("(Intercept)")
    setNames(as.numeric(re[[1]]), rownames(re))
  }
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
fit_twopart_corr <- function(prep, lambda, verbose,
                              prob_engine = "glmmPQL", start = NULL) {

  if (verbose) message("Fitting two-part correlated model (profile likelihood over rho)...")

  # Step 1: Fit uncorrelated model to get all parameters
  if (verbose) message("  Step 1: Fitting uncorrelated model...")
  uncorr <- fit_twopart_uncorr(prep, lambda, verbose = FALSE,
                                prob_engine = prob_engine, start = start)

  if (!uncorr$converged || is.na(uncorr$sigma2_v1) || uncorr$sigma2_v1 <= 0) {
    if (!uncorr$converged) {
      warning("Uncorrelated model failed; falling back to rho = 0.")
    } else {
      warning("Probability variance (sigma2_v1) is NA or zero; falling back to rho = 0.")
    }
    uncorr$rho <- 0
    uncorr$rho_profile <- data.frame(rho = numeric(0), loglik = numeric(0))
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
  # We approximate v1_i by the empirical Bayes estimate, then standardise.
  # IMPORTANT: align by subject ID — amt_fit is fit only on consumed rows, so
  # some subjects (never-consumers) may be absent from amt_re. Misalignment
  # would corrupt the profile likelihood calculation.
  all_subjects <- as.character(sort(unique(work$subject)))

  prob_re_named <- tryCatch(
    prob_ranef(uncorr$prob_fit),
    error = function(e) setNames(rep(0, length(all_subjects)), all_subjects)
  )
  amt_re_named <- tryCatch({
    re_df <- nlme::ranef(uncorr$amt_fit)[[1]]
    setNames(as.numeric(re_df[[1]]), rownames(re_df))
  }, error = function(e) setNames(rep(0, length(all_subjects)), all_subjects))

  # Build aligned vectors (never-consumers get amt_re = 0)
  prob_re <- as.numeric(prob_re_named[all_subjects])
  amt_re  <- numeric(length(all_subjects))
  names(amt_re) <- all_subjects
  common_subjects <- intersect(all_subjects, names(amt_re_named))
  if (length(common_subjects) > 0) {
    amt_re[common_subjects] <- amt_re_named[common_subjects]
  }

  # Standardise to unit variance (guard against NA/NaN/zero)
  sd_v1 <- sqrt(uncorr$sigma2_v1)
  sd_v2 <- sqrt(uncorr$sigma2_v2)
  r1 <- if (isTRUE(sd_v1 > 0)) prob_re / sd_v1 else prob_re
  r2 <- if (isTRUE(sd_v2 > 0)) amt_re  / sd_v2 else amt_re

  # Replace any residual NAs (e.g. from subject lookup gaps) with 0 so that
  # the profile log-likelihood never receives NA inputs and cannot return NA.
  n_na_r1 <- sum(is.na(r1))
  n_na_r2 <- sum(is.na(r2))
  r1[is.na(r1)] <- 0
  r2[is.na(r2)] <- 0

  # Guard against degenerate residuals: when too many residuals are zero
  # (e.g. from NA replacement or failed random-effect extraction), the profile
  # log-likelihood reduces to -0.5*n*log(1-rho^2) which increases with |rho|,
  # causing the grid search to select a boundary value (rho ~ +-1).
  # Detect this and fall back to rho = 0 (conservative but correct).
  n_total <- length(r1)
  n_both_nonzero <- sum(r1 != 0 & r2 != 0)
  if (n_both_nonzero < 5 || n_na_r1 > 0.5 * n_total || n_na_r2 > 0.5 * n_total) {
    if (verbose) {
      message(sprintf(
        "  Residuals degenerate (non-zero pairs: %d/%d, NAs replaced: r1=%d, r2=%d); falling back to rho = 0.",
        n_both_nonzero, n_total, n_na_r1, n_na_r2
      ))
    }
    warning("Too few valid residual pairs for profile likelihood; falling back to rho = 0.")
    uncorr$rho <- 0
    uncorr$rho_profile <- data.frame(rho = numeric(0), loglik = numeric(0))
    return(uncorr)
  }

  # Profile log-likelihood over rho:
  # Given rho, the bivariate density of (v1, v2) changes. We evaluate how well
  # the standardised residuals (r1, r2) fit a BVN(0, [[1, rho],[rho, 1]]) model
  # using the bivariate normal log-likelihood.
  profile_loglik <- function(rho) {
    if (!is.finite(rho)) return(-Inf)
    n <- length(r1)
    det_val <- 1 - rho^2
    if (!is.finite(det_val) || det_val <= 0) return(-Inf)
    # Sum of bivariate normal log-densities for (r1_i, r2_i)
    ll <- -0.5 * n * log(det_val) -
      0.5 / det_val * sum(r1^2 - 2 * rho * r1 * r2 + r2^2) +
      # subtract the univariate terms already counted (constants)
      0.5 * sum(r1^2) + 0.5 * sum(r2^2)
    if (!is.finite(ll)) return(-Inf)
    ll
  }

  # Helper: safe which.max that falls back to index 1 on NA/empty result
  safe_which_max <- function(x) {
    idx <- which.max(x)
    if (length(idx) == 0L || is.na(idx)) 1L else idx
  }

  # Coarse grid search then refine around maximum
  rho_grid_coarse <- seq(-0.9, 0.9, by = 0.1)
  ll_coarse <- vapply(rho_grid_coarse, profile_loglik, numeric(1))
  best_coarse <- rho_grid_coarse[safe_which_max(ll_coarse)]

  # Fine grid within ±0.15 of coarse best
  rho_grid_fine <- seq(
    max(-0.99, best_coarse - 0.15),
    min( 0.99, best_coarse + 0.15),
    by = 0.01
  )
  ll_fine <- vapply(rho_grid_fine, profile_loglik, numeric(1))
  rho_hat <- rho_grid_fine[safe_which_max(ll_fine)]

  # Sanity check: if rho_hat is at the grid boundary, the profile likelihood
  # likely degenerated (e.g. insufficient residual signal). Cross-check with
  # the sample correlation of the residuals.
  if (abs(rho_hat) >= 0.98) {
    r_sample <- tryCatch(stats::cor(r1, r2), error = function(e) NA_real_)
    if (is.na(r_sample) || sign(r_sample) != sign(rho_hat) || abs(r_sample) < 0.3) {
      if (verbose) {
        message(sprintf(
          "    Profile rho (%.4f) at grid boundary; sample cor = %.4f. Using sample cor.",
          rho_hat, ifelse(is.na(r_sample), 0, r_sample)
        ))
      }
      rho_hat <- if (is.na(r_sample)) 0 else max(-0.95, min(0.95, r_sample))
    }
  }

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


# ---------------------------------------------------------------------------
# GHQ engine for correlated two-part model (v0.2)
# ---------------------------------------------------------------------------

#' Fit two-part correlated model via Gauss-Hermite Quadrature
#'
#' Integrates out the bivariate random effects (v1_i, v2_i) using tensor-
#' product Gauss-Hermite quadrature and optimises the marginal log-likelihood
#' over the variance components (sigma2_v1, sigma2_v2, sigma2_e) and the
#' correlation rho.
#'
#' Fixed effects (alpha for the probability sub-model and beta for the amount
#' sub-model) are held at the estimates from the uncorrelated model; this
#' mirrors the approach taken by the profile likelihood engine and is
#' substantially faster than full joint optimisation over all parameters.
#'
#' @param prep Prepared data list from prepare_mixtran_data()
#' @param lambda Box-Cox lambda (already determined)
#' @param verbose Print progress
#' @param prob_engine "glmmPQL" or "glmer" for the uncorr initialisation step
#' @param n_nodes Number of GH nodes per dimension (default 5)
#' @param start Optional mixtran_fit for warm-starting variance components
#' @keywords internal
fit_twopart_corr_ghq <- function(prep, lambda, verbose,
                                  prob_engine = "glmmPQL",
                                  n_nodes = 5L,
                                  start = NULL) {

  if (verbose) message(sprintf(
    "Fitting two-part correlated model (GHQ, n_nodes=%d)...", n_nodes))

  # --- Step 1: starting values from uncorrelated model ---
  if (verbose) message("  Step 1: Fitting uncorrelated model for starting values...")
  uncorr <- fit_twopart_uncorr(prep, lambda, verbose = FALSE,
                                prob_engine = prob_engine, start = start)

  if (!uncorr$converged) {
    warning("Uncorrelated model did not converge; falling back to profile_rho.")
    return(fit_twopart_corr(prep, lambda, verbose,
                             prob_engine = prob_engine, start = start))
  }

  work      <- prep$data
  cov_names <- prep$cov_names

  # Population-level linear predictors (no random effects)
  eta_pop <- as.numeric(prob_popn_predict(uncorr$prob_fit, newdata = work))

  work_pos     <- work[work$consumed == 1, ]
  mu_pop       <- as.numeric(stats::predict(uncorr$amt_fit, level = 0))
  y_pos        <- work_pos$intake
  t_y_pos      <- boxcox_transform(y_pos, lambda)
  log_y_pos    <- log(y_pos)

  # Integer subject indices (1 : n_subj) for fast tapply-style aggregation
  subjects     <- unique(work$subject)
  subj_map     <- stats::setNames(seq_along(subjects), as.character(subjects))
  person_id    <- as.integer(subj_map[as.character(work$subject)])
  person_id_pos <- as.integer(subj_map[as.character(work_pos$subject)])
  n_subj       <- prep$n_subjects

  consumed <- work$consumed

  # --- Step 2: marginal log-likelihood via GHQ ---
  ghq_loglik <- function(log_sv1, log_sv2, log_se, atanh_rho) {
    sigma_v1 <- exp(log_sv1)
    sigma_v2 <- exp(log_sv2)
    sigma_e  <- exp(log_se)
    rho      <- tanh(atanh_rho)

    nodes_biv <- gh_nodes_bivariate(n_nodes, rho, sigma_v1, sigma_v2)
    n_q       <- nrow(nodes_biv)

    # log-lik contribution matrix: n_subj × n_q
    ll_mat <- matrix(0, nrow = n_subj, ncol = n_q)

    for (q in seq_len(n_q)) {
      v1_q <- nodes_biv$v1[q]
      v2_q <- nodes_biv$v2[q]

      # Probability contributions (all observations)
      eta_q <- eta_pop + v1_q
      log_p <- ifelse(consumed == 1L,
                      stats::plogis(eta_q, log.p = TRUE),
                      stats::plogis(eta_q, log.p = TRUE, lower.tail = FALSE))
      log_prob_subj <- as.numeric(tapply(log_p, person_id, sum))

      # Amount contributions (positive observations only)
      log_amt_subj <- numeric(n_subj)
      if (length(t_y_pos) > 0L) {
        log_f <- stats::dnorm(t_y_pos, mean = mu_pop + v2_q,
                              sd = sigma_e, log = TRUE) +
                 (lambda - 1) * log_y_pos
        amt_by_subj <- tapply(log_f, person_id_pos, sum)
        log_amt_subj[as.integer(names(amt_by_subj))] <-
          as.numeric(amt_by_subj)
      }

      ll_mat[, q] <- log_prob_subj + log_amt_subj
    }

    # Log-sum-exp over quadrature points: log L_i = log Σ_q w_q exp(ll_iq)
    lw <- nodes_biv$log_weight   # length n_q
    ll_plus_lw <- ll_mat + matrix(lw, nrow = n_subj, ncol = n_q, byrow = TRUE)
    row_max    <- apply(ll_plus_lw, 1, max)
    log_marg   <- row_max + log(rowSums(exp(ll_plus_lw - row_max)))
    sum(log_marg)
  }

  # --- Step 3: optimise variance components + rho ---
  if (verbose) message("  Step 2: Optimising variance components via GHQ...")

  # Starting values from uncorr fit (or prior start if supplied)
  if (!is.null(start) && inherits(start, "mixtran_fit") &&
      !is.null(start$sigma2_v1)) {
    sv1_init <- log(sqrt(start$sigma2_v1))
    sv2_init <- log(sqrt(start$sigma2_v2))
    se_init  <- log(sqrt(start$sigma2_e))
    rho_init <- atanh(min(max(start$rho, -0.95), 0.95))
  } else {
    sv1_init <- log(sqrt(uncorr$sigma2_v1))
    sv2_init <- log(sqrt(uncorr$sigma2_v2))
    se_init  <- log(sqrt(uncorr$sigma2_e))
    rho_init <- atanh(0)
  }

  neg_ll <- function(par) {
    tryCatch(
      -ghq_loglik(par[1], par[2], par[3], par[4]),
      error = function(e) 1e10
    )
  }

  opt <- tryCatch(
    stats::optim(
      par     = c(sv1_init, sv2_init, se_init, rho_init),
      fn      = neg_ll,
      method  = "L-BFGS-B",
      lower   = c(-6, -6, -6, -4),
      upper   = c( 6,  6,  6,  4),
      control = list(maxit = 300, factr = 1e8)
    ),
    error = function(e) {
      if (verbose) message("    L-BFGS-B failed (", e$message,
                           "); retrying with Nelder-Mead...")
      tryCatch(
        stats::optim(
          par     = c(sv1_init, sv2_init, se_init, rho_init),
          fn      = neg_ll,
          method  = "Nelder-Mead",
          control = list(maxit = 1000, reltol = 1e-6)
        ),
        error = function(e2) NULL
      )
    }
  )

  if (is.null(opt) || opt$convergence > 1) {
    if (verbose) message("    GHQ optimisation did not converge; ",
                         "falling back to profile_rho.")
    return(fit_twopart_corr(prep, lambda, verbose,
                             prob_engine = prob_engine, start = start))
  }

  p       <- opt$par
  sigma_v1 <- exp(p[1]);  sigma_v2 <- exp(p[2])
  sigma_e  <- exp(p[3]);  rho_hat  <- tanh(p[4])

  if (verbose) {
    message(sprintf("    sigma_v1 = %.4f  sigma_v2 = %.4f  sigma_e = %.4f  rho = %.4f",
                    sigma_v1, sigma_v2, sigma_e, rho_hat))
    message(sprintf("    GHQ log-likelihood = %.3f  (convergence code: %d)",
                    -opt$value, opt$convergence))
  }

  list(
    prob_fit    = uncorr$prob_fit,
    amt_fit     = uncorr$amt_fit,
    alpha       = uncorr$alpha,
    beta        = uncorr$beta,
    sigma2_v1   = sigma_v1^2,
    sigma2_v2   = sigma_v2^2,
    sigma2_e    = sigma_e^2,
    rho         = rho_hat,
    ghq_n_nodes = n_nodes,
    ghq_loglik  = -opt$value,
    predicted   = uncorr$predicted,
    converged   = (opt$convergence == 0)
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
