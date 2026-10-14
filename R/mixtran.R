#' MIXTRAN: Fit NCI Usual Intake Models
#'
#' R port of the NCI MIXTRAN SAS macro (v2.1). Fits nonlinear mixed models
#' to 24-hour dietary recall data to estimate the parameters of the usual
#' intake distribution, separating within-person from between-person
#' variation.
#'
#' Three model types are supported, matching MIXTRAN's `modeltype=`:
#' - "amount": amount-only model for ubiquitously consumed components
#'             (MIXTRAN `modeltype=amount`)
#' - "uncorr": two-part model with uncorrelated random effects
#'             (MIXTRAN `modeltype=nocorr`)
#' - "corr":   two-part model with correlated random effects
#'             (MIXTRAN `modeltype=corr`, the default in SAS)
#'
#' @name mixtran
NULL

#' Fit an NCI usual intake model (MIXTRAN equivalent)
#'
#' All models are fitted the way MIXTRAN fits them with PROC NLMIXED:
#'
#' * **Survey weights** (`weight_var`, MIXTRAN's `replicate_var`) multiply each
#'   subject's marginal log-likelihood, as the NLMIXED `REPLICATE` statement
#'   does. The weight must be constant within subject. (SAS additionally
#'   requires integer weights; non-integer weights are accepted here.)
#' * **Box-Cox lambda** is a free parameter of the amount model, estimated by
#'   maximum likelihood with the bound `lambda >= 0.01` (`A_LAMBDA > 0.01`).
#'   Its starting value comes from a grid search over `lambda_grid`
#'   (0.05 to 1 by 0.05, as in MIXTRAN's PROC TRANSREG step) using the amount
#'   covariates, or 0.4 when there are none. Supply `lambda` to fix it.
#' * **Two-part models**: the uncorrelated model fits the probability and
#'   amount parts separately; the correlated model then maximises the joint
#'   likelihood starting from the uncorrelated estimates, with the random
#'   effects integrated out (adaptive Gauss-Hermite quadrature over the
#'   probability random effect; the amount random effect is integrated in
#'   closed form).
#'
#' @param data Data frame in long format (one row per recall per person)
#' @param intake_var Name of the intake variable (character)
#' @param subject_var Name of the subject ID variable (character)
#' @param repeat_var Name of the recall day indicator (character), e.g. 1 or 2
#' @param model_type One of "amount", "uncorr", "corr"
#' @param covariates Character vector of numeric covariate names for the
#'   amount part (MIXTRAN `covars_amt`) and, unless `covariates_prob` is given,
#'   also for the probability part. Categorical variables must be dummy-coded.
#' @param subgroup_var Optional name of a subgroup variable (character) to
#'   carry through to the subject-level predicted data for use by `distrib()`.
#'   It is not used in model fitting.
#' @param weekend_var Optional name of a 0/1 weekend (Fri-Sun) indicator. As in
#'   MIXTRAN it enters both parts of the model; `distrib()` then averages
#'   weekday and weekend usual intake (default weights 4/7 and 3/7).
#' @param weight_var Optional name of the survey weight variable (MIXTRAN
#'   `replicate_var`). Used in the likelihood and in `distrib()`.
#' @param lambda Box-Cox lambda. If `NULL` (default) it is estimated by
#'   maximum likelihood. A supplied value (including 0 for log) is held fixed.
#' @param lambda_grid Grid used to find the starting value for lambda.
#' @param min_positive Value used to replace zero intakes in the amount-only
#'   model. Default `NULL` uses half the smallest positive intake (MIXTRAN).
#' @param covariates_prob Covariates for the probability part (MIXTRAN
#'   `covars_prob`). Defaults to `covariates`.
#' @param include_seq Logical (default `TRUE`). Add a sequence indicator
#'   `seq_num` (1 for recalls after the first) to both parts, equivalent to
#'   passing `seq=day2` to MIXTRAN. Its effect is removed in `distrib()`.
#' @param n_quad Number of adaptive Gauss-Hermite points used to integrate
#'   over the probability random effect (two-part models). Default 21.
#' @param start Optional `mixtran_fit` from a previous call (e.g. the
#'   full-sample fit in BRR). Its parameter estimates are used as starting
#'   values (MIXTRAN `start_val1-3`).
#' @param skip_if_empty Logical (default `TRUE`). When the data are
#'   insufficient to fit the requested model, warn and return `NULL`; when
#'   `FALSE`, stop with an error instead.
#' @param verbose Print model fitting progress
#' @param prob_engine,corr_engine,ghq_n_nodes Deprecated and ignored. Earlier
#'   versions offered approximate engines (glmmPQL, profile likelihood over
#'   rho, non-adaptive GHQ); all models are now fitted by weighted maximum
#'   likelihood as in SAS.
#' @return A `mixtran_fit` object, or `NULL` (with a warning) when
#'   `skip_if_empty = TRUE` and the data are insufficient.
#' @export
mixtran <- function(data,
                    intake_var,
                    subject_var,
                    repeat_var,
                    model_type = c("amount", "uncorr", "corr"),
                    covariates = NULL,
                    subgroup_var = NULL,
                    weekend_var = NULL,
                    weight_var = NULL,
                    lambda = NULL,
                    lambda_grid = seq(0.05, 1.0, by = 0.05),
                    min_positive = NULL,
                    covariates_prob = covariates,
                    include_seq = TRUE,
                    n_quad = 21L,
                    start = NULL,
                    skip_if_empty = TRUE,
                    verbose = TRUE,
                    prob_engine = NULL,
                    corr_engine = NULL,
                    ghq_n_nodes = NULL) {

  if (!is.null(prob_engine) || !is.null(corr_engine) || !is.null(ghq_n_nodes)) {
    warning("'prob_engine', 'corr_engine' and 'ghq_n_nodes' are deprecated and ",
            "ignored: all models are fitted by weighted maximum likelihood, ",
            "as MIXTRAN does with PROC NLMIXED.", call. = FALSE)
  }
  if (!is.null(start) && !inherits(start, "mixtran_fit")) {
    stop("'start' must be a mixtran_fit object or NULL.")
  }
  model_type <- match.arg(model_type)
  data <- as.data.frame(data)

  # --- Validate inputs ---
  required_vars <- unique(c(intake_var, subject_var, repeat_var, weekend_var,
                            weight_var, covariates, covariates_prob, subgroup_var))
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  for (v in unique(c(covariates, covariates_prob, weekend_var, weight_var))) {
    if (!is.numeric(data[[v]]) && !is.logical(data[[v]])) {
      stop("Variable '", v, "' must be numeric (dummy-code categorical covariates).")
    }
  }
  reserved <- c("seq_num", "weekend", "subject", "intake", "weight", "consumed")
  if (any(c(covariates, covariates_prob) %in% reserved)) {
    stop("Covariate names may not be any of: ", paste(reserved, collapse = ", "))
  }

  insufficient <- function(msg) {
    if (isTRUE(skip_if_empty)) {
      warning(msg, " Returning NULL.", call. = FALSE)
      return(NULL)
    }
    stop(msg, call. = FALSE)
  }

  # --- Early exit when data are insufficient (MIXTRAN stops in these cases) ---
  vals <- data[[intake_var]]
  if (model_type == "amount") {
    if (!any(vals > 0, na.rm = TRUE)) {
      return(insufficient(sprintf(
        "Skipping '%s' (model_type = \"amount\"): no positive intake values found.",
        intake_var)))
    }
    n_rep <- sum(tapply(!is.na(vals), data[[subject_var]], sum) >= 2)
    if (n_rep < 2L) {
      return(insufficient(sprintf(
        "Skipping '%s': fewer than two subjects have repeat recalls (MIXTRAN stops in this case).",
        intake_var)))
    }
  } else {
    pos_recalls  <- tapply(vals > 0 & !is.na(vals), data[[subject_var]], sum)
    n_paired_pos <- sum(pos_recalls >= 2, na.rm = TRUE)
    if (n_paired_pos < 2L) {
      return(insufficient(sprintf(
        "Skipping '%s' (model_type = \"%s\"): fewer than two subjects have 2+ positive recalls.",
        intake_var, model_type)))
    }
    if (n_paired_pos <= 10L) {
      warning(sprintf(
        "Only %d subject(s) have 2+ positive recalls for '%s'; estimates might be unstable (MIXTRAN asks for at least 10).",
        n_paired_pos, intake_var), call. = FALSE)
    }
    pct_pos <- mean(vals > 0, na.rm = TRUE)
    if (pct_pos > 0.9) {
      message(sprintf(
        "Note: %.0f%% of person-days have positive '%s' intake; model_type = \"amount\" may be more suitable.",
        pct_pos * 100, intake_var))
    }
  }

  # --- Prepare data ---
  prep <- prepare_mixtran_data(
    data = data, intake_var = intake_var, subject_var = subject_var,
    repeat_var = repeat_var, covariates = covariates,
    covariates_prob = covariates_prob, weekend_var = weekend_var,
    weight_var = weight_var, min_positive = min_positive,
    model_type = model_type, include_seq = include_seq, verbose = verbose
  )
  if (is.null(prep)) return(insufficient(sprintf(
    "Skipping '%s': no usable records after removing missing values.", intake_var)))

  result <- fit_nci_model(prep, model_type = model_type, lambda = lambda,
                          lambda_grid = lambda_grid, n_quad = n_quad,
                          start = start, verbose = verbose)

  # Carry the subgroup variable through to the subject-level predictions
  if (!is.null(subgroup_var)) {
    lk <- data[!duplicated(data[[subject_var]]), c(subject_var, subgroup_var), drop = FALSE]
    n_distinct <- tapply(data[[subgroup_var]], data[[subject_var]],
                         function(z) length(unique(z)))
    if (any(n_distinct > 1, na.rm = TRUE)) {
      warning("'", subgroup_var, "' varies within subject; using each subject's first value.",
              call. = FALSE)
    }
    result$predicted[[subgroup_var]] <-
      lk[[subgroup_var]][match(result$predicted$subject, lk[[subject_var]])]
  }

  result$model_type    <- model_type
  result$intake_var    <- intake_var
  result$subject_var   <- subject_var
  result$covariates    <- covariates
  result$covariates_prob <- if (model_type == "amount") NULL else covariates_prob
  result$subgroup_var  <- subgroup_var
  result$weekend_var   <- weekend_var
  result$weight_var    <- weight_var
  result$n_subjects    <- prep$n_subjects
  result$n_recalls     <- nrow(prep$data)
  result$n_positive    <- prep$n_positive_raw
  result$pct_consuming <- prep$n_positive_raw / nrow(prep$data) * 100
  result$min_amt       <- prep$min_amt
  result$min_pos_raw   <- prep$min_pos_raw

  class(result) <- "mixtran_fit"
  result
}


#' Prepare data for MIXTRAN model fitting
#' @keywords internal
prepare_mixtran_data <- function(data, intake_var, subject_var, repeat_var,
                                 covariates, weekend_var, weight_var,
                                 min_positive, model_type,
                                 covariates_prob = covariates,
                                 include_seq = TRUE, verbose = TRUE) {

  work <- data.frame(
    subject    = data[[subject_var]],
    repeat_num = data[[repeat_var]],
    intake     = data[[intake_var]],
    stringsAsFactors = FALSE
  )
  if (!is.null(weekend_var)) work$weekend <- as.numeric(data[[weekend_var]])
  for (cv in unique(c(covariates, covariates_prob))) work[[cv]] <- as.numeric(data[[cv]])
  work$weight <- if (!is.null(weight_var)) as.numeric(data[[weight_var]]) else 1

  # NLMIXED drops records with missing values in any model variable
  model_cols <- intersect(c("intake", "repeat_num", "weekend", "weight",
                            covariates, covariates_prob), names(work))
  keep <- stats::complete.cases(work[, model_cols, drop = FALSE])
  if (any(!keep) && verbose) {
    message(sprintf("  Removing %d record(s) with missing intake, weight or covariates.",
                    sum(!keep)))
  }
  work <- work[keep, , drop = FALSE]
  if (nrow(work) == 0L) return(NULL)
  if (any(work$intake < 0)) stop("Intake values must be non-negative.")
  if (any(work$weight < 0)) stop("Weights must be non-negative.")

  # Sort by subject and repeat (MIXTRAN: proc sort by subject repeat)
  work <- work[order(work$subject, work$repeat_num), , drop = FALSE]
  rownames(work) <- NULL

  wt_var <- tapply(work$weight, work$subject, function(z) length(unique(z)))
  if (any(wt_var > 1)) {
    warning("The weight varies within subject; each subject's first-record weight ",
            "is used (as with the NLMIXED REPLICATE statement).", call. = FALSE)
  }

  # Sequence indicator (MIXTRAN seq=): 0 for the first recall, 1 afterwards
  if (include_seq) work$seq_num <- as.numeric(work$repeat_num > 1)

  n_positive_raw <- sum(work$intake > 0)
  min_pos <- if (n_positive_raw > 0) min(work$intake[work$intake > 0]) else NA_real_

  # Amount-only model: zeros -> half the smallest positive amount (MIXTRAN)
  if (model_type == "amount" && any(work$intake == 0)) {
    if (is.null(min_positive)) min_positive <- min_pos / 2
    if (verbose) message(sprintf("  Replacing %d zero values with %.6g (half-min)",
                                 sum(work$intake == 0), min_positive))
    work$intake[work$intake == 0] <- min_positive
  }
  work$consumed <- as.numeric(work$intake > 0)

  # min_amt passed to DISTRIB: PROC UNIVARIATE minimum of the response in the
  # MIXTRAN working data (zeros already replaced for amount models; zeros set
  # to missing for two-part models).
  min_amt <- min(work$intake[work$intake > 0])

  base_covs <- c(if (include_seq) "seq_num", if (!is.null(weekend_var)) "weekend")
  list(
    data = work,
    cov_names = c(base_covs, covariates),
    cov_names_prob = c(base_covs, covariates_prob),
    n_subjects = length(unique(work$subject)),
    n_positive_raw = n_positive_raw,
    min_positive = min_positive,
    min_amt = min_amt,
    min_pos_raw = min_pos
  )
}


#' Drop covariates that are constant within a data slice
#'
#' A covariate with fewer than 2 distinct non-NA values produces a singular
#' design matrix. This helper removes such columns and emits a message.
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


#' Build the design list used by the likelihood engine
#' @keywords internal
.nci_design <- function(work, prob_covs, amt_covs, lambda_fixed, n_quad,
                        has_prob = !is.null(prob_covs)) {
  id  <- match(work$subject, unique(work$subject))
  N   <- max(id)
  pos <- work$consumed == 1
  X2  <- cbind(1, as.matrix(work[pos, amt_covs, drop = FALSE]))
  colnames(X2) <- c("(Intercept)", as.character(amt_covs))
  X1 <- NULL
  if (has_prob) {
    X1 <- cbind(1, as.matrix(work[, prob_covs, drop = FALSE]))
    colnames(X1) <- c("(Intercept)", as.character(prob_covs))
  }
  w <- work$weight[!duplicated(id)]
  list(
    N = N, id = id, w = w, sum_w = sum(w),
    yn = work$consumed, sgn = 2 * work$consumed - 1,
    X1 = X1, X2 = X2,
    y_pos = work$intake[pos], logy_pos = log(work$intake[pos]),
    id_pos = id[pos], uid_pos = sort(unique(id[pos])),
    lambda_fixed = lambda_fixed, n_quad = as.integer(n_quad)
  )
}


#' Fit the MIXTRAN model by weighted maximum likelihood
#' @keywords internal
fit_nci_model <- function(prep, model_type, lambda = NULL,
                          lambda_grid = seq(0.05, 1, by = 0.05),
                          n_quad = 21L, start = NULL, verbose = TRUE) {
  work <- prep$data
  pos  <- work$consumed == 1

  amt_covs  <- drop_constant_covs(prep$cov_names, work[pos, , drop = FALSE])
  prob_covs <- if (model_type == "amount") NULL else
    drop_constant_covs(prep$cov_names_prob, work)

  d <- .nci_design(work, prob_covs, amt_covs,
                   lambda_fixed = if (is.null(lambda)) NA_real_ else lambda,
                   n_quad = n_quad, has_prob = model_type != "amount")
  row_w <- d$w[d$id] / mean(d$w)   # normalised: GLM/LM start fits only

  # --- Starting values (MIXTRAN: TRANSREG grid for lambda, GENMOD fits) ---
  y_pos <- d$y_pos
  if (is.null(lambda)) {
    lam0 <- if (length(amt_covs) > 0) {
      find_optimal_lambda(y_pos, lambda_grid = lambda_grid,
                          covariates = as.matrix(work[pos, amt_covs, drop = FALSE]),
                          use_log_rule = FALSE)$lambda
    } else 0.4
  } else lam0 <- lambda
  z0  <- boxcox_transform(y_pos, lam0)
  lmf <- stats::lm.wfit(d$X2, z0, w = row_w[pos])
  v0  <- sum(row_w[pos] * lmf$residuals^2) / sum(row_w[pos])
  p_amt <- c(stats::setNames(lmf$coefficients, paste0("A.", colnames(d$X2))),
             A_LOGSDE = 0.5 * log(v0 / 2), A_LOGSDU2 = 0.5 * log(v0 / 2))
  p_amt[is.na(p_amt)] <- 0
  if (is.null(lambda)) p_amt["A_LAMBDA"] <- lam0

  p_prob <- NULL
  if (model_type != "amount") {
    gf <- suppressWarnings(stats::glm.fit(d$X1, d$yn, weights = row_w,
                                          family = stats::binomial()))
    p_prob <- c(stats::setNames(gf$coefficients, paste0("P.", colnames(d$X1))),
                P_LOGSDU1 = 0.5)
    p_prob[is.na(p_prob)] <- 0
  }

  use_start <- function(p0) {
    if (!is.null(start) && !is.null(start$theta)) {
      sh <- intersect(names(p0), names(start$theta))
      p0[sh] <- start$theta[sh]
    }
    p0
  }

  if (verbose) {
    message(sprintf("Fitting %s model by weighted maximum likelihood (%d subjects)...",
                    model_type, d$N))
  }

  # --- Amount part (amount-only model, or the amount part of the two-part model) ---
  fit_a <- .nlmixed_optim(use_start(p_amt), d, "amount", verbose)
  theta <- fit_a$par
  conv  <- fit_a$converged
  m2ll  <- 2 * d$sum_w * fit_a$value

  if (model_type != "amount") {
    fit_p <- .nlmixed_optim(use_start(p_prob), d, "prob", verbose)
    theta <- c(fit_p$par, theta)
    conv  <- conv && fit_p$converged
    m2ll  <- m2ll + 2 * d$sum_w * fit_p$value
    theta_uncorr <- theta
    m2ll_uncorr  <- m2ll
  }

  if (model_type == "corr") {
    p_corr <- use_start(c(theta, Z_U = 0))
    fit_c  <- .nlmixed_optim(p_corr, d, "corr", verbose)
    theta  <- fit_c$par
    conv   <- fit_c$converged
    m2ll   <- 2 * d$sum_w * fit_c$value
  }

  q <- .unpack(theta, d)
  alpha <- if (model_type == "amount") NULL else
    stats::setNames(as.numeric(q$alpha), colnames(d$X1))
  beta  <- stats::setNames(as.numeric(q$beta), colnames(d$X2))
  lam   <- q$lambda

  # --- Empirical Bayes estimates of the random effects ---
  if (model_type == "amount") {
    re <- .ebe(theta, d, "amount")
  } else if (model_type == "uncorr") {
    re <- list(u1 = .ebe(theta, d, "prob")$u1,
               u2 = { r <- .ebe(theta, d, "amount")$u2
                      r[!(seq_len(d$N) %in% d$uid_pos)] <- 0; r })
  } else {
    re <- .ebe(theta, d, "corr")
  }

  # --- Subject-level predictions for DISTRIB (MIXTRAN _pred data set) ---
  # Linear predictors from each subject's first record with the sequence
  # effect removed (seq = 0); with a weekend variable, both weekday (_0, the
  # default columns) and weekend (_1) versions, as MIXTRAN's x1b1_0/x1b1_1.
  first <- work[!duplicated(work$subject), , drop = FALSE]
  xb <- function(coef, covs, wk) {
    if (is.null(coef)) return(rep(NA_real_, nrow(first)))
    df <- first
    if ("seq_num" %in% names(df)) df$seq_num <- 0
    if (!is.null(wk) && "weekend" %in% names(df)) df$weekend <- wk
    as.numeric(cbind(1, as.matrix(df[, covs, drop = FALSE])) %*% coef)
  }
  has_wk <- "weekend" %in% c(amt_covs, prob_covs)
  predicted <- data.frame(subject = first$subject, weight = first$weight,
                          stringsAsFactors = FALSE)
  predicted$x1b1 <- xb(alpha, prob_covs, if (has_wk) 0 else NULL)
  predicted$x2b2 <- xb(beta,  amt_covs,  if (has_wk) 0 else NULL)
  if (has_wk) {
    predicted$x1b1_1 <- xb(alpha, prob_covs, 1)
    predicted$x2b2_1 <- xb(beta,  amt_covs,  1)
  }
  predicted$u1 <- re$u1
  predicted$u2 <- re$u2

  if (verbose) {
    message(sprintf("  lambda = %.4f%s", lam, if (is.null(lambda)) " (estimated)" else " (fixed)"))
    if (model_type == "amount") {
      message(sprintf("  var_u2 (between) = %.4f   var_e (within) = %.4f", q$su2, q$se2))
    } else {
      message(sprintf("  var_u1 = %.4f   var_u2 = %.4f   var_e = %.4f   rho = %.4f",
                      q$sig1^2, q$su2, q$se2, q$rho))
    }
    message(sprintf("  -2 log likelihood (weighted) = %.6g", m2ll))
  }
  if (!conv) {
    warning("The likelihood maximisation did not fully converge; check the estimates.",
            call. = FALSE)
  }

  out <- list(
    lambda    = lam,
    lambda_estimated = is.null(lambda),
    beta      = beta,
    alpha     = alpha,
    theta     = theta,
    m2loglik  = m2ll,
    converged = conv,
    predicted = predicted,
    data      = work,
    prob_covs_used = prob_covs,
    amt_covs_used  = amt_covs,
    n_quad    = n_quad
  )
  if (model_type == "amount") {
    out$sigma2_b  <- q$su2
    out$sigma2_w  <- q$se2
    out$var_ratio <- q$se2 / q$su2
  } else {
    out$sigma2_v1 <- q$sig1^2
    out$sigma2_v2 <- q$su2
    out$sigma2_e  <- q$se2
    out$rho       <- q$rho
    out$cov_u1u2  <- q$rho * q$sig1 * sqrt(q$su2)
    out$theta_uncorr <- theta_uncorr
    out$m2loglik_uncorr <- m2ll_uncorr
  }
  out
}


#' Print method for mixtran_fit objects
#' @export
print.mixtran_fit <- function(x, ...) {
  cat("NCI Usual Intake Model (MIXTRAN)\n")
  cat(sprintf("  Model type:      %s\n", x$model_type))
  cat(sprintf("  Intake variable: %s\n", x$intake_var))
  cat(sprintf("  N subjects:      %d\n", x$n_subjects))
  cat(sprintf("  N recalls:       %d\n", x$n_recalls))
  cat(sprintf("  N positive:      %d (%.1f%%)\n", x$n_positive, x$pct_consuming))
  cat(sprintf("  Box-Cox lambda:  %.4f%s\n", x$lambda,
              if (isTRUE(x$lambda_estimated)) " (estimated)" else " (fixed)"))
  cat(sprintf("  Converged:       %s\n", x$converged))

  if (x$model_type == "amount") {
    cat(sprintf("  var_u2 (between): %.4f\n", x$sigma2_b))
    cat(sprintf("  var_e  (within):  %.4f\n", x$sigma2_w))
    cat("  Fixed effects (amount):\n")
    for (nm in names(x$beta)) cat(sprintf("    %s: %.4f\n", nm, x$beta[nm]))
  } else {
    cat(sprintf("  var_u1 (prob):    %.4f\n", x$sigma2_v1))
    cat(sprintf("  var_u2 (amount):  %.4f\n", x$sigma2_v2))
    cat(sprintf("  var_e (within):   %.4f\n", x$sigma2_e))
    cat(sprintf("  rho:              %.4f\n", x$rho))
    cat("  Fixed effects (probability):\n")
    for (nm in names(x$alpha)) cat(sprintf("    %s: %.4f\n", nm, x$alpha[nm]))
    cat("  Fixed effects (amount):\n")
    for (nm in names(x$beta)) cat(sprintf("    %s: %.4f\n", nm, x$beta[nm]))
  }
  invisible(x)
}
