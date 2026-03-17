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
#' Uses a two-step approach:
#' 1. Fit uncorrelated model to get starting values
#' 2. Maximize joint log-likelihood with Gauss-Hermite quadrature
#'
#' @keywords internal
fit_twopart_corr <- function(prep, lambda, verbose) {

  if (verbose) message("Fitting two-part correlated model...")

  # Step 1: Get starting values from uncorrelated model
  if (verbose) message("  Step 1: Getting starting values from uncorrelated model...")
  uncorr <- fit_twopart_uncorr(prep, lambda, verbose = FALSE)

  if (!uncorr$converged) {
    warning("Uncorrelated model failed. Cannot fit correlated model.")
    uncorr$model_type <- "corr"
    uncorr$rho <- NA
    return(uncorr)
  }

  # Step 2: Joint optimization
  if (verbose) message("  Step 2: Joint optimization with Gauss-Hermite quadrature...")

  work <- prep$data
  cov_names <- prep$cov_names

  # Build design matrices
  if (length(cov_names) > 0) {
    X_full <- stats::model.matrix(
      stats::as.formula(paste("~", paste(cov_names, collapse = " + "))),
      data = work
    )
  } else {
    X_full <- matrix(1, nrow = nrow(work), ncol = 1)
  }

  # Gauss-Hermite quadrature nodes and weights (for bivariate integration)
  gh <- gauss_hermite_2d(n_nodes = 10)

  # Pack starting values
  n_beta <- length(uncorr$beta)
  n_alpha <- length(uncorr$alpha)

  start <- c(
    uncorr$alpha,                              # alpha (probability fixed effects)
    uncorr$beta,                               # beta (amount fixed effects)
    log(uncorr$sigma2_v1),                     # log(var of prob random effect)
    log(uncorr$sigma2_v2),                     # log(var of amount random effect)
    log(uncorr$sigma2_e),                      # log(within-person variance)
    0                                           # atanh(rho) = 0 (start uncorrelated)
  )

  # Organize data by subject for the likelihood
  subjects <- unique(work$subject)
  subj_data <- lapply(subjects, function(s) {
    idx <- which(work$subject == s)
    list(
      consumed = work$consumed[idx],
      intake = work$intake[idx],
      X = X_full[idx, , drop = FALSE],
      weight = work$weight[idx][1]
    )
  })

  # Joint negative log-likelihood
  neg_loglik <- function(par) {
    alpha <- par[1:n_alpha]
    beta <- par[(n_alpha + 1):(n_alpha + n_beta)]
    log_s2_v1 <- par[n_alpha + n_beta + 1]
    log_s2_v2 <- par[n_alpha + n_beta + 2]
    log_s2_e <- par[n_alpha + n_beta + 3]
    atanh_rho <- par[n_alpha + n_beta + 4]

    s2_v1 <- exp(log_s2_v1)
    s2_v2 <- exp(log_s2_v2)
    s2_e <- exp(log_s2_e)
    rho <- tanh(atanh_rho)

    sd_v1 <- sqrt(s2_v1)
    sd_v2 <- sqrt(s2_v2)

    total_ll <- 0

    for (sd in subj_data) {
      # For this subject, integrate over bivariate random effects (v1, v2)
      # using Gauss-Hermite quadrature
      subj_ll_vals <- vapply(seq_len(nrow(gh$nodes)), function(q) {
        # Transform standard normal nodes to correlated random effects
        z1 <- gh$nodes[q, 1]
        z2 <- gh$nodes[q, 2]
        v1 <- sd_v1 * z1
        v2 <- sd_v2 * (rho * z1 + sqrt(1 - rho^2) * z2)

        ll_q <- 0
        for (j in seq_along(sd$consumed)) {
          # Probability component
          eta_p <- sum(sd$X[j, ] * alpha) + v1
          p_consume <- 1 / (1 + exp(-eta_p))
          p_consume <- pmin(pmax(p_consume, 1e-10), 1 - 1e-10)

          if (sd$consumed[j] == 0) {
            ll_q <- ll_q + log(1 - p_consume)
          } else {
            ll_q <- ll_q + log(p_consume)
            # Amount component
            t_y <- boxcox_transform(sd$intake[j], lambda)
            mu_a <- sum(sd$X[j, ] * beta) + v2
            ll_q <- ll_q - 0.5 * log(2 * pi * s2_e) -
              0.5 * (t_y - mu_a)^2 / s2_e
            # Jacobian of Box-Cox: (lambda - 1) * log(y)
            if (abs(lambda) > 1e-10) {
              ll_q <- ll_q + (lambda - 1) * log(sd$intake[j])
            } else {
              ll_q <- ll_q - log(sd$intake[j])
            }
          }
        }

        # Return on regular scale (not log) for quadrature summation
        exp(ll_q) * gh$weights[q]
      }, numeric(1))

      marginal_lik <- sum(subj_ll_vals)
      if (marginal_lik <= 0) marginal_lik <- 1e-300
      total_ll <- total_ll + sd$weight * log(marginal_lik)
    }

    return(-total_ll)
  }

  # Optimize
  opt_result <- tryCatch({
    stats::optim(
      par = start,
      fn = neg_loglik,
      method = "BFGS",
      control = list(maxit = 500, reltol = 1e-8, trace = if (verbose) 1 else 0)
    )
  }, error = function(e) {
    if (verbose) message("    BFGS failed: ", e$message, ". Trying Nelder-Mead...")
    stats::optim(
      par = start,
      fn = neg_loglik,
      method = "Nelder-Mead",
      control = list(maxit = 5000, reltol = 1e-8)
    )
  })

  # Unpack results
  par <- opt_result$par
  alpha_hat <- par[1:n_alpha]
  names(alpha_hat) <- names(uncorr$alpha)
  beta_hat <- par[(n_alpha + 1):(n_alpha + n_beta)]
  names(beta_hat) <- names(uncorr$beta)
  sigma2_v1_hat <- exp(par[n_alpha + n_beta + 1])
  sigma2_v2_hat <- exp(par[n_alpha + n_beta + 2])
  sigma2_e_hat <- exp(par[n_alpha + n_beta + 3])
  rho_hat <- tanh(par[n_alpha + n_beta + 4])

  converged <- opt_result$convergence == 0

  if (verbose) {
    message(sprintf("  Convergence: %s (code %d)", converged, opt_result$convergence))
    message(sprintf("  sigma2_v1 (prob):   %.4f", sigma2_v1_hat))
    message(sprintf("  sigma2_v2 (amount): %.4f", sigma2_v2_hat))
    message(sprintf("  sigma2_e (within):  %.4f", sigma2_e_hat))
    message(sprintf("  rho:                %.4f", rho_hat))
  }

  # Build predicted values (linear predictors)
  pred_df <- work
  pred_df$prob_linpred <- as.numeric(X_full %*% alpha_hat)
  pred_df$amt_linpred <- as.numeric(X_full %*% beta_hat)

  list(
    alpha = alpha_hat,
    beta = beta_hat,
    sigma2_v1 = sigma2_v1_hat,
    sigma2_v2 = sigma2_v2_hat,
    sigma2_e = sigma2_e_hat,
    rho = rho_hat,
    optim_result = opt_result,
    predicted = pred_df,
    converged = converged
  )
}


#' Generate 2D Gauss-Hermite quadrature nodes and weights
#'
#' Produces the tensor product of 1D GH quadrature for bivariate integration
#' over standard normal random effects.
#'
#' @param n_nodes Number of nodes per dimension
#' @return List with nodes (matrix) and weights (vector)
#' @keywords internal
gauss_hermite_2d <- function(n_nodes = 10) {
  # 1D Gauss-Hermite nodes and weights
  gh1 <- gauss_hermite_1d(n_nodes)

  # Tensor product for 2D
  grid <- expand.grid(i = seq_len(n_nodes), j = seq_len(n_nodes))
  nodes <- cbind(gh1$nodes[grid$i], gh1$nodes[grid$j])
  weights <- gh1$weights[grid$i] * gh1$weights[grid$j]

  list(nodes = nodes, weights = weights)
}


#' Generate 1D Gauss-Hermite quadrature nodes and weights
#'
#' For integration of f(x)*exp(-x^2), converted to standard normal
#' weights (i.e., f(x)*dnorm(x)) via rescaling.
#'
#' @param n Number of nodes
#' @return List with nodes and weights for standard normal integration
#' @keywords internal
gauss_hermite_1d <- function(n) {
  # Eigenvalue method for Gauss-Hermite quadrature
  # The nodes are eigenvalues of the companion matrix
  i <- seq_len(n - 1)
  b <- sqrt(i / 2)
  cm <- diag(0, n)
  for (k in seq_along(b)) {
    cm[k, k + 1] <- b[k]
    cm[k + 1, k] <- b[k]
  }
  eig <- eigen(cm, symmetric = TRUE)
  nodes <- eig$values  # These are for exp(-x^2) weighting

  # Weights for exp(-x^2) weighting
  weights <- (eig$vectors[1, ])^2 * sqrt(pi)

  # Convert to standard normal: x_new = x * sqrt(2), w_new = w / sqrt(pi)
  nodes <- nodes * sqrt(2)
  weights <- weights / sqrt(pi)

  # Sort by nodes
  ord <- order(nodes)
  list(nodes = nodes[ord], weights = weights[ord])
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
