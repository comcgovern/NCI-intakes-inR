#' INDIVINT: Individual Predicted Usual Intake (BLUP-based)
#'
#' R port of the NCI INDIVINT SAS macro. Computes individual-level predictions
#' of usual (habitual) dietary intake using the Best Linear Unbiased Predictor
#' (BLUP / empirical Bayes estimate) from the mixed model.
#'
#' The primary use case is **regression calibration**: using predicted usual
#' intake as a covariate in a health outcome model to correct for measurement
#' error in observed dietary recalls.
#'
#' @name indivint
NULL


#' Predict individual usual intake (INDIVINT equivalent)
#'
#' Uses empirical Bayes (BLUP) random-effect estimates from the fitted
#' `mixtran_fit` object to produce subject-level predictions of usual intake
#' on the original (back-transformed) scale.
#'
#' **Amount-only model:** `BLUP_i = X_ref_i'β̂ + û_i`, where `û_i` is the
#' estimated random intercept from `nlme::lme`. The prediction is evaluated at
#' "reference" covariate values: the sequence indicator (`seq_num`) is set to 0
#' (first recall), and the weekend indicator is set to 0, following NCI
#' convention. Other subject-level covariates retain their observed (averaged)
#' values.
#'
#' **Two-part model:** `prob_i = logistic(Z_ref_i'α̂ + v̂1_i)` and
#' `amt_i = T⁻¹(X_ref_i'β̂ + v̂2_i)`. Individual usual intake is
#' `prob_i × amt_i`. Subjects who never consumed the food have `v̂2_i`
#' imputed to 0 (population mean on the transformed scale).
#'
#' @param mixtran_obj A `mixtran_fit` object returned by [mixtran()].
#' @param zero_seq Logical. If `TRUE` (default), set the sequence number
#'   covariate to 0 when computing reference predictions (NCI convention for
#'   "usual day" predictions). Has no effect if `seq_num` was not in the model.
#' @param zero_weekend Logical. If `TRUE` (default), set the weekend indicator
#'   to 0 when computing reference predictions. Has no effect if no weekend
#'   variable was in the model.
#' @return A data frame with one row per subject containing:
#'   \describe{
#'     \item{subject}{Subject identifier}
#'     \item{predicted_usual}{Predicted usual intake on the **transformed**
#'       scale (Box-Cox transformed)}
#'     \item{predicted_usual_orig}{Predicted usual intake on the **original**
#'       scale (back-transformed via Box-Cox inverse)}
#'     \item{prob_usual}{(Two-part only) Predicted probability of consumption}
#'     \item{amt_usual_orig}{(Two-part only) Predicted amount given consumption,
#'       original scale}
#'   }
#' @references
#' Kipnis V, et al. (2009). Modeling data with excess zeros and measurement
#' error. *Biometrics* 65(4):1003-10.
#' @export
indivint <- function(mixtran_obj, zero_seq = TRUE, zero_weekend = TRUE) {

  stopifnot(inherits(mixtran_obj, "mixtran_fit"))

  model_type <- mixtran_obj$model_type
  lambda     <- mixtran_obj$lambda
  pred       <- mixtran_obj$predicted

  if (model_type == "amount") {
    .indivint_amount(mixtran_obj, pred, lambda, zero_seq, zero_weekend)
  } else {
    .indivint_twopart(mixtran_obj, pred, lambda, zero_seq, zero_weekend)
  }
}


# ---------------------------------------------------------------------------
# Internal: amount-only BLUP prediction
# ---------------------------------------------------------------------------

#' @keywords internal
.indivint_amount <- function(obj, pred, lambda, zero_seq, zero_weekend) {

  fit      <- obj$model_fit
  beta     <- obj$beta
  subjects <- unique(pred$subject)

  # predict(fit, level = 1) returns X_ij'beta + u_i for each observation.
  # We want the prediction for a "usual day" (seq_num = 0, weekend = 0).
  # Strategy: take level-1 prediction, subtract the unwanted covariate effects,
  # then average per subject.

  level1 <- stats::predict(fit, level = 1)   # length = nrow(pred)

  # Adjust for reference covariate values
  level1 <- .subtract_ref_effects(level1, pred, beta, zero_seq, zero_weekend)

  # Average over observations within each subject
  blup_t <- tapply(level1, pred$subject, mean)  # named by subject

  result <- data.frame(
    subject             = subjects,
    predicted_usual     = as.numeric(blup_t[as.character(subjects)]),
    stringsAsFactors    = FALSE
  )
  result$predicted_usual_orig <- pmax(
    boxcox_inverse(result$predicted_usual, lambda), 0
  )

  result
}


# ---------------------------------------------------------------------------
# Internal: two-part BLUP prediction
# ---------------------------------------------------------------------------

#' @keywords internal
.indivint_twopart <- function(obj, pred, lambda, zero_seq, zero_weekend) {

  prob_fit <- obj$prob_fit
  amt_fit  <- obj$amt_fit
  alpha    <- obj$alpha
  beta     <- obj$beta
  subjects <- unique(pred$subject)

  # ---- Probability sub-model ----
  # Population-level linear predictor (no RE) per observation
  prob_linpred_popn <- prob_popn_predict(prob_fit, newdata = pred)

  # Adjust to reference covariate values
  prob_linpred_popn <- .subtract_ref_effects(
    prob_linpred_popn, pred, alpha, zero_seq, zero_weekend
  )

  # Subject-level random effects from the prob model (named by subject ID)
  v1_re <- prob_ranef(prob_fit)  # named numeric vector

  # Per-subject prob BLUP: average population linpred + subject RE
  prob_blup <- vapply(as.character(subjects), function(s) {
    popn_mu <- mean(prob_linpred_popn[pred$subject == s], na.rm = TRUE)
    re_val  <- if (s %in% names(v1_re)) v1_re[s] else 0
    popn_mu + re_val
  }, numeric(1))

  prob_usual <- stats::plogis(prob_blup)

  # ---- Amount sub-model ----
  # Population-level predictions for consumed rows only
  pos_mask <- pred$consumed == 1
  pos_pred <- pred[pos_mask, ]

  amt_linpred_popn <- stats::predict(amt_fit, level = 0)  # consumed rows only
  amt_linpred_popn <- .subtract_ref_effects(
    amt_linpred_popn, pos_pred, beta, zero_seq, zero_weekend
  )

  # Subject-level RE from amount model (named by subject ID)
  v2_re_df <- nlme::ranef(amt_fit)[[1]]  # data frame: rows=subjects
  v2_re    <- setNames(as.numeric(v2_re_df[[1]]), rownames(v2_re_df))

  # Population-mean transformed amount (for subjects who never consumed)
  global_amt_mu <- mean(amt_linpred_popn, na.rm = TRUE)

  amt_blup_t <- vapply(as.character(subjects), function(s) {
    if (s %in% names(v2_re)) {
      consumed_rows <- pos_pred$subject == s
      popn_mu <- if (any(consumed_rows)) {
        mean(amt_linpred_popn[consumed_rows], na.rm = TRUE)
      } else {
        global_amt_mu
      }
      popn_mu + v2_re[s]
    } else {
      # Never-consumer: BLUP amount is the population mean (v2 shrinks to 0)
      global_amt_mu
    }
  }, numeric(1))

  amt_usual_orig <- pmax(boxcox_inverse(amt_blup_t, lambda), 0)

  # ---- Combine ----
  data.frame(
    subject             = subjects,
    predicted_usual     = prob_blup * amt_blup_t,  # on transformed scale (approx)
    predicted_usual_orig = prob_usual * amt_usual_orig,
    prob_usual          = prob_usual,
    amt_usual_orig      = amt_usual_orig,
    stringsAsFactors    = FALSE
  )
}


# ---------------------------------------------------------------------------
# Helper: subtract fixed effects for reference-covariate adjustment
# ---------------------------------------------------------------------------

#' Subtract fixed-effect contributions of seq_num and/or weekend from a
#' linear predictor vector, so the result represents a reference ("usual") day.
#'
#' @keywords internal
.subtract_ref_effects <- function(linpred, obs_df, coef_vec,
                                   zero_seq, zero_weekend) {

  # seq_num: always the first covariate added by prepare_mixtran_data
  if (zero_seq && "seq_num" %in% names(coef_vec) &&
      "seq_num" %in% names(obs_df)) {
    linpred <- linpred - coef_vec["seq_num"] * obs_df$seq_num
  }

  # weekend
  if (zero_weekend && "weekend" %in% names(coef_vec) &&
      "weekend" %in% names(obs_df)) {
    linpred <- linpred - coef_vec["weekend"] * obs_df$weekend
  }

  linpred
}


#' Print method for indivint results
#'
#' @param x Data frame returned by [indivint()].
#' @param n Number of rows to display. Default 10.
#' @param ... Ignored.
#' @export
print_indivint <- function(x, n = 10, ...) {
  cat("Individual Usual Intake Predictions (INDIVINT)\n")
  cat(sprintf("  N subjects: %d\n", nrow(x)))
  cat(sprintf("  Mean predicted usual intake: %.2f\n",
              mean(x$predicted_usual_orig, na.rm = TRUE)))
  cat(sprintf("  Range: [%.2f, %.2f]\n",
              min(x$predicted_usual_orig, na.rm = TRUE),
              max(x$predicted_usual_orig, na.rm = TRUE)))
  cat("\n")
  print(utils::head(x, n))
  invisible(x)
}
