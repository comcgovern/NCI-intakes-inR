#' INDIVINT: Individual Predicted Usual Intake
#'
#' R port of the NCI INDIVINT SAS macro (v2.1). Computes, for each subject,
#' the conditional expectation of usual intake given that subject's observed
#' 24-hour recalls (Kipnis et al., 2009), for use in regression calibration.
#'
#' @name indivint
NULL


#' Predict individual usual intake (INDIVINT equivalent)
#'
#' For each subject INDIVINT computes
#' \deqn{E[T_i \mid \mathrm{recalls}_i] = \frac{\int T_i(u) f(\mathrm{recalls}_i \mid u)\,\phi(u)\,du}
#'       {\int f(\mathrm{recalls}_i \mid u)\,\phi(u)\,du}}
#' where `u` are the subject's random effects and
#' `T_i(u) = p_i(u1) * g*(x2b2_i + u2)` is usual intake: `p_i` is the
#' consumption probability (1 for amount-only models) and `g*` is the
#' back-transformation integrated over the within-person error (exact when
#' lambda = 0, otherwise INDIVINT's 9-point rule, each point bounded below by
#' half the minimum amount). The INDIVINT macro evaluates both integrals by
#' adaptive Gauss-Hermite quadrature. Here, conditional on `u1` the posterior
#' of `u2` is normal, so the `u2` integral uses Gauss-Hermite nodes on that
#' normal and the `u1` integral uses adaptive Gauss-Hermite quadrature centred
#' at the posterior mode.
#'
#' The recall likelihood uses each recall's own linear predictor (including
#' sequence and weekend effects). The usual-intake function `T_i` uses
#' sequence = 0 and, when the model has a weekend indicator, averages weekday
#' and weekend usual intake with weights `1 - wkend_prop` and `wkend_prop`, as
#' in [distrib()]. (The SAS macro takes a single user-supplied `xbeta` per
#' subject; for models without sequence or weekend effects the two coincide.)
#'
#' @param mixtran_obj A `mixtran_fit` object returned by [mixtran()].
#' @param lamt Optional Box-Cox parameter (non-zero) for predicting usual
#'   intake on a transformed scale (INDIVINT `boxcox_t_lamt=y, lamt=`). When
#'   given, `indusint` is `E[(T^lamt - 1) / lamt | recalls]`.
#' @param min_amt Minimum amount; half of it bounds the back-transformed
#'   amount (INDIVINT `min_amt`). Defaults to the smallest positive recall.
#' @param wkend_prop Weekend weight used when the model includes a weekend
#'   indicator. Default 3/7.
#' @param n_quad Quadrature points per dimension. Default 15.
#' @param zero_seq,zero_weekend Deprecated; ignored.
#' @return A data frame with one row per subject:
#'   \describe{
#'     \item{subject}{Subject identifier}
#'     \item{indusint}{INDIVINT's predicted value: `E[T | recalls]`, or
#'       `E[BoxCox(T, lamt) | recalls]` when `lamt` is given}
#'     \item{predicted_usual}{Same as `indusint` (kept for compatibility)}
#'     \item{predicted_usual_orig}{`E[T | recalls]` on the original scale}
#'     \item{prob_usual}{(Two-part only) `E[p | recalls]`}
#'     \item{amt_usual_orig}{(Two-part only) `E[g* | recalls]`}
#'   }
#' @references
#' Kipnis V, et al. (2009). Modeling data with excess zeros and measurement
#' error. *Biometrics* 65(4):1003-10.
#' @export
indivint <- function(mixtran_obj, lamt = NULL, min_amt = NULL,
                     wkend_prop = 3 / 7, n_quad = 15L,
                     zero_seq = NULL, zero_weekend = NULL) {

  stopifnot(inherits(mixtran_obj, "mixtran_fit"))
  if (!is.null(zero_seq) || !is.null(zero_weekend)) {
    warning("'zero_seq' and 'zero_weekend' are deprecated and ignored.", call. = FALSE)
  }
  if (!is.null(lamt) && (!is.numeric(lamt) || lamt == 0)) {
    stop("'lamt' must be a non-zero number (INDIVINT does not allow 0).")
  }
  fit <- mixtran_obj
  if (is.null(min_amt)) min_amt <- fit$min_pos_raw %||% fit$min_amt
  floor_amt <- 0.5 * min_amt

  work <- fit$data
  q <- .unpack(fit$theta, list(lambda_fixed = fit$lambda))
  lambda <- fit$lambda
  pred <- fit$predicted
  N    <- nrow(pred)
  id   <- match(work$subject, pred$subject)
  pos  <- work$consumed == 1
  has_wk <- "x2b2_1" %in% names(pred)
  wd <- 1 - wkend_prop

  # per-recall linear predictors and transformed amounts
  X2 <- cbind(1, as.matrix(work[pos, fit$amt_covs_used, drop = FALSE]))
  r  <- boxcox_transform(work$intake[pos], lambda) - as.numeric(X2 %*% fit$beta)
  S  <- .amt_sums(r, log(work$intake[pos]), lambda, id[pos],
                  sort(unique(id[pos])), N)

  gh <- gh_nodes_normal(n_quad)
  se2 <- q$se2; su2 <- q$su2

  # g*(mu): expected amount over within-person error (INDIVINT backtransform)
  gstar <- function(mu) .backtransform_indivint(mu, lambda, se2, floor_amt)
  tr <- function(t) if (is.null(lamt)) t else (t^lamt - 1) / lamt

  if (fit$model_type == "amount") {
    # posterior of u2 is N(m, v)
    m <- su2 * S$R / (se2 + S$k * su2)
    v <- su2 * se2 / (se2 + S$k * su2)
    U2 <- m + sqrt(v) %o% gh$nodes                          # N x Q
    A0 <- gstar(pred$x2b2 + U2)
    A  <- if (has_wk) wd * A0 + wkend_prop * gstar(pred$x2b2_1 + U2) else A0
    wq <- matrix(gh$weights, N, n_quad, byrow = TRUE)
    t_mean <- rowSums(wq * A)
    ind    <- rowSums(wq * tr(A))
    return(data.frame(subject = pred$subject, indusint = ind,
                      predicted_usual = ind, predicted_usual_orig = t_mean,
                      stringsAsFactors = FALSE))
  }

  # ---- two-part model ----
  X1  <- cbind(1, as.matrix(work[, fit$prob_covs_used, drop = FALSE]))
  eta <- as.numeric(X1 %*% fit$alpha)
  yn  <- work$consumed
  sig1 <- q$sig1; rho <- q$rho
  cc  <- rho * sqrt(su2) / sig1
  s2  <- su2 * (1 - rho^2)
  amt <- S

  # marginal posterior of u1 (u2 integrated out): adaptive GH around the mode
  md <- .prob_agq(eta, 2 * yn - 1, yn, id, N, sig1, amt = amt, c_shift = cc,
                  se2 = se2, s2 = s2, mode_only = TRUE)
  shat <- 1 / sqrt(-md$curv)
  gh1  <- .gh_std(n_quad)
  U1   <- md$mode + sqrt(2) * outer(shat, gh1$x)            # N x Q1
  lb   <- .by_subj(stats::plogis((2 * yn - 1) * (eta + U1[id, , drop = FALSE]),
                                 log.p = TRUE), id)
  G    <- lb + .amt_logf(amt, cc * U1, se2, s2) - 0.5 * U1^2 / sig1^2 +
          matrix(gh1$lw + gh1$x^2, N, n_quad, byrow = TRUE)
  W1   <- exp(G - apply(G, 1, max)); W1 <- W1 / rowSums(W1)  # posterior weights

  # u2 | u1, recalls ~ N(m(u1), v)
  den <- se2 + amt$k * s2
  v   <- s2 * se2 / den
  E_t <- E_ind <- E_p <- E_a <- numeric(N)
  for (a in seq_len(n_quad)) {
    u1 <- U1[, a]
    m  <- cc * u1 + s2 * (amt$R - amt$k * cc * u1) / den
    U2 <- m + sqrt(v) %o% gh$nodes                          # N x Q2
    p0 <- stats::plogis(pred$x1b1 + u1); A0 <- gstar(pred$x2b2 + U2)
    if (has_wk) {
      p1 <- stats::plogis(pred$x1b1_1 + u1); A1 <- gstar(pred$x2b2_1 + U2)
      Tm <- wd * p0 * A0 + wkend_prop * p1 * A1
      pm <- wd * p0 + wkend_prop * p1
      Am <- wd * A0 + wkend_prop * A1
    } else {
      Tm <- p0 * A0; pm <- p0; Am <- A0
    }
    wq <- matrix(gh$weights, N, n_quad, byrow = TRUE)
    E_t   <- E_t   + W1[, a] * rowSums(wq * Tm)
    E_ind <- E_ind + W1[, a] * rowSums(wq * tr(Tm))
    E_p   <- E_p   + W1[, a] * pm
    E_a   <- E_a   + W1[, a] * rowSums(wq * Am)
  }

  data.frame(subject = pred$subject, indusint = E_ind, predicted_usual = E_ind,
             predicted_usual_orig = E_t, prob_usual = E_p, amt_usual_orig = E_a,
             stringsAsFactors = FALSE)
}


# INDIVINT backtransform(): 9-point rule with normalised weights; every point
# is bounded below by half the minimum amount (gbc_inverse); exact
# exp(mu + var_e/2), then bounded, when lambda = 0.
.backtransform_indivint <- function(mu, lambda, var_e, floor_amt) {
  if (abs(lambda) < 1e-10) return(pmax(exp(mu + var_e / 2), floor_amt))
  x9 <- c(-2.1, -1.3, -0.8, -0.5, 0, 0.5, 0.8, 1.3, 2.1)
  w9 <- c(0.063345, 0.080255, 0.070458, 0.159698, 0.252489,
          0.159698, 0.070458, 0.080255, 0.063345)
  w9 <- w9 / sum(w9)
  out <- 0
  for (j in seq_along(x9)) {
    z <- mu + x9[j] * sqrt(var_e)
    out <- out + w9[j] * pmax(pmax(1 + lambda * z, 0)^(1 / lambda), floor_amt)
  }
  out
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
  cat(sprintf("  Mean predicted usual intake: %.4f\n",
              mean(x$predicted_usual_orig, na.rm = TRUE)))
  cat(sprintf("  Range: [%.4f, %.4f]\n",
              min(x$predicted_usual_orig, na.rm = TRUE),
              max(x$predicted_usual_orig, na.rm = TRUE)))
  cat("\n")
  print(utils::head(x, n))
  invisible(x)
}
