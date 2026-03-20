#' Box-Cox Transformation Utilities
#'
#' The NCI method applies a Box-Cox transformation to consumption-day amounts
#' to approximate normality before fitting the mixed model. These functions
#' handle the transformation, back-transformation, and optimal lambda selection.

#' Apply Box-Cox transformation
#'
#' @param y Numeric vector of positive values
#' @param lambda Box-Cox parameter. lambda=0 gives log transform.
#' @return Transformed values
#' @export
boxcox_transform <- function(y, lambda) {
  stopifnot(all(y[!is.na(y)] > 0))
  if (abs(lambda) < 1e-10) {
    return(log(y))
  } else {
    return((y^lambda - 1) / lambda)
  }
}

#' Inverse Box-Cox transformation
#'
#' @param z Numeric vector of transformed values
#' @param lambda Box-Cox parameter
#' @return Back-transformed values on original scale
#' @export
boxcox_inverse <- function(z, lambda) {
  if (abs(lambda) < 1e-10) {
    return(exp(z))
  } else {
    # Clamp to avoid negative values inside power
    inner <- lambda * z + 1
    inner <- pmax(inner, 1e-10)
    return(inner^(1 / lambda))
  }
}

#' Find optimal Box-Cox lambda via profile likelihood
#'
#' Searches over a grid of lambda values to find the transformation that
#' best normalizes the positive consumption-day amounts. Mirrors the
#' approach used in SAS MIXTRAN, which evaluates lambda on a grid from
#' 0.01 to 1.0 (or uses log if lambda < 0.15).
#'
#' @param y Numeric vector of positive consumption-day amounts
#' @param weights Optional survey weights
#' @param lambda_grid Grid of lambda values to search
#' @param covariates Optional model matrix of covariates (for residual normality)
#' @return List with optimal lambda, log-likelihoods, and diagnostics
#' @export
find_optimal_lambda <- function(y,
                                 weights = NULL,
                                 lambda_grid = seq(0.01, 1.0, by = 0.01),
                                 covariates = NULL) {
  y_pos <- y[!is.na(y) & y > 0]
  n <- length(y_pos)

  if (n == 0L) {
    stop("No positive values in 'y'; cannot estimate Box-Cox lambda. ",
         "Ensure the intake variable has at least some positive observations, ",
         "or supply a fixed lambda.")
  }

  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }
  w_pos <- weights[!is.na(y) & y > 0]

  # Profile log-likelihood for Box-Cox (includes Jacobian)
  log_liks <- vapply(lambda_grid, function(lam) {
    z <- boxcox_transform(y_pos, lam)

    if (!is.null(covariates)) {
      # Residuals after regressing on covariates
      cov_pos <- covariates[!is.na(y) & y > 0, , drop = FALSE]
      fit <- tryCatch(
        stats::lm.wfit(x = cbind(1, cov_pos), y = z, w = w_pos),
        error = function(e) NULL
      )
      if (is.null(fit)) return(-Inf)
      resid <- fit$residuals
    } else {
      resid <- z - stats::weighted.mean(z, w_pos)
    }

    # Weighted variance of residuals
    sigma2 <- sum(w_pos * resid^2) / sum(w_pos)
    if (sigma2 <= 0) return(-Inf)

    # Profile log-likelihood
    ll <- -0.5 * n * log(sigma2) + (lam - 1) * sum(w_pos * log(y_pos)) / sum(w_pos) * n
    return(ll)
  }, numeric(1))

  best_idx <- which.max(log_liks)
  best_lambda <- lambda_grid[best_idx]

  # NCI convention: if lambda < 0.15, use log transform (lambda = 0)
  # because the method may not perform well with very small lambda
  use_log <- best_lambda < 0.15

  list(
    lambda = if (use_log) 0 else best_lambda,
    lambda_raw = best_lambda,
    used_log = use_log,
    log_likelihoods = data.frame(lambda = lambda_grid, loglik = log_liks),
    n_positive = n
  )
}
