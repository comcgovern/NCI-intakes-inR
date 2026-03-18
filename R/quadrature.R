#' Gauss-Hermite Quadrature Utilities
#'
#' Nodes and weights for Gauss-Hermite quadrature, used by the GHQ engine
#' of the correlated two-part model. Supports the physicist's convention
#' (integrating against exp(-x²)) and the probabilist's convention
#' (integrating against the standard normal density).
#'
#' Hardcoded high-precision node/weight tables for n = 3, 5, 7, 9 are
#' provided. For other values of n the nodes and weights are computed via
#' the Golub-Welsch eigenvalue algorithm.
#'
#' @name quadrature
NULL


# ---------------------------------------------------------------------------
# Hardcoded tables (physicist convention: integral of f(x) exp(-x²) dx)
# ---------------------------------------------------------------------------

.GH_TABLES <- list(
  `3` = list(
    nodes   = c(-1.2247448713916, 0.0, 1.2247448713916),
    weights = c(0.2954089751509, 1.1816359006037, 0.2954089751509)
  ),
  `5` = list(
    nodes   = c(-2.0201828704561, -0.9585724646138, 0.0,
                 0.9585724646138,  2.0201828704561),
    weights = c(0.0199532420590, 0.3942467176458, 0.9453087204829,
                0.3942467176458, 0.0199532420590)
  ),
  `7` = list(
    nodes   = c(-2.6519613568352, -1.6735516287675, -0.8162877718535, 0.0,
                 0.8162877718535,  1.6735516287675,  2.6519613568352),
    weights = c(0.0009717812450, 0.0545155828191, 0.4256072526101,
                0.8102646175568, 0.4256072526101, 0.0545155828191,
                0.0009717812450)
  ),
  `9` = list(
    nodes   = c(-3.1909932017815, -2.2665805084533, -1.4685532892167,
                -0.7235510187528,  0.0,
                 0.7235510187528,  1.4685532892167,  2.2665805084533,
                 3.1909932017815),
    weights = c(0.0000724044820, 0.0065745384742, 0.1060263195642,
                0.6587818792132, 0.9450677499676,
                0.6587818792132, 0.1060263195642, 0.0065745384742,
                0.0000724044820)
  )
)


#' Gauss-Hermite quadrature nodes and weights
#'
#' Returns nodes x_q and weights w_q for the physicist's GH rule:
#'   `integral f(x) exp(-x^2) dx  ≈  sum_q  w_q  f(x_q)`
#'
#' To integrate against N(0, 1) use [gh_nodes_normal()].
#'
#' @param n Number of quadrature points (3, 5, 7, or 9 use pre-computed
#'   tables; other positive integers are computed via eigenvalue method).
#' @return List with components `nodes` and `weights`.
#' @export
gh_nodes_weights <- function(n) {
  key <- as.character(n)
  if (!is.null(.GH_TABLES[[key]])) {
    return(.GH_TABLES[[key]])
  }
  # Golub-Welsch: tridiagonal symmetric eigenproblem
  # Jacobi matrix for Hermite polynomials: off-diagonal = sqrt(i/2), i=1,...,n-1
  off <- sqrt(seq_len(n - 1) / 2)
  J   <- diag(0, n)
  J[cbind(seq_len(n - 1), seq_len(n - 1) + 1)] <- off
  J[cbind(seq_len(n - 1) + 1, seq_len(n - 1))] <- off
  ev  <- eigen(J, symmetric = TRUE)
  # Weights = (first component of eigenvector)^2 * sqrt(pi)
  w   <- (ev$vectors[1, ])^2 * sqrt(pi)
  x   <- ev$values
  ord <- order(x)
  list(nodes = x[ord], weights = w[ord])
}


#' Gauss-Hermite nodes / weights for N(0, 1) integration
#'
#' Transforms the physicist GH rule to the probabilist convention:
#'   `integral f(x) phi(x) dx  ≈  sum_q  w_q^*  f(x_q^*)`
#'
#' Transformation:  x_q^* = sqrt(2) * x_q,  w_q^* = w_q / sqrt(pi)
#'
#' The weights w_q^* sum to 1.
#'
#' @param n Number of quadrature points.
#' @return List with components `nodes` (on normal scale) and `weights`
#'   (sum to 1).
#' @export
gh_nodes_normal <- function(n) {
  gh <- gh_nodes_weights(n)
  list(
    nodes   = sqrt(2) * gh$nodes,
    weights = gh$weights / sqrt(pi)
  )
}


#' Bivariate Gauss-Hermite quadrature nodes
#'
#' Generates the tensor-product bivariate quadrature nodes for integrating
#' against BVN(0, Sigma), where Sigma has variances (sigma2_v1, sigma2_v2)
#' and off-diagonal rho * sqrt(sigma2_v1 * sigma2_v2).
#'
#' The Cholesky decomposition of Sigma is used to transform the standard
#' bivariate normal nodes:
#'   v1 = sigma_v1 * z1
#'   v2 = sigma_v2 * (rho * z1 + sqrt(1 - rho^2) * z2)
#'
#' @param n Number of points per dimension (total n^2 bivariate nodes).
#' @param rho Correlation parameter in (-1, 1).
#' @param sigma_v1 Standard deviation of the first random effect.
#' @param sigma_v2 Standard deviation of the second random effect.
#' @return Data frame with columns v1, v2, log_weight (log of the bivariate
#'   quadrature weight, i.e. log(w1 * w2)).
#' @export
gh_nodes_bivariate <- function(n, rho, sigma_v1, sigma_v2) {
  gh   <- gh_nodes_normal(n)  # nodes on N(0,1) scale, weights sum to 1
  z    <- gh$nodes
  w    <- gh$weights
  rho  <- min(max(rho, -0.9999), 0.9999)
  sqrt_1r2 <- sqrt(1 - rho^2)

  # Tensor product: all pairs (z1, z2)
  n_total <- n^2
  idx1 <- rep(seq_len(n), each = n)
  idx2 <- rep(seq_len(n), times = n)

  z1  <- z[idx1]
  z2  <- z[idx2]
  lw  <- log(w[idx1]) + log(w[idx2])

  v1 <- sigma_v1 * z1
  v2 <- sigma_v2 * (rho * z1 + sqrt_1r2 * z2)

  data.frame(v1 = v1, v2 = v2, log_weight = lw, stringsAsFactors = FALSE)
}


#' Adaptive Gauss-Hermite quadrature nodes for a single subject
#'
#' Places quadrature nodes at the mode and scales them by the posterior
#' curvature (Laplace-like centering), improving accuracy for subjects
#' with many observations where the posterior is concentrated.
#'
#' For each subject i, the mode (v1*, v2*) of the integrand is found by
#' Newton steps, then the standard nodes are shifted and scaled:
#'   z_q_adapted = mode + scale * z_q_standard
#'
#' This mirrors the adaptive GH used in SAS PROC NLMIXED.
#'
#' @param n Number of quadrature points per dimension.
#' @param eta_i Numeric vector of probability linear predictors for subject i's recalls.
#' @param mu_i Numeric vector of amount linear predictors (positive recalls only).
#' @param consumed_i Binary vector (0/1) for subject i's recalls.
#' @param t_y_i Numeric vector of Box-Cox transformed intake for positive recalls.
#' @param sigma_v1,sigma_v2,sigma_e,rho Variance component parameters.
#' @param lambda Box-Cox lambda.
#' @param max_iter Maximum Newton iterations for mode finding (default 20).
#' @param tol Convergence tolerance (default 1e-6).
#' @return Data frame with v1, v2, log_weight columns (adapted nodes).
#' @export
gh_nodes_adaptive <- function(n, eta_i, mu_i, consumed_i, t_y_i,
                               sigma_v1, sigma_v2, sigma_e, rho, lambda,
                               max_iter = 20, tol = 1e-6) {

  rho      <- min(max(rho, -0.9999), 0.9999)
  sqrt_1r2 <- sqrt(1 - rho^2)

  # Precision matrix of the bivariate normal prior
  inv_sigma <- matrix(
    c(1/sigma_v1^2, -rho/(sigma_v1*sigma_v2*sqrt_1r2^2 * sigma_v1 * sigma_v2),
      -rho/(sigma_v1*sigma_v2*sqrt_1r2^2 * sigma_v1 * sigma_v2), 1/sigma_v2^2),
    2, 2
  )
  # Correct precision matrix for BVN(0, Sigma)
  sigma_mat <- matrix(
    c(sigma_v1^2, rho*sigma_v1*sigma_v2,
      rho*sigma_v1*sigma_v2, sigma_v2^2),
    2, 2
  )
  inv_sigma <- solve(sigma_mat)

  # Log-posterior (unnormalised) and its gradient/Hessian for Newton steps
  log_post <- function(v) {
    v1 <- v[1]; v2 <- v[2]
    # Prior
    lp_prior <- -0.5 * as.numeric(t(v) %*% inv_sigma %*% v)
    # Likelihood — probability part
    eta_q  <- eta_i + v1
    lp_bin <- sum(ifelse(consumed_i == 1,
                         stats::plogis(eta_q, log.p = TRUE),
                         stats::plogis(eta_q, log.p = TRUE, lower.tail = FALSE)))
    # Likelihood — amount part
    lp_amt <- if (length(t_y_i) > 0) {
      sum(stats::dnorm(t_y_i, mean = mu_i + v2, sd = sigma_e, log = TRUE))
    } else 0
    lp_prior + lp_bin + lp_amt
  }

  # Find mode via Nelder-Mead (robust, no gradient needed for the adapt step)
  opt <- tryCatch(
    stats::optim(c(0, 0), function(v) -log_post(v),
                 method  = "Nelder-Mead",
                 control = list(maxit = max_iter * 50, reltol = tol)),
    error = function(e) NULL
  )
  mode_v <- if (!is.null(opt) && opt$convergence == 0) opt$par else c(0, 0)

  # Approximate Hessian at mode via finite differences → scale matrix
  eps <- 1e-4
  H <- matrix(0, 2, 2)
  for (j in 1:2) {
    e_j <- rep(0, 2); e_j[j] <- eps
    H[j, j] <- (log_post(mode_v + e_j) - 2*log_post(mode_v) +
                   log_post(mode_v - e_j)) / eps^2
  }
  e12 <- c(eps, eps)
  H[1,2] <- H[2,1] <- (log_post(mode_v + e12) - log_post(mode_v + c(eps,-eps)) -
                          log_post(mode_v + c(-eps,eps)) + log_post(mode_v - e12)) /
                        (4*eps^2)

  # Scale matrix: Cholesky of -H^{-1} (covariance of Gaussian approx)
  neg_inv_H <- tryCatch(-solve(H), error = function(e) diag(c(sigma_v1^2, sigma_v2^2)))
  # Ensure positive definite
  eig <- eigen(neg_inv_H, symmetric = TRUE)
  eig$values <- pmax(eig$values, 1e-8)
  scale_mat  <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  log_det_scale <- sum(log(pmax(eig$values, 1e-8))) / 2

  # Standard bivariate GH nodes (unit normal)
  gh  <- gh_nodes_normal(n)
  z   <- gh$nodes
  w   <- gh$weights
  idx1 <- rep(seq_len(n), each = n)
  idx2 <- rep(seq_len(n), times = n)
  z_std <- cbind(z[idx1], z[idx2])

  # Shift and scale nodes
  v_nodes <- t(scale_mat %*% t(z_std)) + matrix(mode_v, nrow = nrow(z_std), ncol = 2, byrow = TRUE)

  # Adjusted log-weights: log(w1*w2) + log|scale_mat|
  lw_std  <- log(w[idx1]) + log(w[idx2])
  lw_adj  <- lw_std + log_det_scale

  data.frame(v1 = v_nodes[, 1], v2 = v_nodes[, 2],
             log_weight = lw_adj, stringsAsFactors = FALSE)
}
