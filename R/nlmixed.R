#' Weighted maximum-likelihood engine mirroring the MIXTRAN PROC NLMIXED fits
#'
#' MIXTRAN v2.1 fits every model with PROC NLMIXED, using the survey weight in
#' the `REPLICATE` statement (each subject's marginal log-likelihood is
#' multiplied by its weight) and treating the Box-Cox parameter `A_LAMBDA` as a
#' free parameter (`bounds A_LAMBDA > 0.01`). This file implements the same
#' likelihoods:
#'
#' * amount part: `T(y; lambda) = X2 beta + u2 + e`, including the Jacobian
#'   `(lambda - 1) log(y)`. Given the shift of `u2`, the subject's amount
#'   likelihood is a closed-form compound-symmetric normal density.
#' * probability part: `logit P(y > 0) = X1 alpha + u1`.
#' * correlated model: `(u1, u2) ~ BVN`. Conditional on `u1`, `u2` is normal,
#'   so `u2` is integrated out exactly and only the one-dimensional integral
#'   over `u1` is computed numerically, by adaptive Gauss-Hermite quadrature
#'   centred at each subject's mode (as NLMIXED does by default).
#'
#' Parameterisation follows NLMIXED: log standard deviations
#' (`A_LOGSDE`, `A_LOGSDU2`, `P_LOGSDU1`) and Fisher's z for rho (`Z_U`).
#'
#' @name nlmixed_engine
#' @keywords internal
NULL


# Sum a row-level vector (or matrix) within subjects 1..N. Rows must be sorted
# by subject and every subject must have at least one row.
.by_subj <- function(x, id) {
  out <- rowsum(x, id, reorder = TRUE)
  if (is.matrix(x)) out else out[, 1]
}

# Per-subject sums of the amount-part residuals, with zero for subjects that
# have no positive recall.
.amt_sums <- function(r, logy, lambda, id_pos, uid_pos, N) {
  R <- Q <- J <- numeric(N)
  k <- integer(N)
  if (length(r) > 0L) {
    R[uid_pos] <- rowsum(r, id_pos, reorder = TRUE)[, 1]
    Q[uid_pos] <- rowsum(r^2, id_pos, reorder = TRUE)[, 1]
    J[uid_pos] <- rowsum((lambda - 1) * logy, id_pos, reorder = TRUE)[, 1]
    k[uid_pos] <- as.integer(table(id_pos))
  }
  list(R = R, Q = Q, J = J, k = k)
}

# Log-density of a subject's positive-recall amounts when u2 ~ N(m, s2) has
# been integrated out (compound-symmetric multivariate normal). `m` may be a
# vector (one value per subject) or an N x Q matrix of node-specific shifts.
.amt_logf <- function(S, m, se2, s2) {
  k   <- S$k
  den <- se2 + k * s2
  base <- -0.5 * k * log(2 * pi) -
    0.5 * ((k - 1) * log(se2) + log(den)) + S$J
  Rm   <- S$R - k * m
  quad <- S$Q - 2 * m * S$R + k * m^2 - s2 * Rm^2 / den
  base - 0.5 * quad / se2
}

.logsumexp_rows <- function(M) {
  mx <- apply(M, 1L, max)
  mx + log(rowSums(exp(M - mx)))
}

.gh_std <- function(n) {
  gh <- gh_nodes_weights(n)
  list(x = gh$nodes, lw = log(gh$weights))
}

# Mode of the per-subject log-integrand in u1 (vectorised Newton with step
# halving) and the marginal log-likelihood by adaptive Gauss-Hermite.
#
#   g(u) = sum_t log Bern(yn_t; plogis(eta_t + u)) + A_i(c u) + log phi(u; 0, sig1^2)
#
# `amt` is NULL for the probability-only model.
.prob_agq <- function(eta, sgn, yn, id, N, sig1, amt = NULL, c_shift = 0,
                      se2 = 1, s2 = 0, n_quad = 11L, mode_only = FALSE) {
  sig1_2 <- sig1^2
  g_at <- function(u) {
    lb <- .by_subj(stats::plogis(sgn * (eta + u[id]), log.p = TRUE), id)
    a  <- if (is.null(amt)) 0 else .amt_logf(amt, c_shift * u, se2, s2)
    lb + a - 0.5 * u^2 / sig1_2
  }
  d_at <- function(u) {
    p  <- stats::plogis(eta + u[id])
    g1 <- .by_subj(yn - p, id) - u / sig1_2
    g2 <- -.by_subj(p * (1 - p), id) - 1 / sig1_2
    if (!is.null(amt)) {
      den <- se2 + amt$k * s2
      g1  <- g1 + c_shift * (amt$R - amt$k * c_shift * u) / den
      g2  <- g2 - c_shift^2 * amt$k / den
    }
    list(g1 = g1, g2 = g2)
  }

  u <- numeric(N)
  g_old <- g_at(u)
  for (it in seq_len(50L)) {
    d    <- d_at(u)
    step <- -d$g1 / d$g2
    u_new <- u + step
    g_new <- g_at(u_new)
    for (h in seq_len(20L)) {
      bad <- !is.finite(g_new) | g_new < g_old - 1e-10
      if (!any(bad)) break
      step[bad]  <- step[bad] / 2
      u_new[bad] <- u[bad] + step[bad]
      g_new[bad] <- g_at(u_new)[bad]
    }
    u <- u_new
    g_old <- g_new
    if (max(abs(step)) < 1e-9) break
  }
  d <- d_at(u)
  if (mode_only) return(list(mode = u, curv = d$g2))

  shat <- 1 / sqrt(-d$g2)
  gh   <- .gh_std(n_quad)
  U    <- u + sqrt(2) * outer(shat, gh$x)          # N x Q
  lb   <- .by_subj(stats::plogis(sgn * (eta + U[id, , drop = FALSE]),
                                 log.p = TRUE), id)
  a    <- if (is.null(amt)) 0 else .amt_logf(amt, c_shift * U, se2, s2)
  G    <- lb + a - 0.5 * U^2 / sig1_2
  G    <- G + matrix(gh$lw + gh$x^2, N, length(gh$x), byrow = TRUE)
  log(sqrt(2) * shat) + .logsumexp_rows(G) - 0.5 * log(2 * pi * sig1_2)
}


# ---------------------------------------------------------------------------
# Model-specific per-subject log-likelihoods. `d` is the design list built by
# .nci_design(); `p` is a named parameter vector.
# ---------------------------------------------------------------------------

.unpack <- function(p, d) {
  nm <- names(p)
  list(
    alpha  = p[startsWith(nm, "P.")],
    beta   = p[startsWith(nm, "A.")],
    se2    = if ("A_LOGSDE" %in% nm) exp(2 * p[["A_LOGSDE"]]) else NA_real_,
    su2    = if ("A_LOGSDU2" %in% nm) exp(2 * p[["A_LOGSDU2"]]) else NA_real_,
    lambda = if ("A_LAMBDA" %in% names(p)) p[["A_LAMBDA"]] else d$lambda_fixed,
    sig1   = if ("P_LOGSDU1" %in% names(p)) exp(p[["P_LOGSDU1"]]) else NA_real_,
    rho    = if ("Z_U" %in% names(p)) tanh(p[["Z_U"]]) else 0
  )
}

.amt_part_sums <- function(q, d) {
  z <- boxcox_transform(d$y_pos, q$lambda)
  r <- z - as.numeric(d$X2 %*% q$beta)
  .amt_sums(r, d$logy_pos, q$lambda, d$id_pos, d$uid_pos, d$N)
}

.ll_subject <- function(p, d, part) {
  q <- .unpack(p, d)
  if (part == "amount") {
    S <- .amt_part_sums(q, d)
    return(.amt_logf(S, 0, q$se2, q$su2))
  }
  eta <- as.numeric(d$X1 %*% q$alpha)
  if (part == "prob") {
    return(.prob_agq(eta, d$sgn, d$yn, d$id, d$N, q$sig1, n_quad = d$n_quad))
  }
  # correlated two-part model
  S  <- .amt_part_sums(q, d)
  su <- sqrt(q$su2)
  .prob_agq(eta, d$sgn, d$yn, d$id, d$N, q$sig1, amt = S,
            c_shift = q$rho * su / q$sig1, se2 = q$se2,
            s2 = q$su2 * (1 - q$rho^2), n_quad = d$n_quad)
}

# Weighted objective (normalised by the total weight so its scale does not
# depend on the magnitude of the survey weights; the maximiser is unchanged).
.neg_wll <- function(p, d, part) {
  ll <- tryCatch(.ll_subject(p, d, part), error = function(e) NA_real_)
  val <- -sum(d$w * ll) / d$sum_w
  if (!is.finite(val)) 1e100 else val
}

.num_grad <- function(f, p, ...) {
  h <- 1e-5 * pmax(1, abs(p))
  vapply(seq_along(p), function(j) {
    e <- numeric(length(p)); e[j] <- h[j]
    (f(p + e, ...) - f(p - e, ...)) / (2 * h[j])
  }, numeric(1))
}

# Maximise the weighted log-likelihood with nlminb (bounded: A_LAMBDA >= 0.01).
.nlmixed_optim <- function(p0, d, part, verbose = FALSE) {
  lower <- rep(-Inf, length(p0)); upper <- rep(Inf, length(p0))
  names(lower) <- names(upper) <- names(p0)
  if ("A_LAMBDA" %in% names(p0)) {
    lower["A_LAMBDA"] <- 0.01
    p0["A_LAMBDA"] <- max(p0["A_LAMBDA"], 0.01)
  }
  # keep log-SDs and Fisher-z in a sane region
  for (nm in intersect(c("A_LOGSDE", "A_LOGSDU2", "P_LOGSDU1"), names(p0))) {
    lower[nm] <- -10; upper[nm] <- 10
  }
  if ("Z_U" %in% names(p0)) { lower["Z_U"] <- -10; upper["Z_U"] <- 10 }

  f0 <- .neg_wll(p0, d, part)
  if (f0 >= 1e100) {
    stop("The ", part, " likelihood is not finite at the starting values.",
         call. = FALSE)
  }
  obj <- function(p) .neg_wll(stats::setNames(p, names(p0)), d, part)
  grd <- function(p) .num_grad(obj, p)
  run <- function(st) stats::nlminb(
    start = st, objective = obj, gradient = grd,
    lower = lower, upper = upper,
    control = list(eval.max = 5000, iter.max = 2000, rel.tol = 1e-14,
                   x.tol = 1e-12)
  )
  opt <- run(p0)
  # nlminb can stop early on flat ridges (e.g. along lambda); restart from
  # the solution until the objective no longer improves.
  for (k in seq_len(5L)) {
    opt2 <- run(opt$par)
    improved <- opt$objective - opt2$objective
    if (opt2$objective <= opt$objective) opt <- opt2
    if (improved < 1e-11) break
  }
  par <- stats::setNames(opt$par, names(p0))
  gr  <- grd(opt$par)
  at_bound <- par <= lower + 1e-8 | par >= upper - 1e-8
  max_grad <- max(abs(gr[!at_bound]), 0)
  # Converged if nlminb reports success, or the (normalised) gradient is
  # negligible relative to the objective's own scale.
  converged <- opt$objective < 1e100 &&
    (opt$convergence == 0 || max_grad < 1e-5 * max(1, abs(opt$objective)))
  if (verbose) {
    message(sprintf("    %s part: max |gradient| = %.2e%s",
                    part, max_grad,
                    if (converged) "" else "  ** NOT CONVERGED **"))
  }
  list(par = par, value = opt$objective, converged = converged,
       max_grad = max_grad, message = opt$message, iterations = opt$iterations)
}


# ---------------------------------------------------------------------------
# Empirical Bayes estimates of the random effects (posterior modes, as in
# NLMIXED `random ... out=`).
# ---------------------------------------------------------------------------
.ebe <- function(p, d, part) {
  q <- .unpack(p, d)
  u1 <- rep(NA_real_, d$N); u2 <- rep(NA_real_, d$N)
  if (part %in% c("amount", "corr")) S <- .amt_part_sums(q, d)
  if (part == "amount") {
    u2 <- q$su2 * S$R / (q$se2 + S$k * q$su2)
  } else if (part == "prob") {
    eta <- as.numeric(d$X1 %*% q$alpha)
    u1 <- .prob_agq(eta, d$sgn, d$yn, d$id, d$N, q$sig1, mode_only = TRUE)$mode
  } else {
    eta <- as.numeric(d$X1 %*% q$alpha)
    su  <- sqrt(q$su2)
    cc  <- q$rho * su / q$sig1
    s2  <- q$su2 * (1 - q$rho^2)
    # Conditional on u1, u2 is Gaussian, so the joint mode in u1 equals the
    # mode of the u2-marginalised integrand; u2 then follows in closed form.
    u1 <- .prob_agq(eta, d$sgn, d$yn, d$id, d$N, q$sig1, amt = S,
                    c_shift = cc, se2 = q$se2, s2 = s2, mode_only = TRUE)$mode
    m  <- cc * u1
    u2 <- m + s2 * (S$R - S$k * m) / (q$se2 + S$k * s2)
  }
  list(u1 = u1, u2 = u2)
}
