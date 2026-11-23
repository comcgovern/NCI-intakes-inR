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
#
# With `moments = TRUE` it also returns, under the normalised quadrature
# weights pi_ik (the posterior of u1 given the subject's data), E[u], E[u^2]
# the record-level E[yn_t - p_t], and E[g'(u)] and E[(u - mode) g'(u)]:
# all that .ll_grad() needs.
.prob_agq <- function(eta, sgn, yn, id, N, sig1, amt = NULL, c_shift = 0,
                      se2 = 1, s2 = 0, n_quad = 11L, mode_only = FALSE,
                      moments = FALSE) {
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
  lse  <- .logsumexp_rows(G)
  ll   <- log(sqrt(2) * shat) + lse - 0.5 * log(2 * pi * sig1_2)
  if (!moments) return(ll)
  Pw <- exp(G - lse)                                # N x Q, rows sum to 1
  P  <- stats::plogis(eta + U[id, , drop = FALSE])  # records x Q
  # g'(u) at the nodes
  Gp <- .by_subj(yn - P, id) - U / sig1_2
  if (!is.null(amt)) {
    Gp <- Gp + c_shift * (amt$R - amt$k * c_shift * U) / (se2 + amt$k * s2)
  }
  list(ll = ll, mode = u, shat = shat,
       Eu = rowSums(Pw * U), Eu2 = rowSums(Pw * U^2),
       resid = yn - rowSums(Pw[id, , drop = FALSE] * P),
       Egp = rowSums(Pw * Gp), EUgp = rowSums(Pw * (U - u) * Gp))
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

# d T(y; lambda) / d lambda for the Box-Cox transform.
.boxcox_dlambda <- function(y, logy, lambda) {
  if (abs(lambda) < 1e-10) return(0.5 * logy^2)
  yl <- y^lambda
  (yl * logy) / lambda - (yl - 1) / lambda^2
}

# Per-subject log-likelihood (`ll`) and its analytic gradient (`grad`, an
# N x length(p) matrix with columns named as `p`).
#
# The amount-only likelihood is closed form. For the two-part models this is
# the exact gradient of the adaptive Gauss-Hermite approximation computed by
# .prob_agq(),
#
#   l~ = log(sqrt(2) s) + log sum_k w_k exp(x_k^2 + g(m + sqrt(2) s x_k))
#        - log(sqrt(2 pi) sig1),
#
# including its dependence on each subject's mode m and scale
# s = (-g''(m))^(-1/2), by the implicit-function theorem:
#   dm/dtheta = -(d g'/dtheta) / g''(m),
#   ds/dtheta = s^3 / 2 * (d g''/dtheta + g'''(m) dm/dtheta),
#   dl~/dtheta = E[d g/dtheta] + (ds/dtheta) (1 + E[(u - m) g'(u)]) / s
#                + (dm/dtheta) E[g'(u)],
# where E[] averages over the nodes with the normalised quadrature weights.
# Each partial derivative of the amount log-density .amt_logf() is a
# polynomial of degree <= 2 in its shift c * u1, so E[d g/dtheta] needs only
# the posterior moments E[u1] and E[u1^2].
.ll_grad <- function(p, d, part) {
  q  <- .unpack(p, d)
  nm <- names(p)
  N  <- d$N
  ia <- startsWith(nm, "P.")
  ib <- startsWith(nm, "A.")
  il <- nm == "A_LAMBDA"
  up <- d$uid_pos
  G  <- matrix(0, N, length(p), dimnames = list(NULL, nm))

  if (part != "prob") {
    z <- boxcox_transform(d$y_pos, q$lambda)
    r <- z - as.numeric(d$X2 %*% q$beta)
    S <- .amt_sums(r, d$logy_pos, q$lambda, d$id_pos, d$uid_pos, N)
    k <- S$k; R <- S$R; se2 <- q$se2
    # dR/dbeta and dR/dlambda, by subject (rows `up`)
    dR_db <- -rowsum(d$X2, d$id_pos, reorder = TRUE)
    if (any(il)) {
      dz <- .boxcox_dlambda(d$y_pos, d$logy_pos, q$lambda)
      dR_dl <- rowsum(dz, d$id_pos, reorder = TRUE)[, 1]
    }
  }

  # Posterior-mean gradient of the amount log-density A(m) for the given
  # within-subject variance s2 of u2 and moments of its shift m. With
  # Rm = R - k m and den = se2 + k s2:
  #   dA/dR   = (m se2 + s2 R) / (den se2),   dA/dQ = -1 / (2 se2),
  #   dA/ds2  = -k / (2 den) + Rm^2 / (2 den^2),
  #   dA/dse2 = -(k-1)/(2 se2) - 1/(2 den) - s2 Rm^2/(2 den^2 se2)
  #             + quad / (2 se2^2).
  # Fills the A.* , A_LAMBDA and A_LOGSDE columns of G; returns E[dA/ds2].
  amount_grad <- function(s2, Em, Em2) {
    den   <- se2 + k * s2
    ERm2  <- R^2 - 2 * k * R * Em + k^2 * Em2
    Equad <- S$Q - 2 * R * Em + k * Em2 - s2 * ERm2 / den
    dR    <- (Em * se2 + s2 * R) / (den * se2)
    # dA/dr_j = dR - r_j / se2 for each positive recall j
    G[up, ib] <<- dR[up] * dR_db +
      rowsum(r * d$X2, d$id_pos, reorder = TRUE) / se2
    if (any(il)) {
      G[up, il] <<- dR[up] * dR_dl -
        rowsum(r * dz, d$id_pos, reorder = TRUE)[, 1] / se2 +
        rowsum(d$logy_pos, d$id_pos, reorder = TRUE)[, 1]
    }
    G[, "A_LOGSDE"] <<- 2 * se2 * (-0.5 * (k - 1) / se2 - 0.5 / den -
      0.5 * s2 * ERm2 / (den^2 * se2) + 0.5 * Equad / se2^2)
    -0.5 * k / den + 0.5 * ERm2 / den^2
  }

  if (part == "amount") {
    G[, "A_LOGSDU2"] <- 2 * q$su2 * amount_grad(q$su2, 0, 0)
    return(list(ll = .amt_logf(S, 0, se2, q$su2), grad = G))
  }

  # --- two-part models ---
  corr <- part == "corr"
  eta  <- as.numeric(d$X1 %*% q$alpha)
  sig1 <- q$sig1
  if (corr) {
    su  <- sqrt(q$su2)
    cc  <- q$rho * su / sig1
    s2  <- q$su2 * (1 - q$rho^2)
    den <- se2 + k * s2
    M <- .prob_agq(eta, d$sgn, d$yn, d$id, N, sig1, amt = S, c_shift = cc,
                   se2 = se2, s2 = s2, n_quad = d$n_quad, moments = TRUE)
  } else {
    M <- .prob_agq(eta, d$sgn, d$yn, d$id, N, sig1, n_quad = d$n_quad,
                   moments = TRUE)
  }
  um <- M$mode
  isg <- nm == "P_LOGSDU1"

  # (a) posterior mean of d g / d theta at fixed nodes
  G[, ia]  <- rowsum(M$resid * d$X1, d$id, reorder = TRUE)
  G[, isg] <- M$Eu2 / sig1^2 - 1

  # (b) D1 = d g'/d theta and D2 = d g''/d theta at the mode, where
  #   g'  = sum(yn - p) + c (R - k c u) / den - u / sig1^2,
  #   g'' = -sum p (1 - p) - c^2 k / den - 1 / sig1^2.
  D1 <- D2 <- matrix(0, N, length(p))
  ph <- stats::plogis(eta + um[d$id])
  w1 <- ph * (1 - ph)
  w2 <- w1 * (1 - 2 * ph)
  D1[, ia]  <- -rowsum(w1 * d$X1, d$id, reorder = TRUE)
  D2[, ia]  <- -rowsum(w2 * d$X1, d$id, reorder = TRUE)
  D1[, isg] <- 2 * um / sig1^2
  D2[, isg] <- 2 / sig1^2
  g3 <- -.by_subj(w2, d$id)

  if (corr) {
    ds2 <- amount_grad(s2, cc * M$Eu, cc^2 * M$Eu2)
    dc  <- (R * M$Eu - k * cc * M$Eu2) / den       # E[dA/dc], A depends on c u
    cRm <- cc * (R - k * cc * um)
    # partials of g' (a1_*) and g'' (a2_*) at the mode
    a1_c  <- (R - 2 * k * cc * um) / den; a2_c  <- -2 * cc * k / den
    a1_s2 <- -cRm * k / den^2;            a2_s2 <- cc^2 * k^2 / den^2
    a1_e2 <- -cRm / den^2;                a2_e2 <- cc^2 * k / den^2
    a1_R  <- cc / den
    D1[up, ib] <- a1_R[up] * dR_db
    if (any(il)) D1[up, il] <- a1_R[up] * dR_dl
    ie <- nm == "A_LOGSDE"
    D1[, ie] <- 2 * se2 * a1_e2; D2[, ie] <- 2 * se2 * a2_e2
    # log sd(u2): c scales with su, s2 with su^2
    iu <- nm == "A_LOGSDU2"
    G[, iu]  <- cc * dc + 2 * s2 * ds2
    D1[, iu] <- cc * a1_c + 2 * s2 * a1_s2
    D2[, iu] <- cc * a2_c + 2 * s2 * a2_s2
    # log sd(u1): c scales with 1 / sig1
    G[, isg]  <- G[, isg]  - cc * dc
    D1[, isg] <- D1[, isg] - cc * a1_c
    D2[, isg] <- D2[, isg] - cc * a2_c
    # Fisher z: drho/dz = 1 - rho^2, dc/drho = su/sig1, ds2/drho = -2 rho su2
    iz <- nm == "Z_U"
    jc <- su / sig1; js <- -2 * q$rho * q$su2; jz <- 1 - q$rho^2
    G[, iz]  <- jz * (jc * dc   + js * ds2)
    D1[, iz] <- jz * (jc * a1_c + js * a1_s2)
    D2[, iz] <- jz * (jc * a2_c + js * a2_s2)
  }

  # (c) dependence of the nodes on theta through the mode and scale
  shat <- M$shat
  dm   <- D1 * shat^2                          # = -D1 / g''(m)
  ds   <- 0.5 * shat^3 * (D2 + g3 * dm)
  G <- G + ds * ((1 + M$EUgp) / shat) + dm * M$Egp
  list(ll = M$ll, grad = G)
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

# Weighted objective and its gradient from one pass of .ll_grad().
.neg_wll_grad <- function(p, d, part) {
  r <- tryCatch(.ll_grad(p, d, part), error = function(e) NULL)
  val <- if (is.null(r)) NA_real_ else -sum(d$w * r$ll) / d$sum_w
  if (!is.finite(val)) return(list(value = 1e100, grad = NULL))
  gr <- -colSums(d$w * r$grad) / d$sum_w
  list(value = val, grad = if (all(is.finite(gr))) gr else NULL)
}

# Maximise the weighted log-likelihood with nlminb (bounded: A_LAMBDA >= 0.01)
# using the analytic gradient (.ll_grad()). With `verbose`, a progress line
# is printed every `progress_every` seconds of a long fit. `control` overrides
# the nlminb control settings.
.nlmixed_optim <- function(p0, d, part, verbose = FALSE, progress_every = 30,
                           control = list()) {
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

  # nlminb asks for the objective and then the gradient at the same point;
  # both come from one likelihood pass, cached here.
  cache <- new.env(parent = emptyenv())
  eval_at <- function(p) {
    if (!identical(p, cache$p)) {
      cache$p <- p
      cache$r <- .neg_wll_grad(stats::setNames(p, names(p0)), d, part)
    }
    cache$r
  }
  obj <- function(p) eval_at(p)$value
  grd <- function(p) {
    g <- eval_at(p)$grad
    if (is.null(g)) {  # not finite: fall back to finite differences
      g <- .num_grad(function(x) .neg_wll(stats::setNames(x, names(p0)), d, part), p)
    }
    track(p)
    g
  }

  # progress reporting (gradient calls ~ iterations)
  t_start <- t_last <- Sys.time()
  n_grad <- 0L
  track <- function(p) {
    n_grad <<- n_grad + 1L
    if (verbose && as.numeric(Sys.time() - t_last, units = "secs") > progress_every) {
      t_last <<- Sys.time()
      rho <- if ("Z_U" %in% names(p0)) {
        sprintf("  rho = %.4f", tanh(p[match("Z_U", names(p0))]))
      } else ""
      message(sprintf("    %s part: %d iterations, %.0f s, -2 log L = %.8g%s",
                      part, n_grad,
                      as.numeric(Sys.time() - t_start, units = "secs"),
                      2 * d$sum_w * cache$r$value, rho))
    }
  }

  ctrl <- utils::modifyList(
    list(eval.max = 5000, iter.max = 2000, rel.tol = 1e-14, x.tol = 1e-12),
    control)
  run <- function(st) stats::nlminb(start = st, objective = obj, gradient = grd,
                                    lower = lower, upper = upper, control = ctrl)
  limit_hit <- function(o) o$convergence != 0 &&
    grepl("limit reached", o$message, fixed = TRUE)

  opt <- run(p0)
  runs <- 1L
  hit_limit <- limit_hit(opt)
  # nlminb can stop early on flat ridges (e.g. along lambda); restart from
  # the solution until the objective no longer improves.
  max_restarts <- 5L
  improved <- Inf
  for (k in seq_len(max_restarts)) {
    opt2 <- run(opt$par)
    runs <- runs + 1L
    hit_limit <- hit_limit || limit_hit(opt2)
    improved <- opt$objective - opt2$objective
    if (opt2$objective <= opt$objective) opt <- opt2
    if (improved < 1e-11) break
  }
  # Still improving materially when the restarts ran out?
  still_improving <- improved >= 1e-8 * max(1, abs(opt$objective))

  par <- stats::setNames(opt$par, names(p0))
  gr  <- eval_at(opt$par)$grad
  if (is.null(gr)) gr <- rep(Inf, length(p0))
  at_bound <- par <= lower + 1e-8 | par >= upper - 1e-8
  max_grad <- max(abs(gr[!at_bound]), 0)
  # Converged if nlminb reports success, or the (normalised) gradient is
  # negligible relative to the objective's own scale.
  converged <- opt$objective < 1e100 &&
    (opt$convergence == 0 || max_grad < 1e-5 * max(1, abs(opt$objective)))
  secs <- as.numeric(Sys.time() - t_start, units = "secs")
  if (verbose) {
    message(sprintf("    %s part: %d iterations in %d run(s), %.1f s, max |gradient| = %.2e%s",
                    part, n_grad, runs, secs, max_grad,
                    if (converged) "" else "  ** NOT CONVERGED **"))
  }
  # (non-convergence itself is reported by mixtran())
  if (hit_limit || still_improving) {
    warning(sprintf(paste0(
      "The %s fit %s after %d iterations (%.0f s); max |gradient| = %.2e. ",
      "Check the estimates, or refit with better starting values via 'start'."),
      part,
      if (hit_limit) "reached nlminb's iteration/evaluation limit"
      else sprintf("was still improving after %d restarts", max_restarts),
      n_grad, secs, max_grad), call. = FALSE)
  }
  list(par = par, value = opt$objective, converged = converged,
       max_grad = max_grad, message = opt$message, iterations = n_grad,
       runs = runs, seconds = secs)
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
