# Helper: fast amount model fit
fit_amount_for_distrib <- function(n_subjects = 600, seed = 42) {
  set.seed(seed)
  n_recalls <- sample(1:2, n_subjects, replace = TRUE, prob = c(0.4, 0.6))
  total     <- sum(n_recalls)
  sids      <- rep(seq_len(n_subjects), times = n_recalls)
  reps      <- unlist(lapply(n_recalls, seq_len))
  u_i       <- rep(stats::rnorm(n_subjects, 0, 2), times = n_recalls)
  t_intake  <- 15 + u_i + stats::rnorm(total, 0, sqrt(2.5))
  intake    <- pmax(boxcox_inverse(t_intake, 0.30), 1)
  wt        <- rep(stats::runif(n_subjects, 0.5, 2.0), times = n_recalls)
  dat       <- data.frame(SEQN = sids, day = reps, sodium_mg = intake,
                           wt = wt, stringsAsFactors = FALSE)

  suppressMessages(mixtran(
    data       = dat,
    intake_var = "sodium_mg",
    subject_var = "SEQN",
    repeat_var  = "day",
    model_type  = "amount",
    weight_var  = "wt",
    lambda      = 0.30,
    verbose     = FALSE
  ))
}

test_that("distrib returns distrib_result with correct structure", {
  fit <- fit_amount_for_distrib()
  res <- suppressMessages(distrib(fit, n_sims = 50, seed = 1))

  expect_s3_class(res, "distrib_result")
  expect_true(!is.null(res$results))
  expect_true("Overall" %in% names(res$results))
  expect_true(!is.null(res$results$Overall$mean))
  expect_true(!is.null(res$results$Overall$sd))
  expect_true(!is.null(res$results$Overall$percentiles))
})

test_that("usual intake distribution is narrower than observed", {
  fit <- fit_amount_for_distrib()
  res <- suppressMessages(distrib(fit, n_sims = 100, seed = 7))

  usual_sd <- res$results$Overall$sd
  obs_sd   <- stats::sd(fit$data$intake)

  expect_lt(usual_sd, obs_sd,
            label = "Usual intake SD < observed SD (within-person variance removed)")
})

test_that("distrib mean is within 10% of theoretical", {
  fit <- fit_amount_for_distrib()
  res <- suppressMessages(distrib(fit, n_sims = 200, seed = 99))

  # Theoretical mean (NCI definition): E[T^{-1}(beta0 + u + e)], with
  # u ~ N(0, sigma2_b) and e ~ N(0, sigma2_w)
  theo_samp  <- boxcox_inverse(
    stats::rnorm(100000, fit$beta[["(Intercept)"]],
                 sqrt(fit$sigma2_b + fit$sigma2_w)),
    fit$lambda
  )
  theo_mean  <- mean(theo_samp)
  est_mean   <- res$results$Overall$mean
  pct_err    <- abs(est_mean - theo_mean) / theo_mean

  expect_lte(pct_err, 0.10, label = "Mean within 10% of theoretical")
})

test_that("distrib cutpoints work correctly", {
  fit <- fit_amount_for_distrib()
  res <- suppressMessages(distrib(
    fit, n_sims = 100, seed = 5,
    cutpoints = c(1500, 2300, 3500)
  ))

  r <- res$results$Overall
  expect_true(!is.null(r$cutpoint_below))
  expect_true(!is.null(r$cutpoint_above))
  expect_equal(length(r$cutpoint_below), 3L)
  # proportions must be in [0,1]
  expect_true(all(r$cutpoint_below >= 0 & r$cutpoint_below <= 1))
  # above + below = 1
  expect_equal(unname(r$cutpoint_below + r$cutpoint_above),
               rep(1, 3), tolerance = 1e-10)
})

test_that("distrib_to_df returns a data frame with one row per subgroup", {
  fit <- fit_amount_for_distrib()
  res <- suppressMessages(distrib(fit, n_sims = 50, seed = 3))
  df  <- distrib_to_df(res)

  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 1L)            # "Overall" only
  expect_true("mean" %in% names(df))
  expect_true("sd"   %in% names(df))
})

test_that("distrib is reproducible with the same seed", {
  fit <- fit_amount_for_distrib()
  r1  <- suppressMessages(distrib(fit, n_sims = 50, seed = 42))
  r2  <- suppressMessages(distrib(fit, n_sims = 50, seed = 42))

  expect_equal(r1$results$Overall$mean, r2$results$Overall$mean)
  expect_equal(r1$results$Overall$percentiles, r2$results$Overall$percentiles)
})

test_that("weighted_quantiles returns correct quantiles", {
  # With equal weights, should match base R quantile type 7
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  w <- rep(1, 10)
  q <- nciusual:::weighted_quantiles(x, w, c(0.1, 0.5, 0.9))

  expect_equal(q[2], 5.5, tolerance = 0.1)   # median
  expect_true(q[1] < q[2] && q[2] < q[3])    # monotone
})

test_that("weighted_quantiles falls back to unweighted when all weights are zero", {
  x <- c(1, 2, 3, 4, 5)
  w_zero <- rep(0, 5)
  q <- expect_warning(
    nciusual:::weighted_quantiles(x, w_zero, c(0.5)),
    "zero or non-finite"
  )
  expect_true(is.finite(q))
  expect_equal(q, 3, tolerance = 0.01)  # unweighted median of 1:5
})

test_that("weighted_quantiles with equal weights matches base R quantile (type 5)", {
  # DISTRIB's mid-weight plotting positions reduce to quantile type 5 for
  # equal weights and distinct values (inside the data range).
  x <- c(3, 1, 4, 7, 5, 9, 2, 6)
  w <- rep(1, length(x))
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  wq  <- nciusual:::weighted_quantiles(x, w, probs)
  ref <- as.numeric(stats::quantile(x, probs, type = 5))
  expect_equal(wq, ref, tolerance = 1e-8)
})

test_that("weighted_quantiles pools ties as the DISTRIB PC macro does", {
  # tied values are one point carrying their combined weight
  x <- c(1, 1, 2, 3); w <- rep(1, 4)
  # positions: value 1 -> 0.25, value 2 -> 0.625, value 3 -> 0.875
  expect_equal(nciusual:::weighted_quantiles(x, w, 0.625), 2)
  expect_equal(nciusual:::weighted_quantiles(x, w, 0.25), 1)
  expect_equal(nciusual:::weighted_quantiles(x, w, 0.75), 2.5)
})

test_that("DISTRIB back-transformation matches the 9-point rule and the exact lognormal mean", {
  # lambda = 0: exact
  expect_equal(nciusual:::backtransform_9pt(1, 0, 0.5, 1e-6), exp(1 + 0.25))
  # lambda > 0: 9-point approximation close to exact quadrature
  mu <- 5; lam <- 0.4; ve <- 0.6
  exact <- integrate(function(e) boxcox_inverse(mu + e, lam) * dnorm(e, 0, sqrt(ve)),
                     -Inf, Inf)$value
  expect_equal(nciusual:::backtransform_9pt(mu, lam, ve, 1e-6), exact, tolerance = 1e-3)
  # bounded below by half the minimum amount
  expect_equal(nciusual:::backtransform_9pt(-50, lam, ve, 2), 1)
})

test_that("distrib averages weekday and weekend usual intake (wkend_prop)", {
  set.seed(3); n <- 400
  sid <- rep(1:n, each = 2); wk <- rbinom(2 * n, 1, 0.4)
  y <- exp(2 + 0.6 * wk + rep(rnorm(n, 0, .5), each = 2) + rnorm(2 * n, 0, .4))
  dat <- data.frame(id = sid, day = rep(1:2, n), y = y, wk = wk)
  fit <- mixtran(dat, "y", "id", "day", weekend_var = "wk", lambda = 0,
                 verbose = FALSE)
  m0 <- distrib(fit, n_sims = 400, seed = 1, wkend_prop = 1e-9)$results$Overall$mean
  m1 <- distrib(fit, n_sims = 400, seed = 1, wkend_prop = 1 - 1e-9)$results$Overall$mean
  mm <- distrib(fit, n_sims = 400, seed = 1)$results$Overall$mean
  expect_equal(mm, 4/7 * m0 + 3/7 * m1, tolerance = 1e-6)
  expect_gt(m1, m0)
})

test_that("factor subgroup variables are handled", {
  fit <- fit_amount_for_distrib(n_subjects = 300)
  fit$predicted$grp <- factor(ifelse(seq_len(nrow(fit$predicted)) %% 2 == 0, "a", "b"))
  res <- distrib(fit, subgroup_var = "grp", n_sims = 20)
  expect_setequal(names(res$results), c("a", "b", "Overall"))
  expect_true(all(is.finite(vapply(res$results, `[[`, numeric(1), "mean"))))
})

test_that("distrib does not change the caller's RNG state", {
  fit <- fit_amount_for_distrib(n_subjects = 100)
  set.seed(42); a <- runif(1)
  set.seed(42); invisible(distrib(fit, n_sims = 5, seed = 1)); b <- runif(1)
  expect_equal(a, b)
})
