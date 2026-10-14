# Correlated two-part model: joint weighted maximum likelihood (MIXTRAN
# modeltype = corr), random effects integrated by adaptive quadrature.

make_corr_data <- function(n_subjects     = 600,
                           true_rho       = 0.5,
                           true_sigma2_v1 = 1.2,
                           true_sigma2_v2 = 1.0,
                           true_lambda    = 0.25,
                           seed           = 11) {
  set.seed(seed)
  n_recalls <- sample(1:2, n_subjects, replace = TRUE, prob = c(0.35, 0.65))
  total     <- sum(n_recalls)
  subject_ids <- rep(seq_len(n_subjects), times = n_recalls)
  repeat_nums <- unlist(lapply(n_recalls, seq_len))
  sd1 <- sqrt(true_sigma2_v1); sd2 <- sqrt(true_sigma2_v2)
  z1s <- stats::rnorm(n_subjects); z2s <- stats::rnorm(n_subjects)
  v1 <- rep(sd1 * z1s, times = n_recalls)
  v2 <- rep(sd2 * (true_rho * z1s + sqrt(1 - true_rho^2) * z2s), times = n_recalls)
  consumed <- stats::rbinom(total, 1, stats::plogis(0.3 + v1))
  amount   <- pmax(boxcox_inverse(3.0 + v2 + stats::rnorm(total, 0, sqrt(0.8)),
                                  true_lambda), 0.1)
  data.frame(SEQN = subject_ids, day = repeat_nums,
             wholegrain_g = consumed * amount, stringsAsFactors = FALSE)
}

test_that("correlated two-part model returns a valid fit", {
  fit <- suppressMessages(mixtran(make_corr_data(), "wholegrain_g", "SEQN", "day",
                                  model_type = "corr", lambda = 0.25, verbose = FALSE))
  expect_s3_class(fit, "mixtran_fit")
  expect_equal(fit$model_type, "corr")
  expect_true(fit$converged)
  expect_true(fit$rho > -1 && fit$rho < 1)
  expect_equal(fit$cov_u1u2, fit$rho * sqrt(fit$sigma2_v1 * fit$sigma2_v2))
  # the correlated fit can only improve on the uncorrelated likelihood
  expect_lte(fit$m2loglik, fit$m2loglik_uncorr + 1e-6)
})

test_that("joint ML recovers a strong positive rho", {
  fit <- suppressMessages(mixtran(make_corr_data(n_subjects = 1500, true_rho = 0.6, seed = 77),
                                  "wholegrain_g", "SEQN", "day",
                                  model_type = "corr", lambda = 0.25, verbose = FALSE))
  expect_gt(fit$rho, 0.35)
  expect_lt(fit$rho, 0.85)
})

test_that("joint ML gives rho near 0 when random effects are independent", {
  fit <- suppressMessages(mixtran(make_corr_data(n_subjects = 1500, true_rho = 0, seed = 5),
                                  "wholegrain_g", "SEQN", "day",
                                  model_type = "corr", lambda = 0.25, verbose = FALSE))
  expect_lt(abs(fit$rho), 0.25)
})

test_that("corr marginal likelihood matches brute-force 2-D integration", {
  dat <- make_corr_data(n_subjects = 80, seed = 3)
  fit <- suppressMessages(mixtran(dat, "wholegrain_g", "SEQN", "day",
                                  model_type = "corr", verbose = FALSE))
  w <- fit$data
  q <- nciusual:::.unpack(fit$theta, list(lambda_fixed = NA))
  S <- matrix(c(q$sig1^2, rep(q$rho * q$sig1 * sqrt(q$su2), 2), q$su2), 2)
  Si <- solve(S)
  ll <- 0
  for (s in unique(w$subject)[1:15]) {
    r <- w[w$subject == s, ]
    eta <- as.numeric(cbind(1, as.matrix(r[, fit$prob_covs_used])) %*% fit$alpha)
    mu  <- as.numeric(cbind(1, as.matrix(r[, fit$amt_covs_used])) %*% fit$beta)
    f <- function(a, b) {
      l <- sum(stats::plogis((2 * r$consumed - 1) * (eta + a), log.p = TRUE))
      p <- r$consumed == 1
      if (any(p)) l <- l + sum(stats::dnorm(boxcox_transform(r$intake[p], fit$lambda),
                                            mu[p] + b, sqrt(q$se2), log = TRUE) +
                               (fit$lambda - 1) * log(r$intake[p]))
      exp(l - 0.5 * (a^2 * Si[1, 1] + 2 * a * b * Si[1, 2] + b^2 * Si[2, 2]) -
            log(2 * pi) - 0.5 * log(det(S)))
    }
    inner <- function(a) vapply(a, function(aa) stats::integrate(
      function(b) vapply(b, function(bb) f(aa, bb), 1), -Inf, Inf, rel.tol = 1e-10)$value, 1)
    ll <- ll + log(stats::integrate(inner, -Inf, Inf, rel.tol = 1e-10)$value)
  }
  d <- nciusual:::.nci_design(w, fit$prob_covs_used, fit$amt_covs_used, NA, 21, has_prob = TRUE)
  eng <- nciusual:::.ll_subject(fit$theta, d, "corr")[1:15]
  expect_equal(sum(eng), ll, tolerance = 1e-6)
})

test_that("corr model with fewer than 2 paired-positive subjects returns NULL (MIXTRAN stops)", {
  dat <- data.frame(id = rep(1:30, each = 2), day = rep(1:2, 30), y = 0)
  dat$y[c(1, 2)] <- c(5, 6)                 # one subject with 2 positive recalls
  dat$y[seq(3, 60, by = 2)] <- 3             # others: one positive recall
  expect_warning(res <- mixtran(dat, "y", "id", "day", model_type = "corr",
                                lambda = 0.3, verbose = FALSE),
                 "fewer than two subjects")
  expect_null(res)
  expect_error(mixtran(dat, "y", "id", "day", model_type = "corr", lambda = 0.3,
                       skip_if_empty = FALSE, verbose = FALSE),
               "fewer than two subjects")
})

test_that("2-10 paired-positive subjects warns (MIXTRAN 'might be unstable') but fits", {
  set.seed(9)
  dat <- make_corr_data(n_subjects = 200, seed = 9)
  ids <- unique(dat$SEQN)
  pos2 <- names(which(tapply(dat$wholegrain_g > 0, dat$SEQN, sum) >= 2))
  drop <- dat$SEQN %in% as.numeric(pos2[-(1:6)]) & dat$day == 2
  dat2 <- dat[!drop, ]
  expect_warning(fit <- suppressMessages(mixtran(dat2, "wholegrain_g", "SEQN", "day",
                                                 model_type = "uncorr", lambda = 0.25,
                                                 verbose = FALSE)),
                 "might be unstable")
  expect_s3_class(fit, "mixtran_fit")
})

test_that("mixtran notes when >90% of person-days are positive for a two-part model", {
  set.seed(1)
  n <- 150
  dat <- data.frame(id = rep(1:n, each = 2), day = rep(1:2, n),
                    y = stats::rlnorm(2 * n, 3, 0.5))
  dat$y[sample(2 * n, 10)] <- 0
  expect_message(suppressWarnings(mixtran(dat, "y", "id", "day", model_type = "uncorr",
                                          lambda = 0.3, verbose = FALSE)),
                 "may be more suitable")
})

test_that("deprecated engine arguments warn and are ignored", {
  dat <- make_corr_data(n_subjects = 200, seed = 4)
  expect_warning(fit <- suppressMessages(mixtran(dat, "wholegrain_g", "SEQN", "day",
                                                 model_type = "corr", lambda = 0.25,
                                                 corr_engine = "ghq", verbose = FALSE)),
                 "deprecated")
  expect_s3_class(fit, "mixtran_fit")
})

test_that("separate probability covariates are supported (covars_prob)", {
  set.seed(2); n <- 500
  dat <- make_corr_data(n_subjects = n, seed = 2)
  dat$x <- stats::rbinom(n, 1, 0.5)[dat$SEQN]
  fit <- suppressMessages(mixtran(dat, "wholegrain_g", "SEQN", "day", model_type = "uncorr",
                                  covariates = NULL, covariates_prob = "x",
                                  lambda = 0.25, verbose = FALSE))
  expect_true("x" %in% names(fit$alpha))
  expect_false("x" %in% names(fit$beta))
})
