# Helper (re-use the episodic generator from uncorr tests, with rho)
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

  # Correlated random effects at subject level
  sd1 <- sqrt(true_sigma2_v1)
  sd2 <- sqrt(true_sigma2_v2)
  z1s <- stats::rnorm(n_subjects)
  z2s <- stats::rnorm(n_subjects)
  v1s <- sd1 * z1s
  v2s <- sd2 * (true_rho * z1s + sqrt(1 - true_rho^2) * z2s)

  v1 <- rep(v1s, times = n_recalls)
  v2 <- rep(v2s, times = n_recalls)

  p_consume <- stats::plogis(0.3 + v1)
  consumed  <- stats::rbinom(total, 1, p_consume)

  t_amount <- 3.0 + v2 + stats::rnorm(total, 0, sqrt(0.8))
  amount   <- pmax(boxcox_inverse(t_amount, true_lambda), 0.1)

  intake <- consumed * amount

  data.frame(
    SEQN         = subject_ids,
    day          = repeat_nums,
    wholegrain_g = intake,
    stringsAsFactors = FALSE
  )
}

test_that("correlated two-part model returns valid rho in (-1, 1)", {
  dat <- make_corr_data()
  fit <- suppressMessages(mixtran(
    data       = dat,
    intake_var = "wholegrain_g",
    subject_var = "SEQN",
    repeat_var  = "day",
    model_type  = "corr",
    lambda      = 0.25,
    verbose     = FALSE
  ))

  expect_s3_class(fit, "mixtran_fit")
  expect_equal(fit$model_type, "corr")
  expect_true(fit$rho > -1 && fit$rho < 1)
})

test_that("correlated model rho estimate has correct sign for positive correlation", {
  dat <- make_corr_data(true_rho = 0.6, seed = 77)
  fit <- suppressMessages(mixtran(
    data       = dat,
    intake_var = "wholegrain_g",
    subject_var = "SEQN",
    repeat_var  = "day",
    model_type  = "corr",
    lambda      = 0.25,
    verbose     = FALSE
  ))

  # Profile likelihood should recover the sign of rho (not necessarily magnitude)
  expect_gt(fit$rho, 0, label = "rho positive when true rho = 0.6")
})

test_that("correlated model stores rho_profile data frame", {
  dat <- make_corr_data(n_subjects = 300, seed = 22)
  fit <- suppressMessages(mixtran(
    data       = dat,
    intake_var = "wholegrain_g",
    subject_var = "SEQN",
    repeat_var  = "day",
    model_type  = "corr",
    lambda      = 0.25,
    verbose     = FALSE
  ))

  expect_true(!is.null(fit$rho_profile))
  expect_true(is.data.frame(fit$rho_profile))
  expect_true(all(c("rho", "loglik") %in% names(fit$rho_profile)))
})
