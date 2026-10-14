# Helper: generate episodically consumed food data with known parameters
make_episodic_data <- function(n_subjects    = 700,
                                true_alpha0   = 0.5,
                                true_sigma2_v1 = 1.0,
                                true_beta0    = 3.0,
                                true_sigma2_v2 = 1.5,
                                true_sigma2_e  = 1.0,
                                true_lambda    = 0.25,
                                seed           = 99) {
  set.seed(seed)

  n_recalls <- sample(1:2, n_subjects, replace = TRUE, prob = c(0.35, 0.65))
  total     <- sum(n_recalls)

  subject_ids <- rep(seq_len(n_subjects), times = n_recalls)
  repeat_nums <- unlist(lapply(n_recalls, seq_len))

  v1 <- rep(stats::rnorm(n_subjects, 0, sqrt(true_sigma2_v1)), times = n_recalls)
  v2 <- rep(stats::rnorm(n_subjects, 0, sqrt(true_sigma2_v2)), times = n_recalls)

  p_consume <- stats::plogis(true_alpha0 + v1)
  consumed  <- stats::rbinom(total, 1, p_consume)

  t_amount <- true_beta0 + v2 + stats::rnorm(total, 0, sqrt(true_sigma2_e))
  amount   <- pmax(boxcox_inverse(t_amount, true_lambda), 0.1)

  intake <- consumed * amount

  data.frame(
    SEQN         = subject_ids,
    day          = repeat_nums,
    wholegrain_g = intake,
    stringsAsFactors = FALSE
  )
}

test_that("uncorrelated two-part model recovers variance components", {
  dat <- make_episodic_data(n_subjects = 1000)
  fit <- suppressMessages(mixtran(
    data       = dat,
    intake_var = "wholegrain_g",
    subject_var = "SEQN",
    repeat_var  = "day",
    model_type  = "uncorr",
    lambda      = 0.25,
    verbose     = FALSE
  ))

  expect_s3_class(fit, "mixtran_fit")
  expect_equal(fit$model_type, "uncorr")
  expect_equal(fit$rho, 0)

  v1_err <- abs(fit$sigma2_v1 - 1.0) / 1.0
  v2_err <- abs(fit$sigma2_v2 - 1.5) / 1.5
  e_err  <- abs(fit$sigma2_e  - 1.0) / 1.0

  # sigma2_v1 is weakly identified with 1-2 binary recalls per person; the
  # MLE for this dataset is 0.66 (identical to lme4::glmer with nAGQ = 25)
  expect_lte(v1_err, 0.40, label = "sigma2_v1 within 40%")
  expect_lte(v2_err, 0.25, label = "sigma2_v2 within 25%")
  expect_lte(e_err,  0.25, label = "sigma2_e within 25%")
})

test_that("uncorrelated two-part model returns required fields", {
  dat <- make_episodic_data(n_subjects = 300)
  fit <- suppressMessages(mixtran(
    data       = dat,
    intake_var = "wholegrain_g",
    subject_var = "SEQN",
    repeat_var  = "day",
    model_type  = "uncorr",
    lambda      = 0.25,
    verbose     = FALSE
  ))

  expect_true(!is.null(fit$alpha))
  expect_true(!is.null(fit$beta))
  expect_true(!is.null(fit$sigma2_v1))
  expect_true(!is.null(fit$sigma2_v2))
  expect_true(!is.null(fit$sigma2_e))
  expect_equal(fit$rho, 0)
  expect_true(!is.null(fit$predicted))
})

test_that("two-part model returns NULL with warning when no subjects have 2+ positive recalls", {
  # All single-recall data with ubiquitous intake — no within-person info
  dat <- data.frame(
    id     = 1:50,
    day    = 1,
    intake = stats::rnorm(50, 100, 10),
    stringsAsFactors = FALSE
  )

  # Make about half zeros (one recall per person — no repeated measures)
  set.seed(99)
  dat$intake <- dat$intake * stats::rbinom(50, 1, 0.5)

  # Default (skip_if_empty = TRUE): should warn and return NULL
  result <- expect_warning(
    mixtran(
      data       = dat,
      intake_var = "intake",
      subject_var = "id",
      repeat_var  = "day",
      model_type  = "uncorr",
      lambda      = 0.3,
      verbose     = FALSE
    ),
    "2\\+ positive recalls"
  )
  expect_null(result)

  # skip_if_empty = FALSE: stops with an error
  expect_error(
    mixtran(
      data          = dat,
      intake_var    = "intake",
      subject_var   = "id",
      repeat_var    = "day",
      model_type    = "uncorr",
      lambda        = 0.3,
      skip_if_empty = FALSE,
      verbose       = FALSE
    ),
    "2\\+ positive recalls"
  )
})

test_that("uncorrelated fit equals separate probability and amount fits (MIXTRAN nocorr)", {
  dat <- make_episodic_data(n_subjects = 300)
  fit <- suppressMessages(mixtran(dat, "wholegrain_g", "SEQN", "day",
                                  model_type = "uncorr", verbose = FALSE))
  # the amount part uses consumption days only and estimates lambda
  pos <- dat[dat$wholegrain_g > 0, ]
  famt <- suppressMessages(mixtran(pos, "wholegrain_g", "SEQN", "day",
                                   model_type = "amount", verbose = FALSE))
  expect_equal(unname(fit$beta), unname(famt$beta), tolerance = 1e-4)
  expect_equal(fit$lambda, famt$lambda, tolerance = 1e-4)
  expect_equal(fit$sigma2_v2, famt$sigma2_b, tolerance = 1e-4)
})


test_that("probability part equals lme4::glmer with adaptive quadrature (same likelihood)", {
  skip_if_not_installed("lme4")
  dat <- make_episodic_data(n_subjects = 600, seed = 7)
  fit <- suppressMessages(mixtran(dat, "wholegrain_g", "SEQN", "day",
                                  model_type = "uncorr", lambda = 0.25, verbose = FALSE))
  dat$yn <- as.numeric(dat$wholegrain_g > 0)
  dat$seq_num <- as.numeric(dat$day > 1)
  g <- lme4::glmer(yn ~ seq_num + (1 | SEQN), data = dat, family = stats::binomial,
                   nAGQ = 25)
  expect_equal(fit$sigma2_v1, as.numeric(lme4::VarCorr(g)$SEQN), tolerance = 1e-3)
  expect_equal(unname(fit$alpha), unname(lme4::fixef(g)), tolerance = 1e-3)
})
