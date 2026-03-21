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
  expect_true(is.numeric(fit$rho) && !is.na(fit$rho))
  expect_true(fit$rho >= -1 && fit$rho <= 1)
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

  # Profile likelihood should recover the sign of rho (not necessarily magnitude).
  # If the probability model failed, rho falls back to 0 — still valid, just conservative.
  expect_true(is.numeric(fit$rho) && !is.na(fit$rho))
  expect_gte(fit$rho, 0, label = "rho non-negative when true rho = 0.6")
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

# ---------------------------------------------------------------------------
# Problem 1: skip_if_empty should check paired positive observations for corr
# ---------------------------------------------------------------------------

test_that("skip_if_empty returns NULL for corr when fewer than 10 subjects have paired positive recalls", {
  # Build data with 200 subjects but only 5 have 2 positive recalls;
  # the rest have at most one positive recall so n_paired_pos < 10.
  set.seed(42)
  n_subj <- 200
  # Each subject gets exactly 2 recalls
  subject_ids <- rep(seq_len(n_subj), each = 2)
  repeat_nums <- rep(1:2, times = n_subj)

  # By default everyone has zero intake ...
  intake <- rep(0, n_subj * 2)
  # ... except 5 subjects who have positive intake on both their recalls
  paired_pos_idx <- seq_len(5)
  pos_rows <- which(subject_ids %in% paired_pos_idx)
  intake[pos_rows] <- stats::runif(length(pos_rows), 1, 10)
  # A few more subjects have positive intake on only one recall (n_paired_pos stays at 5)
  single_pos_subjs <- 6:30
  for (s in single_pos_subjs) {
    row1 <- which(subject_ids == s)[1]
    intake[row1] <- stats::runif(1, 1, 10)
  }

  dat <- data.frame(
    SEQN         = subject_ids,
    day          = repeat_nums,
    wholegrain_g = intake,
    stringsAsFactors = FALSE
  )

  expect_warning(
    result <- mixtran(
      data        = dat,
      intake_var  = "wholegrain_g",
      subject_var = "SEQN",
      repeat_var  = "day",
      model_type  = "corr",
      lambda      = 0.25,
      verbose     = FALSE,
      skip_if_empty = TRUE
    ),
    regexp = "only 5 subject"
  )
  expect_null(result)
})

test_that("skip_if_empty corr check does NOT fire when >=10 paired-positive subjects exist", {
  # make_corr_data() produces ~65% of subjects with 2 recalls and realistic
  # consumption — well above the 10-subject threshold.
  dat <- make_corr_data(n_subjects = 200, seed = 55)
  # Capture all warnings and assert none contain "only N subject(s)"
  caught_warnings <- character(0)
  withCallingHandlers(
    suppressMessages(mixtran(
      data        = dat,
      intake_var  = "wholegrain_g",
      subject_var = "SEQN",
      repeat_var  = "day",
      model_type  = "corr",
      lambda      = 0.25,
      verbose     = FALSE,
      skip_if_empty = TRUE
    )),
    warning = function(w) {
      caught_warnings <<- c(caught_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(
    any(grepl("only.*subject", caught_warnings, ignore.case = TRUE)),
    info = "paired-positive threshold warning should not fire when n_paired_pos >= 10"
  )
})

# ---------------------------------------------------------------------------
# Problem 2: all-NA prob_ranef() detected and reported with clear warning
# ---------------------------------------------------------------------------

test_that("fit_twopart_corr emits specific warning and sets convergence_warning when prob_ranef is all NA", {
  # Patch prob_ranef to return all NAs by temporarily replacing it in the
  # package namespace so fit_twopart_corr() sees the broken version.
  ns <- asNamespace("nciusual")
  orig_fn <- get("prob_ranef", envir = ns)
  on.exit(assignInNamespace("prob_ranef", orig_fn, ns = ns), add = TRUE)
  assignInNamespace("prob_ranef", function(prob_fit) {
    re <- orig_fn(prob_fit)
    re[] <- NA_real_
    re
  }, ns = ns)

  dat <- make_corr_data(n_subjects = 200, seed = 33)
  expect_warning(
    fit <- suppressMessages(mixtran(
      data        = dat,
      intake_var  = "wholegrain_g",
      subject_var = "SEQN",
      repeat_var  = "day",
      model_type  = "corr",
      lambda      = 0.25,
      verbose     = FALSE
    )),
    regexp = "probability sub-model returned all-NA"
  )
  # Should fall back to the uncorrelated result with rho = 0 and a flag
  expect_equal(fit$rho, 0)
  expect_equal(fit$convergence_warning, "prob_model_all_na_ranef")
})

# ---------------------------------------------------------------------------
# Problem 3: >90% positive intake emits a diagnostic message
# ---------------------------------------------------------------------------

test_that("mixtran emits diagnostic message when >90% of person-days are positive for two-part model", {
  # Create data where nearly everyone consumes the food on every recall
  set.seed(77)
  n_subj <- 150
  subject_ids <- rep(seq_len(n_subj), each = 2)
  repeat_nums <- rep(1:2, times = n_subj)
  # 95% positive rate
  intake <- ifelse(stats::runif(n_subj * 2) < 0.95, stats::runif(n_subj * 2, 1, 50), 0)
  dat <- data.frame(
    SEQN         = subject_ids,
    day          = repeat_nums,
    wholegrain_g = intake,
    stringsAsFactors = FALSE
  )

  # uncorr model
  expect_message(
    suppressWarnings(mixtran(
      data        = dat,
      intake_var  = "wholegrain_g",
      subject_var = "SEQN",
      repeat_var  = "day",
      model_type  = "uncorr",
      lambda      = 0.25,
      verbose     = FALSE
    )),
    regexp = "model_type.*amount.*more suitable|amount.*more appropriate",
    ignore.case = TRUE
  )

  # corr model also triggers the message
  expect_message(
    suppressWarnings(mixtran(
      data        = dat,
      intake_var  = "wholegrain_g",
      subject_var = "SEQN",
      repeat_var  = "day",
      model_type  = "corr",
      lambda      = 0.25,
      verbose     = FALSE
    )),
    regexp = "model_type.*amount.*more suitable|amount.*more appropriate",
    ignore.case = TRUE
  )
})

test_that("mixtran does NOT emit the >90% diagnostic for amount model", {
  set.seed(88)
  n_subj <- 100
  subject_ids <- rep(seq_len(n_subj), each = 2)
  repeat_nums <- rep(1:2, times = n_subj)
  # 95% positive but using amount model — no diagnostic expected
  intake <- ifelse(stats::runif(n_subj * 2) < 0.95, stats::runif(n_subj * 2, 1, 50), 0)
  dat <- data.frame(
    SEQN         = subject_ids,
    day          = repeat_nums,
    wholegrain_g = intake,
    stringsAsFactors = FALSE
  )
  # Capture messages; none should mention the two-part diagnostic
  caught_messages <- character(0)
  withCallingHandlers(
    suppressWarnings(mixtran(
      data        = dat,
      intake_var  = "wholegrain_g",
      subject_var = "SEQN",
      repeat_var  = "day",
      model_type  = "amount",
      lambda      = 0.25,
      verbose     = FALSE
    )),
    message = function(m) {
      caught_messages <<- c(caught_messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(
    any(grepl("amount.*more suitable|more appropriate", caught_messages, ignore.case = TRUE)),
    info = "amount model should not trigger the two-part ubiquitous-food diagnostic"
  )
})
