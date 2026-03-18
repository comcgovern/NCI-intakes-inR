# BRR tests use a small dataset and just K=4 replicates to keep runtime low.

make_brr_data <- function(n_subjects = 300, seed = 55) {
  set.seed(seed)
  n_recalls <- sample(1:2, n_subjects, replace = TRUE, prob = c(0.4, 0.6))
  total     <- sum(n_recalls)
  sids      <- rep(seq_len(n_subjects), times = n_recalls)
  reps      <- unlist(lapply(n_recalls, seq_len))
  u_i       <- rep(stats::rnorm(n_subjects, 0, 2), times = n_recalls)
  t_intake  <- 15 + u_i + stats::rnorm(total, 0, sqrt(2.5))
  intake    <- pmax(boxcox_inverse(t_intake, 0.30), 1)

  data.frame(SEQN = sids, day = reps, sodium_mg = intake,
             stringsAsFactors = FALSE)
}

make_toy_brr_weights <- function(data, subject_var = "SEQN", K = 4, seed = 1) {
  set.seed(seed)
  subj_ids <- unique(data[[subject_var]])
  n        <- length(subj_ids)
  # Simple perturbation: randomly scale weights near 1.0
  wt_mat <- matrix(
    stats::runif(n * K, 0.7, 1.3),
    nrow = n, ncol = K
  )
  df <- as.data.frame(wt_mat)
  colnames(df) <- paste0("brr_wt_", seq_len(K))
  cbind(data.frame(SEQN = subj_ids, stringsAsFactors = FALSE), df)
}

test_that("brr_usual_intake returns a brr_result object", {
  dat <- make_brr_data()
  rw  <- make_toy_brr_weights(dat)

  res <- suppressMessages(brr_usual_intake(
    data              = dat,
    replicate_weights = rw,
    subject_var       = "SEQN",
    mixtran_args      = list(
      intake_var  = "sodium_mg",
      subject_var = "SEQN",
      repeat_var  = "day",
      model_type  = "amount",
      lambda      = 0.30,
      verbose     = FALSE
    ),
    distrib_args = list(n_sims = 30, seed = 7),
    verbose      = FALSE
  ))

  expect_s3_class(res, "brr_result")
  expect_equal(res$K, 4L)
  expect_true(!is.null(res$estimates))
  expect_true(!is.null(res$se))
  expect_true(all(res$se >= 0))
})

test_that("brr_to_df returns a data frame with expected columns", {
  dat <- make_brr_data(n_subjects = 200)
  rw  <- make_toy_brr_weights(dat, K = 4)

  res <- suppressMessages(brr_usual_intake(
    data              = dat,
    replicate_weights = rw,
    subject_var       = "SEQN",
    mixtran_args      = list(
      intake_var  = "sodium_mg",
      subject_var = "SEQN",
      repeat_var  = "day",
      model_type  = "amount",
      lambda      = 0.30,
      verbose     = FALSE
    ),
    distrib_args = list(n_sims = 20, seed = 8),
    verbose      = FALSE
  ))

  df <- brr_to_df(res)
  expect_true(is.data.frame(df))
  expect_true(all(c("statistic", "estimate", "se", "ci_lower", "ci_upper") %in% names(df)))
  expect_true(nrow(df) > 0)
})

test_that("brr_usual_intake errors on mismatched replicate weight rows", {
  dat <- make_brr_data(n_subjects = 100)
  # Wrong number of rows
  bad_wt <- matrix(runif(50 * 4), nrow = 50, ncol = 4)

  expect_error(
    suppressMessages(brr_usual_intake(
      data              = dat,
      replicate_weights = bad_wt,
      subject_var       = "SEQN",
      mixtran_args      = list(
        intake_var  = "sodium_mg",
        subject_var = "SEQN",
        repeat_var  = "day",
        model_type  = "amount",
        lambda      = 0.30,
        verbose     = FALSE
      ),
      distrib_args = list(n_sims = 10, seed = 1),
      verbose      = FALSE
    )),
    regexp = "unique subjects"
  )
})
