# Tests for combine_nhanes_cycles()

make_fake_nhanes_cycle <- function(n = 100, seed = 1, weight_base = 5000) {
  set.seed(seed)
  data.frame(
    SEQN    = seq_len(n),
    WTDRD1  = stats::runif(n, weight_base * 0.5, weight_base * 1.5),
    DRXILINE = sample(1:2, n, replace = TRUE),
    DR1TKCAL = stats::rnorm(n, 2000, 400),
    SDMVSTRA = sample(1:15, n, replace = TRUE),
    SDMVPSU  = sample(1:2, n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# -------------------------------------------------------------------------

test_that("combine_nhanes_cycles returns a data frame", {
  c1 <- make_fake_nhanes_cycle(80,  seed = 1)
  c2 <- make_fake_nhanes_cycle(90,  seed = 2)
  c3 <- make_fake_nhanes_cycle(100, seed = 3)
  combined <- combine_nhanes_cycles(
    list(c1, c2, c3),
    weight_var  = "WTDRD1",
    subject_var = "SEQN"
  )
  expect_s3_class(combined, "data.frame")
  expect_equal(nrow(combined), 80 + 90 + 100)
})

test_that("combine_nhanes_cycles divides weights by n_cycles", {
  c1 <- make_fake_nhanes_cycle(50, seed = 10)
  c2 <- make_fake_nhanes_cycle(50, seed = 11)

  original_weights_c1 <- c1$WTDRD1

  combined <- combine_nhanes_cycles(
    list(c1, c2),
    weight_var  = "WTDRD1",
    subject_var = "SEQN"
  )

  # First 50 rows correspond to cycle 1 (row order preserved)
  cycle1_rows <- combined[combined$cycle == "cycle1", ]
  expect_equal(cycle1_rows$WTDRD1, original_weights_c1 / 2, tolerance = 1e-10)
})

test_that("combine_nhanes_cycles makes unique subject IDs by default", {
  c1 <- make_fake_nhanes_cycle(50, seed = 20)
  c2 <- make_fake_nhanes_cycle(50, seed = 21)

  # Both cycles have SEQN 1..50 — would collide without unique IDs
  combined <- combine_nhanes_cycles(
    list(c1, c2),
    weight_var  = "WTDRD1",
    subject_var = "SEQN"
  )
  expect_equal(length(unique(combined$SEQN)), 100)
})

test_that("combine_nhanes_cycles with make_unique_ids=FALSE keeps original IDs", {
  c1 <- make_fake_nhanes_cycle(30, seed = 30)
  c1$SEQN <- seq(1001, 1030)
  c2 <- make_fake_nhanes_cycle(30, seed = 31)
  c2$SEQN <- seq(2001, 2030)

  combined <- combine_nhanes_cycles(
    list(c1, c2),
    weight_var     = "WTDRD1",
    subject_var    = "SEQN",
    make_unique_ids = FALSE
  )
  expect_true("1001" %in% as.character(combined$SEQN) ||
              1001    %in% combined$SEQN)
})

test_that("combine_nhanes_cycles adds cycle column", {
  c1 <- make_fake_nhanes_cycle(40, seed = 40)
  c2 <- make_fake_nhanes_cycle(40, seed = 41)
  combined <- combine_nhanes_cycles(
    list(c1, c2),
    weight_var  = "WTDRD1",
    subject_var = "SEQN"
  )
  expect_true("cycle" %in% names(combined))
  expect_setequal(unique(combined$cycle), c("cycle1", "cycle2"))
})

test_that("combine_nhanes_cycles respects custom cycle_names", {
  c1 <- make_fake_nhanes_cycle(30, seed = 50)
  c2 <- make_fake_nhanes_cycle(30, seed = 51)
  combined <- combine_nhanes_cycles(
    list(c1, c2),
    weight_var  = "WTDRD1",
    subject_var = "SEQN",
    cycle_names = c("0102", "0304")
  )
  expect_setequal(unique(combined$cycle), c("0102", "0304"))
})

test_that("combine_nhanes_cycles errors on missing weight_var", {
  c1 <- make_fake_nhanes_cycle(20, seed = 60)
  expect_error(
    combine_nhanes_cycles(list(c1), weight_var = "NONEXISTENT", subject_var = "SEQN"),
    regexp = "NONEXISTENT"
  )
})

test_that("combine_nhanes_cycles errors on empty list", {
  expect_error(
    combine_nhanes_cycles(list(), weight_var = "WTDRD1", subject_var = "SEQN"),
    regexp = "non-empty"
  )
})
