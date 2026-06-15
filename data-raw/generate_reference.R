#!/usr/bin/env Rscript
# =====================================================================
# generate_reference.R
# ---------------------------------------------------------------------
# Generate ANALYTIC GROUND-TRUTH reference outputs for parity testing.
#
# We cannot run the proprietary NCI SAS macros (MIXTRAN / DISTRIB) in
# this environment, so instead of comparing R against another estimator
# we compare R against the *truth*: synthetic 24-hour recall datasets are
# generated from KNOWN parameters, and the true usual-intake distribution
# implied by those parameters is computed independently by large-sample
# Monte Carlo. A correct SAS run would be estimating exactly this same
# target distribution, so recovering it is the meaningful parity check.
#
# Outputs (all committed to the repo so the reference values are visible
# and reviewable, and the datasets are frozen across R/RNG versions):
#   tests/testthat/fixtures/amount_recall_data.csv
#   tests/testthat/fixtures/uncorr_recall_data.csv
#   tests/testthat/fixtures/corr_recall_data.csv
#   tests/testthat/fixtures/reference_params.csv     (true generating params)
#   tests/testthat/fixtures/reference_outputs.csv    (true distribution stats)
#
# Run with:  Rscript data-raw/generate_reference.R
# =====================================================================

suppressMessages(library(nciusual))  # for boxcox_inverse()

FIX <- "tests/testthat/fixtures"
PERCENTILES <- c(5, 10, 25, 50, 75, 90, 95)
REF_MC_N    <- 4e6   # Monte Carlo size for the analytic ground-truth target

# Percentiles use base R quantile type = 5 (Hazen / midpoint), which matches
# the package's weighted_quantiles() convention for equal weights.
ref_stats <- function(ui) {
  ui <- ui[is.finite(ui)]
  qs <- as.numeric(stats::quantile(ui, PERCENTILES / 100, type = 5))
  out <- c(mean = mean(ui), sd = stats::sd(ui),
           setNames(qs, paste0("p", PERCENTILES)))
  out
}

reference_rows <- list()
param_rows     <- list()

add_reference <- function(model, ui) {
  s <- ref_stats(ui)
  reference_rows[[model]] <<- data.frame(
    model = model, statistic = names(s), value = as.numeric(s),
    stringsAsFactors = FALSE
  )
}

# =====================================================================
# 1. AMOUNT MODEL
#    T(Y) = beta0 + beta_we * weekend + u + e ,  u ~ N(0, s2_b), e ~ N(0, s2_w)
#    Usual intake (weekend removed) = boxcox_inverse( N(beta0, s2_b) )
# =====================================================================
set.seed(1001)
p_amt <- list(lambda = 0.30, beta0 = 15.0, beta_weekend = 0.5,
              sigma2_b = 4.0, sigma2_w = 2.5, n_subjects = 4000)
param_rows$amount <- data.frame(model = "amount",
  parameter = names(p_amt), value = unlist(p_amt), stringsAsFactors = FALSE)

nr   <- sample(1:2, p_amt$n_subjects, replace = TRUE, prob = c(0.3, 0.7))
sid  <- rep(seq_len(p_amt$n_subjects), nr)
day  <- unlist(lapply(nr, seq_len))
u_i  <- rep(rnorm(p_amt$n_subjects, 0, sqrt(p_amt$sigma2_b)), nr)
wknd <- rbinom(length(sid), 1, 0.28)
t_y  <- p_amt$beta0 + p_amt$beta_weekend * wknd + u_i +
        rnorm(length(sid), 0, sqrt(p_amt$sigma2_w))
y    <- pmax(boxcox_inverse(t_y, p_amt$lambda), 0.01)
amount_data <- data.frame(SEQN = sid, day = day, intake = y, weekend = wknd)
write.csv(amount_data, file.path(FIX, "amount_recall_data.csv"), row.names = FALSE)

# Analytic ground truth: marginal usual intake distribution
set.seed(20001)
ui_true <- pmax(boxcox_inverse(
  rnorm(REF_MC_N, p_amt$beta0, sqrt(p_amt$sigma2_b)), p_amt$lambda), 0)
add_reference("amount", ui_true)

# =====================================================================
# 2. TWO-PART UNCORRELATED MODEL
#    P(consume) = plogis(alpha0 + v1),  v1 ~ N(0, s2_v1)
#    T(amount)  = beta0 + v2 + e,       v2 ~ N(0, s2_v2),  e ~ N(0, s2_e)
#    Usual intake = plogis(alpha0 + v1) * boxcox_inverse(beta0 + v2)
# =====================================================================
set.seed(1002)
p_unc <- list(lambda = 0.25, alpha0 = 0.3, sigma2_v1 = 1.0,
              beta0 = 3.0, sigma2_v2 = 1.0, sigma2_e = 0.8,
              rho = 0.0, n_subjects = 4000)
param_rows$uncorr <- data.frame(model = "uncorr",
  parameter = names(p_unc), value = unlist(p_unc), stringsAsFactors = FALSE)

nr  <- sample(1:2, p_unc$n_subjects, replace = TRUE, prob = c(0.35, 0.65))
sid <- rep(seq_len(p_unc$n_subjects), nr)
day <- unlist(lapply(nr, seq_len))
v1s <- rnorm(p_unc$n_subjects, 0, sqrt(p_unc$sigma2_v1))
v2s <- rnorm(p_unc$n_subjects, 0, sqrt(p_unc$sigma2_v2))
v1  <- rep(v1s, nr); v2 <- rep(v2s, nr)
consumed <- rbinom(length(sid), 1, plogis(p_unc$alpha0 + v1))
amt <- pmax(boxcox_inverse(p_unc$beta0 + v2 +
              rnorm(length(sid), 0, sqrt(p_unc$sigma2_e)), p_unc$lambda), 0.01)
uncorr_data <- data.frame(SEQN = sid, day = day, intake = consumed * amt)
write.csv(uncorr_data, file.path(FIX, "uncorr_recall_data.csv"), row.names = FALSE)

set.seed(20002)
v1m <- rnorm(REF_MC_N, 0, sqrt(p_unc$sigma2_v1))
v2m <- rnorm(REF_MC_N, 0, sqrt(p_unc$sigma2_v2))
ui_true <- plogis(p_unc$alpha0 + v1m) *
           pmax(boxcox_inverse(p_unc$beta0 + v2m, p_unc$lambda), 0)
add_reference("uncorr", ui_true)

# =====================================================================
# 3. TWO-PART CORRELATED MODEL (rho = 0.5)
# =====================================================================
set.seed(1003)
p_cor <- list(lambda = 0.25, alpha0 = 0.3, sigma2_v1 = 1.2,
              beta0 = 3.0, sigma2_v2 = 1.0, sigma2_e = 0.8,
              rho = 0.5, n_subjects = 4000)
param_rows$corr <- data.frame(model = "corr",
  parameter = names(p_cor), value = unlist(p_cor), stringsAsFactors = FALSE)

nr  <- sample(1:2, p_cor$n_subjects, replace = TRUE, prob = c(0.35, 0.65))
sid <- rep(seq_len(p_cor$n_subjects), nr)
day <- unlist(lapply(nr, seq_len))
sd1 <- sqrt(p_cor$sigma2_v1); sd2 <- sqrt(p_cor$sigma2_v2)
z1  <- rnorm(p_cor$n_subjects); z2 <- rnorm(p_cor$n_subjects)
v1s <- sd1 * z1
v2s <- sd2 * (p_cor$rho * z1 + sqrt(1 - p_cor$rho^2) * z2)
v1  <- rep(v1s, nr); v2 <- rep(v2s, nr)
consumed <- rbinom(length(sid), 1, plogis(p_cor$alpha0 + v1))
amt <- pmax(boxcox_inverse(p_cor$beta0 + v2 +
              rnorm(length(sid), 0, sqrt(p_cor$sigma2_e)), p_cor$lambda), 0.01)
corr_data <- data.frame(SEQN = sid, day = day, intake = consumed * amt)
write.csv(corr_data, file.path(FIX, "corr_recall_data.csv"), row.names = FALSE)

set.seed(20003)
z1m <- rnorm(REF_MC_N); z2m <- rnorm(REF_MC_N)
v1m <- sd1 * z1m
v2m <- sd2 * (p_cor$rho * z1m + sqrt(1 - p_cor$rho^2) * z2m)
ui_true <- plogis(p_cor$alpha0 + v1m) *
           pmax(boxcox_inverse(p_cor$beta0 + v2m, p_cor$lambda), 0)
add_reference("corr", ui_true)

# =====================================================================
# Write reference tables
# =====================================================================
reference_outputs <- do.call(rbind, reference_rows)
reference_params  <- do.call(rbind, param_rows)
rownames(reference_outputs) <- NULL
rownames(reference_params)  <- NULL

write.csv(reference_outputs, file.path(FIX, "reference_outputs.csv"), row.names = FALSE)
write.csv(reference_params,  file.path(FIX, "reference_params.csv"),  row.names = FALSE)

cat("Wrote reference fixtures to", FIX, "\n")
print(reference_outputs)
