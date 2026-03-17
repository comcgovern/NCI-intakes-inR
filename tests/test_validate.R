#!/usr/bin/env Rscript
# test_nciusual.R
# -----------------------------------------------------------------
# Validate the nciusual R package against known properties of the
# NCI method using synthetic dietary recall data.
#
# Test strategy:
# 1. Generate synthetic 24-hour recall data with KNOWN usual intake
#    distribution parameters (between-person and within-person variance)
# 2. Fit the NCI model and recover those parameters
# 3. Verify that the estimated usual intake distribution matches the
#    known generating distribution
# -----------------------------------------------------------------

# Source the package files directly (for testing without install)
source("R/boxcox.R")
source("R/mixtran.R")
source("R/distrib.R")

cat(paste(rep("=", 66), collapse = ""), "\n")
cat("NCI Usual Intake R Port - Validation Tests\n")
cat(paste(rep("=", 66), collapse = ""), "\n\n")

# =====================================================================
# TEST 1: Amount-Only Model (ubiquitously consumed nutrient)
# =====================================================================
cat("TEST 1: Amount-Only Model\n")
cat("  Simulating sodium intake (consumed daily by everyone)\n")

set.seed(42)

# Known parameters
true_n_subjects <- 2000
true_lambda <- 0.30        # Box-Cox parameter
true_beta0 <- 15.0         # Intercept on transformed scale
true_beta_weekend <- 0.5   # Weekend effect
true_beta_seq <- -0.2      # Sequence (day 2) effect
true_sigma2_b <- 4.0       # Between-person variance
true_sigma2_w <- 2.5       # Within-person variance

# Generate synthetic data
# Each subject has 1-2 recalls
n_recalls_per <- sample(1:2, true_n_subjects, replace = TRUE, prob = c(0.4, 0.6))
total_recalls <- sum(n_recalls_per)

subject_ids <- rep(seq_len(true_n_subjects), times = n_recalls_per)
repeat_nums <- unlist(lapply(n_recalls_per, seq_len))

# Between-person random effects
u_i <- rep(rnorm(true_n_subjects, 0, sqrt(true_sigma2_b)),
           times = n_recalls_per)

# Covariates
weekend <- rbinom(total_recalls, 1, 0.28)  # ~28% of recalls on weekends
seq_num <- as.numeric(repeat_nums > 1)

# Survey weights (simplified)
weights <- rep(runif(true_n_subjects, 0.5, 2.0), times = n_recalls_per)

# Transformed intake
t_intake <- true_beta0 + true_beta_weekend * weekend +
  true_beta_seq * seq_num + u_i +
  rnorm(total_recalls, 0, sqrt(true_sigma2_w))

# Back-transform to original scale
intake_orig <- boxcox_inverse(t_intake, true_lambda)
intake_orig <- pmax(intake_orig, 1)  # Floor at 1 mg

sim_data <- data.frame(
  SEQN = subject_ids,
  day = repeat_nums,
  sodium_mg = intake_orig,
  weekend = weekend,
  wt = weights,
  stringsAsFactors = FALSE
)

cat(sprintf("  Generated %d recalls from %d subjects\n",
            nrow(sim_data), true_n_subjects))
cat(sprintf("  Subjects with 2 recalls: %d\n",
            sum(n_recalls_per == 2)))
cat(sprintf("  True parameters:\n"))
cat(sprintf("    lambda=%.2f, beta0=%.2f, sigma2_b=%.2f, sigma2_w=%.2f\n",
            true_lambda, true_beta0, true_sigma2_b, true_sigma2_w))
cat(sprintf("    Variance ratio (w/b): %.2f\n\n",
            true_sigma2_w / true_sigma2_b))

# Fit model
cat("  Fitting model...\n")
fit1 <- mixtran(
  data = sim_data,
  intake_var = "sodium_mg",
  subject_var = "SEQN",
  repeat_var = "day",
  model_type = "amount",
  covariates = NULL,
  weekend_var = "weekend",
  weight_var = NULL,  # Skip weights for cleaner test
  lambda = true_lambda,  # Fix lambda to test variance recovery
  verbose = TRUE
)

cat("\n  --- Parameter Recovery ---\n")
cat(sprintf("  sigma2_b: true=%.2f, estimated=%.4f\n",
            true_sigma2_b, fit1$sigma2_b))
cat(sprintf("  sigma2_w: true=%.2f, estimated=%.4f\n",
            true_sigma2_w, fit1$sigma2_w))
cat(sprintf("  var_ratio: true=%.2f, estimated=%.2f\n",
            true_sigma2_w / true_sigma2_b, fit1$var_ratio))

# Check parameter recovery (within reasonable tolerance)
b_pct_err <- abs(fit1$sigma2_b - true_sigma2_b) / true_sigma2_b * 100
w_pct_err <- abs(fit1$sigma2_w - true_sigma2_w) / true_sigma2_w * 100
cat(sprintf("  sigma2_b error: %.1f%%\n", b_pct_err))
cat(sprintf("  sigma2_w error: %.1f%%\n", w_pct_err))

if (b_pct_err < 15 && w_pct_err < 15) {
  cat("  PASS: Variance components recovered within 15%\n\n")
} else {
  cat("  WARN: Variance recovery outside 15% tolerance\n\n")
}

# Estimate usual intake distribution
cat("  Running DISTRIB...\n")
dist1 <- distrib(
  fit1,
  percentiles = c(5, 10, 25, 50, 75, 90, 95),
  n_sims = 200,
  seed = 123
)

print(dist1)

# Compare estimated vs theoretical usual intake distribution
# Theoretical: the usual intake distribution is the marginal over u_i
# T(UI) ~ N(beta0, sigma2_b), so UI ~ boxcox_inverse(N(beta0, sigma2_b))
cat("  --- Distribution Comparison ---\n")
theoretical_samples <- boxcox_inverse(
  rnorm(100000, true_beta0, sqrt(true_sigma2_b)),
  true_lambda
)
theo_mean <- mean(theoretical_samples)
theo_median <- median(theoretical_samples)

est <- dist1$results$Overall
cat(sprintf("  Mean:   theoretical=%.1f, estimated=%.1f\n",
            theo_mean, est$mean))
cat(sprintf("  Median: theoretical=%.1f, estimated=%.1f\n",
            theo_median, est$percentiles["p50"]))

mean_err <- abs(est$mean - theo_mean) / theo_mean * 100
cat(sprintf("  Mean error: %.1f%%\n", mean_err))
if (mean_err < 10) {
  cat("  PASS: Mean within 10% of theoretical\n\n")
} else {
  cat("  WARN: Mean outside 10% tolerance\n\n")
}


# =====================================================================
# TEST 2: Two-Part Uncorrelated Model (episodically consumed food)
# =====================================================================
cat("TEST 2: Two-Part Uncorrelated Model\n")
cat("  Simulating whole grain intake (episodic consumption)\n\n")

set.seed(99)

true_n <- 1500
true_alpha0 <- 0.5      # Prob intercept (logit scale) -> ~62% consume on any day
true_sigma2_v1 <- 1.0   # Between-person var (prob)
true_beta0_ep <- 3.0    # Amount intercept (transformed)
true_sigma2_v2 <- 1.5   # Between-person var (amount)
true_sigma2_e <- 1.0    # Within-person var (amount)
true_lambda_ep <- 0.25

# Generate data
n_recalls_ep <- sample(1:2, true_n, replace = TRUE, prob = c(0.35, 0.65))
total_ep <- sum(n_recalls_ep)
subj_ep <- rep(seq_len(true_n), times = n_recalls_ep)
rep_ep <- unlist(lapply(n_recalls_ep, seq_len))

# Random effects
v1 <- rep(rnorm(true_n, 0, sqrt(true_sigma2_v1)), times = n_recalls_ep)
v2 <- rep(rnorm(true_n, 0, sqrt(true_sigma2_v2)), times = n_recalls_ep)

# Probability of consumption
p_consume <- 1 / (1 + exp(-(true_alpha0 + v1)))
consumed <- rbinom(total_ep, 1, p_consume)

# Amount (conditional on consumption)
t_amount <- true_beta0_ep + v2 + rnorm(total_ep, 0, sqrt(true_sigma2_e))
amount <- boxcox_inverse(t_amount, true_lambda_ep)
amount <- pmax(amount, 0.1)

# Final intake = consumed * amount
intake_ep <- consumed * amount

sim_ep <- data.frame(
  SEQN = subj_ep,
  day = rep_ep,
  wholegrain_g = intake_ep,
  stringsAsFactors = FALSE
)

cat(sprintf("  Generated %d recalls, %.1f%% with zero intake\n",
            nrow(sim_ep), mean(intake_ep == 0) * 100))

cat("  Fitting uncorrelated two-part model...\n")
fit2 <- mixtran(
  data = sim_ep,
  intake_var = "wholegrain_g",
  subject_var = "SEQN",
  repeat_var = "day",
  model_type = "uncorr",
  lambda = true_lambda_ep,
  verbose = TRUE
)

cat("\n  --- Parameter Recovery ---\n")
cat(sprintf("  sigma2_v1 (prob):   true=%.2f, est=%.4f\n",
            true_sigma2_v1, fit2$sigma2_v1))
cat(sprintf("  sigma2_v2 (amount): true=%.2f, est=%.4f\n",
            true_sigma2_v2, fit2$sigma2_v2))
cat(sprintf("  sigma2_e  (within): true=%.2f, est=%.4f\n",
            true_sigma2_e, fit2$sigma2_e))

v1_err <- abs(fit2$sigma2_v1 - true_sigma2_v1) / true_sigma2_v1 * 100
v2_err <- abs(fit2$sigma2_v2 - true_sigma2_v2) / true_sigma2_v2 * 100
e_err <- abs(fit2$sigma2_e - true_sigma2_e) / true_sigma2_e * 100
cat(sprintf("  Errors: v1=%.1f%%, v2=%.1f%%, e=%.1f%%\n", v1_err, v2_err, e_err))

if (max(v1_err, v2_err, e_err) < 25) {
  cat("  PASS: All variance components within 25%\n\n")
} else {
  cat("  WARN: Some variance components outside 25% tolerance\n\n")
}

# Estimate distribution
cat("  Running DISTRIB for two-part model...\n")
dist2 <- distrib(fit2, n_sims = 200, seed = 456)
print(dist2)


# =====================================================================
# TEST 3: Box-Cox Lambda Selection
# =====================================================================
cat("TEST 3: Box-Cox Lambda Optimization\n")

set.seed(77)

# Generate data with known lambda = 0.5 (square root transform)
y_test <- (rnorm(5000, mean = 5, sd = 2))^2  # chi-sq-like, best lambda ~0.5
y_test <- y_test[y_test > 0]

bc_result <- find_optimal_lambda(y_test)
cat(sprintf("  True lambda ≈ 0.5, Estimated: %.3f (raw: %.3f)\n",
            bc_result$lambda, bc_result$lambda_raw))
cat(sprintf("  Used log: %s\n", bc_result$used_log))

lambda_err <- abs(bc_result$lambda_raw - 0.5)
if (lambda_err < 0.15) {
  cat("  PASS: Lambda within 0.15 of true value\n\n")
} else {
  cat("  WARN: Lambda error > 0.15\n\n")
}


# =====================================================================
# TEST 4: Verify within-person variance removal
# =====================================================================
cat("TEST 4: Within-Person Variance Removal\n")
cat("  The estimated usual intake distribution should be NARROWER than\n")
cat("  the observed single-day distribution.\n\n")

obs_sd <- sd(sim_data$sodium_mg)
usual_sd <- dist1$results$Overall$sd

cat(sprintf("  Observed (single-day) SD: %.1f\n", obs_sd))
cat(sprintf("  Estimated usual intake SD: %.1f\n", usual_sd))
cat(sprintf("  Ratio (usual/observed):    %.3f\n", usual_sd / obs_sd))

if (usual_sd < obs_sd) {
  cat("  PASS: Usual intake distribution is narrower (as expected)\n\n")
} else {
  cat("  FAIL: Usual intake distribution should be narrower!\n\n")
}


# =====================================================================
# SUMMARY
# =====================================================================
cat(paste(rep("=", 66), collapse = ""), "\n")
cat("Validation complete. Review results above for any warnings.\n")
cat("Key checks:\n")
cat("  1. Variance component recovery (amount model)\n")
cat("  2. Variance component recovery (two-part model)\n")
cat("  3. Box-Cox lambda optimization\n")
cat("  4. Within-person variance removal (distribution narrowing)\n")
cat(paste(rep("=", 66), collapse = ""), "\n")
