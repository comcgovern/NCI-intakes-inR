# `nciusual`: R Port of NCI Usual Dietary Intake Estimation
## Package Specification v1.0

---

## 1. Purpose

This package provides an R-native implementation of the National Cancer Institute (NCI) method for estimating usual dietary intake distributions from short-term dietary assessment data (24-hour recalls). The method is currently implemented exclusively in SAS macros (MIXTRAN v2.21, DISTRIB v2.2, INDIVINT) distributed by NCI at https://epi.grants.cancer.gov/diet/usualintakes/macros.html.

No maintained R implementation exists. This package fills that gap.

**Primary use case:** Estimating population-level usual intake distributions of nutrients and foods from NHANES 24-hour dietary recall data, with proper handling of complex survey designs, within-person measurement error, and episodic consumption.

**Target users:** Researchers analyzing NHANES dietary data in R who currently must either run SAS or use simplified survey-weighted means.

---

## 2. Statistical Background

### 2.1 The Core Problem

A single 24-hour dietary recall captures one day of eating, but dietary recommendations and diet-disease relationships are based on *usual* (long-run average) intake. The distribution of single-day intakes is wider than the distribution of usual intakes because it conflates two sources of variation:

- **Between-person variation** (σ²_b): Real differences in habitual diets across individuals
- **Within-person variation** (σ²_w): Day-to-day fluctuation in what a single person eats

The NCI method uses repeated 24-hour recalls (≥2 per person, in at least a subset) to estimate and remove the within-person component, recovering the narrower distribution of usual intakes.

### 2.2 Model Specifications

The NCI method uses three model types depending on the consumption pattern of the dietary component:

#### 2.2.1 Amount-Only Model

For nutrients/foods consumed daily by nearly everyone (energy, sodium, added sugars, total protein).

**Model:**

```
T(Y_ij) = X_ij'β + u_i + ε_ij
```

Where:
- Y_ij = observed intake for person i on recall j
- T(·) = Box-Cox transformation with parameter λ
- X_ij = covariate vector (intercept, sequence indicator, weekend indicator, demographics)
- β = fixed effect coefficients
- u_i ~ N(0, σ²_b) = between-person random effect
- ε_ij ~ N(0, σ²_w) = within-person error
- u_i ⊥ ε_ij

**Usual intake for person i:**

```
T(UI_i) = X_i'β + u_i    (no ε term — that's the whole point)
UI_i = T⁻¹(X_i'β + u_i)
```

The *population distribution* of usual intake is obtained by integrating over u_i.

#### 2.2.2 Two-Part Model (Episodically Consumed Foods)

For foods consumed on some days but not others (whole grains, fish, specific fruits). The model separates the probability of consumption from the amount consumed:

**Part 1 — Probability of consumption:**

```
logit(P(Y_ij > 0)) = Z_ij'α + v1_i
```

Where:
- Z_ij = covariate vector for probability model
- α = fixed effect coefficients (probability)
- v1_i ~ N(0, σ²_v1) = between-person random effect for consumption probability

**Part 2 — Amount consumed on consumption days:**

```
T(Y_ij | Y_ij > 0) = X_ij'β + v2_i + ε_ij
```

Where:
- v2_i ~ N(0, σ²_v2) = between-person random effect for amount
- ε_ij ~ N(0, σ²_e) = within-person error in amount

**Correlation structure:**

- **Uncorrelated model:** Cov(v1_i, v2_i) = 0. Fit each part independently.
- **Correlated model:** (v1_i, v2_i) ~ BVN(0, Σ), where Σ has off-diagonal ρ·σ_v1·σ_v2. Requires joint estimation.

**Usual intake for person i:**

```
UI_i = P_i × A_i
     = logistic(Z_i'α + v1_i) × T⁻¹(X_i'β + v2_i)
```

### 2.3 Box-Cox Transformation

The NCI method applies a Box-Cox transformation to consumption-day amounts to approximate normality on the transformed scale:

```
T(y; λ) = (y^λ - 1) / λ     if λ ≠ 0
         = log(y)              if λ = 0
```

Inverse:

```
T⁻¹(z; λ) = (λz + 1)^(1/λ)  if λ ≠ 0
           = exp(z)            if λ = 0
```

**Lambda selection:** Profile likelihood over a grid (typically 0.01 to 1.00 by 0.01). NCI convention: if the optimal λ < 0.15, use log transform (λ = 0) because the method may perform poorly with very small λ values (documented in simulation studies by Dodd et al.).

### 2.4 DISTRIB: Monte Carlo Estimation of Usual Intake Distribution

Given the fitted model parameters, the distribution of usual intake is estimated via simulation:

**Amount-only model:**
1. For each of `n_sims` Monte Carlo repetitions r = 1, ..., R:
   - For each person i in the sample:
     - Draw u_i^(r) ~ N(0, σ̂²_b)
     - Compute T(UI_i^(r)) = X_i'β̂ + u_i^(r)
     - Back-transform: UI_i^(r) = T⁻¹(T(UI_i^(r)))
2. Compute weighted percentiles, means, and cutpoint proportions from the simulated usual intakes, using survey weights.

**Two-part model:**
1. For each repetition r, for each person i:
   - Draw (v1_i^(r), v2_i^(r)) from BVN(0, Σ̂) (or independently if uncorrelated)
   - Compute P_i^(r) = logistic(Z_i'α̂ + v1_i^(r))
   - Compute A_i^(r) = T⁻¹(X_i'β̂ + v2_i^(r))
   - Compute UI_i^(r) = P_i^(r) × A_i^(r)
2. Compute weighted distribution statistics as above.

**Key property:** Within-person error ε is NOT drawn in the simulation. This is what makes the estimated distribution narrower than the observed distribution.

### 2.5 Standard Errors via Balanced Repeated Replication (BRR)

For NHANES complex survey data, standard errors for percentiles and means of the usual intake distribution are obtained via BRR:

1. Generate K sets of BRR replicate weights (NHANES provides these, or they can be computed; typically K = 32–72 for multi-cycle analyses)
2. Run the full MIXTRAN → DISTRIB pipeline K times, each time using a different replicate weight
3. Compute the BRR variance of each statistic:
   ```
   Var_BRR(θ̂) = (1/K) × Σ_k (θ̂_k - θ̂_0)²
   ```
   where θ̂_0 is the full-sample estimate and θ̂_k is the k-th replicate estimate
4. With Fay's adjustment factor F (NHANES standard: F = 0.3):
   ```
   Var_BRR(θ̂) = 1/(K × (1-F)²) × Σ_k (θ̂_k - θ̂_0)²
   ```

---

## 3. Package API Design

### 3.1 Core Functions (Public)

```r
# --- Model Fitting (MIXTRAN equivalent) ---

mixtran(
  data,                    # data.frame, long format (1 row per recall per person)
  intake_var,              # character: column name of dietary intake
  subject_var,             # character: column name of subject ID
  repeat_var,              # character: column name of recall sequence (1, 2, ...)
  model_type = "amount",   # "amount" | "uncorr" | "corr"
  covariates = NULL,       # character vector: covariate column names
  weekend_var = NULL,      # character: weekend indicator column name
  weight_var = NULL,       # character: survey weight column name
  lambda = NULL,           # numeric or NULL (auto-estimate)
  lambda_grid = seq(0.01, 1.0, by = 0.01),
  min_positive = NULL,     # numeric: replacement for zero values
  verbose = TRUE
) -> mixtran_fit object


# --- Distribution Estimation (DISTRIB equivalent) ---

distrib(
  mixtran_obj,             # mixtran_fit from mixtran()
  subgroup_var = NULL,     # character: optional subgroup variable name
  cutpoints = NULL,        # numeric vector: intake cutpoints for % above/below
  percentiles = c(5, 10, 25, 50, 75, 90, 95),
  n_sims = 100,            # integer: Monte Carlo repetitions per person
  seed = 12345
) -> distrib_result object


# --- Individual Predictions (INDIVINT equivalent) ---

indivint(
  mixtran_obj,             # mixtran_fit from mixtran()
  newdata = NULL           # optional: predict for new data
) -> data.frame with predicted usual intakes


# --- BRR Standard Errors ---

brr_usual_intake(
  data,                    # full dataset
  replicate_weights,       # matrix: K columns of BRR replicate weights
  mixtran_args,            # list: arguments to pass to mixtran()
  distrib_args,            # list: arguments to pass to distrib()
  fay_factor = 0.3         # Fay's adjustment factor
) -> brr_result object with SEs for all statistics


# --- Convenience Wrappers ---

nci_usual(
  data,                    # shortcut: runs mixtran() then distrib()
  intake_var,
  subject_var,
  repeat_var,
  ...                      # passed to mixtran() and distrib()
) -> distrib_result object


# --- Box-Cox Utilities ---

boxcox_transform(y, lambda) -> numeric vector
boxcox_inverse(z, lambda) -> numeric vector
find_optimal_lambda(y, weights = NULL, lambda_grid, covariates = NULL) -> list
```

### 3.2 Return Object Structures

#### `mixtran_fit`

```r
# Common fields (all model types):
$model_type        # "amount" | "uncorr" | "corr"
$lambda            # Box-Cox parameter used
$beta              # named numeric: fixed effects for amount model
$intake_var        # character
$subject_var       # character
$covariates        # character vector
$weekend_var       # character or NULL
$n_subjects        # integer
$n_recalls         # integer
$n_positive        # integer (recalls with intake > 0)
$pct_consuming     # numeric (% of recalls with consumption)
$predicted         # data.frame: subject-level linear predictors + weights
$converged         # logical

# Amount-only specific:
$sigma2_b          # between-person variance
$sigma2_w          # within-person variance
$var_ratio         # sigma2_w / sigma2_b
$model_fit         # nlme::lme object (for diagnostics)

# Two-part specific (uncorr and corr):
$alpha             # named numeric: fixed effects for probability model
$sigma2_v1         # between-person variance (probability)
$sigma2_v2         # between-person variance (amount)
$sigma2_e          # within-person variance (amount)
$rho               # correlation between v1 and v2 (0 for uncorr)
```

#### `distrib_result`

```r
$results           # named list, one entry per subgroup + "Overall"
  $<subgroup>
    $mean          # weighted mean usual intake
    $sd            # weighted SD of usual intake
    $percentiles   # named numeric vector (p5, p10, p25, p50, p75, p90, p95)
    $cutpoint_below  # named numeric: proportion below each cutpoint
    $cutpoint_above  # named numeric: proportion above each cutpoint
    $n_simulated     # total simulated values
$simulated         # list with raw simulated usual intake values
$mixtran_obj       # the mixtran_fit used
$n_sims            # MC repetitions
$seed              # random seed
```

### 3.3 Covariate Conventions

Matching SAS MIXTRAN behavior:

- **Sequence indicator:** Automatically created from `repeat_var`. First recall = 0, subsequent = 1. Named `seq_num` internally. Effect removed from usual intake predictions in DISTRIB.
- **Weekend indicator:** User-supplied binary variable via `weekend_var`. Effect removed from usual intake predictions in DISTRIB.
- **User covariates:** Must be binary (0/1) or continuous. Categorical variables with >2 levels must be dummy-coded by the user before calling `mixtran()`. These effects are RETAINED in usual intake predictions (they represent real population structure).
- **Subgroup variables:** NOT included as model covariates. Used only in DISTRIB to stratify output. The `subgroup_var` in `distrib()` must be present in the predicted data.

**Effect removal in DISTRIB:** When computing usual intake, sequence and weekend effects are removed by setting `seq_num = 0` and `weekend = 0` (reference levels). This means the usual intake estimate represents intake on a weekday, first-recall basis — the NCI convention. User covariates (e.g., age, sex, race) are kept at their observed values.

---

## 4. Implementation Strategy

### 4.1 Amount-Only Model

**R implementation:** `nlme::lme()` with formula-based interface.

```r
# Model: T(intake) ~ covariates, random = ~ 1 | subject
nlme::lme(
  fixed = t_intake ~ seq_num + weekend + <user covariates>,
  random = ~ 1 | subject,
  data = work,
  method = "REML"
)
```

**Variance extraction:**
```r
vc <- nlme::VarCorr(fit)
sigma2_b <- as.numeric(vc["(Intercept)", "Variance"])
sigma2_w <- as.numeric(vc["Residual", "Variance"])
```

**Survey weights:** The SAS MIXTRAN uses `PROC NLMIXED` with BRR replicate weights (not direct incorporation of weights into the likelihood). In R, this means:
- The base model fit does NOT use survey weights directly
- Survey weights are used in DISTRIB for weighted percentile/mean computation
- BRR replicate weights drive standard error estimation

This is a key design decision: `nlme::lme` does support a `weights` argument (`varFixed`), but NCI's approach is to use BRR rather than pseudo-likelihood weighting for the NLMIXED step. Our implementation should follow the NCI convention: **unweighted model fitting, weighted distribution estimation, BRR for SEs.**

**Prototype validation result:** σ²_b recovered within 5.3%, σ²_w within 1.0% of true values. Mean usual intake within 1.4% of theoretical. Distribution correctly narrower than observed (SD ratio 0.79).

### 4.2 Two-Part Uncorrelated Model

**Part 1 (probability):** Logistic GLMM with random intercept.

```r
MASS::glmmPQL(
  fixed = consumed ~ seq_num + weekend + <covariates>,
  random = ~ 1 | subject,
  family = binomial(link = "logit"),
  data = work
)
```

`MASS::glmmPQL` uses penalized quasi-likelihood, which is faster and more stable than full ML for binary GLMMs, especially with unbalanced data. SAS MIXTRAN uses `PROC NLMIXED` which does full ML via adaptive Gaussian quadrature — this is a known approximation difference. For the purposes of usual intake estimation (which primarily needs the variance component, not precise fixed effects), PQL is an accepted approximation.

**Alternative:** `lme4::glmer()` with Laplace or adaptive GH quadrature provides full ML. Consider offering this as an option for users who want closer SAS parity, accepting slower runtime.

**Part 2 (amount):** Same as amount-only model, but fit only on rows where `intake > 0`.

```r
nlme::lme(
  fixed = t_intake ~ seq_num + weekend + <covariates>,
  random = ~ 1 | subject,
  data = work_positive_only
)
```

**Prototype validation result:** All three variance components recovered within 6% of true values.

### 4.3 Two-Part Correlated Model

This is the hardest piece. The SAS implementation uses `PROC NLMIXED` to jointly maximize the marginal likelihood with adaptive Gaussian quadrature over the bivariate random effects. R has no direct equivalent for a joint binary/continuous mixed model with correlated random effects.

#### 4.3.1 Approach A: Profile Likelihood over ρ (Recommended Default)

1. Fit the uncorrelated model to get starting values for (α, β, σ²_v1, σ²_v2, σ²_e)
2. For each ρ on a grid (e.g., -0.9 to 0.9 by 0.1):
   - Transform the data using the Cholesky decomposition implied by ρ
   - Re-fit the two parts with adjusted residuals
   - Compute the approximate profile log-likelihood
3. Select ρ that maximizes the profile likelihood
4. Re-fit final model at the selected ρ

**Advantages:** Uses existing fast `nlme`/`glmmPQL` machinery. No custom optimization. Stable.
**Disadvantages:** Approximate. Grid resolution limits precision of ρ̂ (can refine with golden section search). Does not jointly account for uncertainty in ρ and other parameters.

**Implementation note:** The profile approach can be refined iteratively — start with a coarse grid, narrow around the maximum, re-estimate on a fine grid.

#### 4.3.2 Approach B: Gauss-Hermite Quadrature with Direct Optimization

Joint marginal log-likelihood:

```
ℓ(θ) = Σ_i log ∫∫ L_i(Y_i | v1, v2; α, β, σ²_e) × φ₂(v1, v2; Σ_v) dv1 dv2
```

where:
- L_i is the product over person i's recalls of [P(consume)^d × P(not consume)^(1-d) × f(T(Y)|consume)]
- φ₂ is the bivariate normal density with covariance Σ_v

The bivariate integral is approximated by Gauss-Hermite quadrature:

```
∫∫ g(v1, v2) φ₂(v1, v2) dv1 dv2 ≈ Σ_q w_q × g(z1_q, z2_q)
```

where nodes/weights come from the tensor product of 1D GH rules, transformed from standard normal to the correlated structure via Cholesky decomposition.

Optimize using `stats::optim()` with BFGS or L-BFGS-B.

**Critical performance issue:** With N subjects and Q quadrature nodes, each likelihood evaluation is O(N × Q² × J_max) where J_max is max recalls per person. In the prototype, n=200 with 10×10 GH nodes timed out at 120 seconds. This needs to be addressed.

**Performance solutions (implement in priority order):**

1. **Vectorize the inner loops.** The prototype iterated over subjects and quadrature nodes in nested R `for` loops. Restructure so the likelihood contribution for all quadrature nodes is computed as matrix operations per subject, and the subject loop is vectorized where possible.

2. **Reduce quadrature nodes.** SAS NLMIXED uses adaptive GH, which places nodes where the integrand is large. Start with n=5 (25 bivariate nodes) rather than n=10 (100 nodes). For smooth integrands this is often sufficient.

3. **Pre-compute constant terms.** Box-Cox transform, design matrix products, and Jacobian terms can all be computed once outside the optimization loop.

4. **C++ inner loop via Rcpp.** If pure-R vectorization is insufficient, write the per-subject log-likelihood contribution in C++ and call via Rcpp. This is the most impactful optimization — expect 50-100x speedup for the inner loop.

5. **Laplace approximation.** Replace GH quadrature with a Laplace approximation (mode + curvature of the integrand). This is what `lme4::glmer` uses by default and is typically sufficient for variance component estimation. Requires computing the mode of the random effects for each subject, which is itself an optimization problem but a small one (2 parameters per subject).

#### 4.3.3 Approach C: Template Model Builder (TMB)

TMB (R package `TMB`) provides automatic differentiation and Laplace approximation for marginal likelihoods with random effects, compiled to C++. This is architecturally the closest R equivalent to what SAS PROC NLMIXED does.

Write the joint model in TMB's C++ template language:
- Define the negative log-likelihood as a function of (α, β, log_σ²_v1, log_σ²_v2, log_σ²_e, atanh_ρ) as fixed parameters and (v1_i, v2_i) as random effects
- TMB handles the Laplace approximation to integrate out random effects
- Derivatives computed via automatic differentiation → fast optimization

**Advantages:** Fast, accurate, principled. Closest to SAS behavior.
**Disadvantages:** Adds TMB as a dependency (requires C++ compilation). More complex to maintain.

#### 4.3.4 Recommendation

- **Phase 1 (release v0.1):** Implement Approach A (profile likelihood over ρ) as default. Fast, robust, uses existing tooling. Document the approximation.
- **Phase 2 (v0.2):** Implement Approach B with Rcpp-accelerated inner loop as `engine = "ghq"` option.
- **Phase 3 (v0.3):** Implement Approach C with TMB as `engine = "tmb"` option for full SAS parity.
- Users who only need uncorrelated models (many practical applications) are already served by v0.1.

### 4.4 DISTRIB Implementation

The DISTRIB step is straightforward Monte Carlo simulation. Key implementation details:

**Per-person simulation (amount-only):**
```r
# For each person i (vectorized over all persons):
u_sim <- matrix(rnorm(n_subjects * n_sims, 0, sqrt(sigma2_b)),
                nrow = n_subjects, ncol = n_sims)
t_usual <- linpred + u_sim   # linpred is the covariate-adjusted prediction
usual_intake <- boxcox_inverse(t_usual, lambda)
```

**Effect removal:** The linear predictor used in DISTRIB should be computed at reference covariate values for sequence and weekend (both = 0). User covariates remain at observed values. In practice: from the MIXTRAN predicted values, subtract β̂_seq × seq_num and β̂_weekend × weekend.

**Weighted percentiles:** Use survey weights for all summary statistics. The weighted quantile calculation uses the standard interpolation approach: sort values, compute cumulative weight proportions, interpolate.

**Two-part simulation:**
```r
# Draw correlated random effects
z1 <- matrix(rnorm(n_subjects * n_sims), n_subjects, n_sims)
z2 <- matrix(rnorm(n_subjects * n_sims), n_subjects, n_sims)
v1 <- sqrt(sigma2_v1) * z1
v2 <- sqrt(sigma2_v2) * (rho * z1 + sqrt(1 - rho^2) * z2)

# Probability and amount
prob <- plogis(prob_linpred + v1)
amount <- boxcox_inverse(amt_linpred + v2, lambda)

# Usual intake
usual_intake <- prob * amount
```

**Prototype validation result:** Mean within 1.4% of theoretical. Distribution correctly narrower than observed.

### 4.5 BRR Standard Errors

This is architecturally a wrapper that:
1. Accepts BRR replicate weights (matrix, n_subjects × K)
2. Runs MIXTRAN + DISTRIB K times, each with a different weight column
3. Collects the target statistics (percentiles, means, cutpoints) from each replicate
4. Computes BRR variance using Fay's formula

**NHANES BRR weight generation** (for multi-cycle analyses):
- NHANES provides masked variance units (SDMVSTRA, SDMVPSU)
- BRR weights can be generated using the `survey` package or manually using the Hadamard matrix approach documented by Herrick et al. (2018)
- Perturbation factor F = 0.3 is the NHANES standard (Fay's modification)
- Number of replicates: typically 2 × (number of strata with 2 PSUs); NHANES 2-year cycle → ~32; multi-cycle → ~72

**Parallelization:** BRR replicates are independent and embarrassingly parallel. Support parallel execution via `future.apply::future_lapply()` or `parallel::mclapply()`.

```r
brr_results <- future.apply::future_lapply(1:K, function(k) {
  # Run mixtran + distrib with replicate weight k
  fit_k <- mixtran(data, ..., weight_var = brr_weight_names[k])
  dist_k <- distrib(fit_k, ...)
  distrib_to_df(dist_k)
})
```

### 4.6 INDIVINT (Individual Predicted Intake)

The INDIVINT macro predicts individual usual intake using the Best Linear Unbiased Predictor (BLUP) from the mixed model. This is primarily used for regression calibration (using usual intake as a predictor of a health outcome).

**Amount-only:** Predicted usual intake = X_i'β̂ + û_i, where û_i is the estimated random effect (shrinkage estimator from the mixed model).

```r
# From nlme::lme fit:
predicted_usual <- predict(fit, level = 1)  # includes random effect
# Remove sequence and weekend effects
predicted_usual <- predicted_usual - beta_seq * seq_num - beta_weekend * weekend
# Back-transform
predicted_usual_orig <- boxcox_inverse(predicted_usual, lambda)
```

**Two-part:** More complex. Requires empirical Bayes estimates of (v1_i, v2_i) and integration. See Kipnis et al. (2009) for the adaptive Gaussian quadrature approach used in the SAS INDIVINT macro.

**Implementation priority:** Lower than MIXTRAN and DISTRIB. Include in v0.2.

---

## 5. Data Preparation Conventions

### 5.1 Input Data Format

The primary input is a long-format data frame where each row represents one 24-hour dietary recall for one person. Required columns:

| Column | Type | Description |
|--------|------|-------------|
| Subject ID | character/integer | Unique person identifier |
| Repeat number | integer | Recall sequence (1 = first, 2 = second) |
| Intake | numeric | Dietary intake of the target component (≥ 0) |

Optional columns:

| Column | Type | Description |
|--------|------|-------------|
| Weekend | binary (0/1) | 1 if recall was collected Fri–Sun |
| Survey weight | numeric | Sampling weight |
| Covariates | binary/continuous | Demographic or other predictors |

### 5.2 NHANES-Specific Helper

Provide a convenience function for NHANES data preparation:

```r
prepare_nhanes(
  individual_food_file,    # DR1IFF / DR2IFF merged
  demo_file,               # DEMO
  dietary_component,       # e.g., "DR1ISODI" for sodium
  food_source_filter = 1,  # 1 = school cafeteria (for school meal analyses)
  cycles = c("2015-2016", "2017-2020", "2021-2023"),
  age_range = c(5, 18)     # children only
) -> data.frame ready for mixtran()
```

This helper should:
- Download data via `nhanesA` package
- Merge Day 1 and Day 2 recalls into long format
- Create the weekend indicator from recall day-of-week (NHANES variable `DR1DAY`)
- Create the sequence indicator
- Compute appropriate multi-cycle survey weights
- Filter to target population and food source
- Handle the NHANES 2017-2020 pre-pandemic subsetting

### 5.3 Zero Handling

- **Amount-only model:** Zeros in ubiquitous nutrients are treated as reporting artifacts. Replace with half the minimum positive observed value. (NCI convention; documented in MIXTRAN v2.1 release notes.)
- **Two-part models:** Zeros are genuine non-consumption days. Handled by Part 1 (probability model).

---

## 6. SAS Macro Parity Matrix

| SAS Feature | R Implementation | Parity Status |
|------------|-----------------|---------------|
| MIXTRAN: amount-only model | `nlme::lme` | Full (validated) |
| MIXTRAN: uncorrelated two-part | `MASS::glmmPQL` + `nlme::lme` | Approximate (PQL vs. full ML) |
| MIXTRAN: correlated two-part | Profile likelihood (v0.1); GHQ/TMB (later) | Approximate initially |
| MIXTRAN: Box-Cox lambda search | Profile likelihood on grid | Full |
| MIXTRAN: covariate effects | Formula-based | Full |
| MIXTRAN: saved parameter datasets | `mixtran_fit` object | Full |
| MIXTRAN: saved predicted datasets | `$predicted` data.frame | Full |
| DISTRIB: Monte Carlo simulation | Vectorized R | Full |
| DISTRIB: weighted percentiles | Custom weighted quantile function | Full |
| DISTRIB: cutpoint proportions | Weighted proportion computation | Full |
| DISTRIB: subgroup estimation | `subgroup_var` argument | Full |
| INDIVINT: individual BLUP predictions | Planned v0.2 | Not yet |
| BRR standard errors | Parallel replicate pipeline | Planned v0.1 |
| Complex survey weight handling | `survey` package integration | Planned v0.1 |
| Starting values from prior run | `start_params` argument | Planned v0.1 |
| Convergence diagnostics | Standard R model diagnostics | Partial |

### Known Parity Gaps

1. **PQL vs. full ML for logistic GLMM.** SAS MIXTRAN uses full ML (PROC NLMIXED with adaptive GH). R's `MASS::glmmPQL` uses penalized quasi-likelihood, which is known to attenuate variance component estimates for binary outcomes, especially with low prevalence. Mitigation: offer `lme4::glmer` as alternative engine.

2. **Survey weight integration in model fitting.** SAS MIXTRAN feeds survey weights to PROC NLMIXED via a replicate statement. The NCI approach is NOT pseudo-ML with weights; it's BRR-based inference. Our implementation should match this: unweighted model fitting, weighted DISTRIB, BRR for SEs. Document clearly that `weight_var` in `mixtran()` affects DISTRIB computations, not the model likelihood.

3. **Balanced Repeated Replication weight generation.** SAS code for generating NHANES BRR weights uses specific Hadamard matrix constructions. The R `survey` package can generate BRR weights from NHANES strata/PSU variables. Provide a helper or document the workflow.

---

## 7. Validation Plan

### 7.1 Synthetic Data Recovery Tests (Unit Tests)

Already prototyped and passing. Formalize as `testthat` tests:

| Test | Pass Criterion | Prototype Result |
|------|---------------|-----------------|
| Amount model: σ²_b recovery | Within 15% of true | 5.3% error |
| Amount model: σ²_w recovery | Within 15% of true | 1.0% error |
| Amount model: β recovery | All coefficients within 2 SE | Pass |
| Amount model: mean usual intake | Within 10% of theoretical | 1.4% error |
| Amount model: distribution narrower than observed | SD_usual < SD_observed | Pass (ratio 0.79) |
| Uncorrelated: σ²_v1 recovery | Within 25% | 5.6% error |
| Uncorrelated: σ²_v2 recovery | Within 25% | 4.4% error |
| Uncorrelated: σ²_e recovery | Within 25% | 5.9% error |
| Box-Cox: λ recovery | Within 0.15 of true | 0.06 error |
| Box-Cox: log convention | λ < 0.15 → use log | Pass |
| Correlated: ρ recovery | Within 0.3 of true | TBD |

**Generate tests across a range of:**
- Sample sizes (n = 200, 500, 1000, 3000)
- Variance ratios (σ²_w/σ²_b = 0.5, 1.0, 2.0, 5.0, 10.0)
- Consumption prevalence for two-part models (30%, 50%, 70%, 90%)
- Correlation values for correlated model (-0.5, 0, 0.3, 0.6, 0.9)

### 7.2 Cross-Validation Against SAS Output

The gold standard. Requires running the SAS macros and our R code on identical NHANES subsets.

**Protocol:**
1. Download NHANES 2017-2020 pre-pandemic dietary data
2. Prepare identical analytical datasets in both SAS and R
3. Run NCI SAS macros (MIXTRAN + DISTRIB) for:
   - Sodium (amount-only, ubiquitous)
   - Added sugars (amount-only, ubiquitous)
   - Whole grains (two-part, episodic)
   - Whole fruit (two-part, episodic)
4. Run `nciusual::mixtran()` + `nciusual::distrib()` on same data
5. Compare:
   - λ estimates (should match within 0.02)
   - Variance components (should match within 10%)
   - Fixed effect coefficients (should match within 0.5 SE)
   - Percentiles of usual intake distribution (should match within 5%)
   - Cutpoint proportions (should match within 2 percentage points)

**NCI provides example datasets and output** at https://epi.grants.cancer.gov/diet/usualintakes/macros_single.html — use these as primary validation targets before moving to custom NHANES analyses.

### 7.3 Published Results Replication

Replicate the baseline intake estimates from Wang et al. (2023) Table 2 using NHANES 2013-2018 data as a secondary validation of the full pipeline (data prep + NCI method + survey weighting).

---

## 8. Package Dependencies

### Required (Imports)

| Package | Purpose | Why Not Optional |
|---------|---------|-----------------|
| `nlme` | Amount model and two-part amount component (lme) | Core model fitting engine |
| `MASS` | Two-part probability model (glmmPQL) | Core model fitting engine |
| `stats` | Optimization, distributions, model formulas | Base R, always available |

### Suggested (Suggests)

| Package | Purpose | When Needed |
|---------|---------|-------------|
| `lme4` | Alternative GLMM engine (glmer) for two-part probability | When `engine = "lme4"` |
| `TMB` | Template Model Builder for correlated model | When `engine = "tmb"` |
| `Rcpp` | C++ acceleration for GH quadrature | When `engine = "ghq"` |
| `survey` | BRR weight generation and complex survey support | For NHANES analyses |
| `nhanesA` | NHANES data download | For NHANES helper functions |
| `haven` | SAS/SPSS data import | For reading NHANES XPT files |
| `future.apply` | Parallel BRR replication | For parallel BRR |
| `testthat` | Testing framework | Development only |

### Dependency Rationale

Keep `Imports` minimal (`nlme`, `MASS`, `stats`) so the package installs easily on any R system without compiled code. All compiled-code dependencies (`TMB`, `Rcpp`, `lme4`) are optional enhancments via `Suggests`. The base package should be fully functional with just `Imports`.

---

## 9. File Structure

```
nciusual/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── README.md
├── R/
│   ├── boxcox.R              # Box-Cox transform, inverse, lambda search
│   ├── mixtran.R             # Core MIXTRAN: amount, uncorr, corr models
│   ├── distrib.R             # DISTRIB: Monte Carlo usual intake distributions
│   ├── indivint.R            # INDIVINT: individual predicted usual intake
│   ├── brr.R                 # BRR standard error computation
│   ├── nhanes_helpers.R      # NHANES data download and prep utilities
│   ├── quadrature.R          # Gauss-Hermite quadrature nodes/weights
│   ├── weighted_stats.R      # Weighted quantiles, means, proportions
│   └── print_methods.R       # print/summary/plot methods for fit objects
├── src/                      # (Phase 2+) C++ code for Rcpp/TMB acceleration
│   └── ghq_loglik.cpp
├── man/                      # roxygen2-generated documentation
├── tests/
│   ├── testthat/
│   │   ├── test-boxcox.R
│   │   ├── test-amount-model.R
│   │   ├── test-twopart-uncorr.R
│   │   ├── test-twopart-corr.R
│   │   ├── test-distrib.R
│   │   ├── test-brr.R
│   │   └── test-nhanes-helpers.R
│   └── testthat.R
├── vignettes/
│   ├── introduction.Rmd      # Package overview and quick start
│   ├── nhanes-example.Rmd    # Full NHANES analysis walkthrough
│   └── sas-comparison.Rmd    # Side-by-side comparison with SAS output
└── data-raw/
    └── nci_example_data.R    # Script to prepare bundled example dataset
```

---

## 10. Development Roadmap

### v0.1.0 — Core Functionality ✓ Released

- [x] Box-Cox transformation utilities (validated)
- [x] Amount-only model via `nlme::lme` (validated)
- [x] Two-part uncorrelated model via `glmmPQL` + `lme` (validated)
- [x] Two-part correlated model via profile likelihood over ρ
- [x] DISTRIB Monte Carlo simulation (validated)
- [x] Weighted percentile/mean computation with proper survey weight handling
- [x] BRR standard error wrapper (parallel-capable)
- [x] `print`, `summary`, `plot` methods for all return objects
- [x] `distrib_to_df()` convenience function for tidy output
- [x] `testthat` test suite covering all synthetic data recovery tests
- [x] NHANES example vignette
- [x] CRAN-compatible package structure (roxygen2, NAMESPACE, etc.)

### v0.2.0 — Enhanced Accuracy ✓ Released

- [x] `lme4::glmer` as alternative engine for probability model (`prob_engine = "glmer"`)
- [x] INDIVINT individual predictions (`indivint()` — BLUP-based, amount and two-part)
- [x] Starting values from prior model fit (`start` parameter to `mixtran()`)
- [x] NHANES multi-cycle weight helper (`combine_nhanes_cycles()`)
- [x] Vectorized DISTRIB simulation (matrix ops replace per-subject R loops; large speedup)
- [x] Pure-R GH quadrature engine for correlated model (`corr_engine = "ghq"`,
      `ghq_n_nodes` controls accuracy vs. speed; optimises σ²_v1, σ²_v2, σ²_e, ρ jointly)
- [x] GH quadrature utilities (`gh_nodes_weights`, `gh_nodes_normal`, `gh_nodes_bivariate`,
      `gh_nodes_adaptive`) in `R/quadrature.R`
- [x] Rcpp C++ stub for GHQ inner loop (`src/ghq_loglik.cpp`) — activatable via adding
      `Rcpp` to `LinkingTo`/`Imports` for a ~50-100× speedup on the GHQ optimisation
- [ ] Cross-validation against SAS macro output on NCI example datasets

### v0.3.0 — Full SAS Parity ✓ Released

- [x] Bivariate usual intake and ratio estimation — `distrib_bivariate()` with joint
      Monte Carlo simulation respecting cross-component random-effect correlation;
      supports arbitrary derived quantities via `fn_ratio` (e.g., sodium density,
      percent energy from fat)
- [x] Adaptive GH quadrature (`gh_nodes_adaptive()`) — mode-finding + Hessian scaling
      for subject-specific node placement, matching SAS PROC NLMIXED behaviour
- [ ] TMB engine for correlated model (`engine = "tmb"`) — requires C++ compilation
- [ ] Full replication of NCI example programs 1–4 (vignettes)

---

## 11. Prototype Code

A working prototype covering the amount-only model, uncorrelated two-part model, and DISTRIB was built and validated. The prototype files are included in this repository as reference for the implementation. All validation tests passed:

- `R/boxcox.R` — Box-Cox utilities (production-ready)
- `R/mixtran.R` — MIXTRAN with all three model types (amount and uncorr validated; corr needs performance work)
- `R/distrib.R` — DISTRIB Monte Carlo simulation (production-ready)
- `tests/test_validate.R` — Synthetic data validation script

These files can be used as the starting point. The correlated model implementation in the prototype (`fit_twopart_corr`) is statistically correct but too slow for production use due to unvectorized R loops in the GH quadrature. Replace with the profile likelihood approach for v0.1.

---

## 12. Key References

1. Tooze JA, et al. (2006). A new statistical method for estimating the usual intake of episodically consumed foods with application to their distribution. *JADA* 106(10):1575-87. — **Original NCI method paper.**

2. Kipnis V, et al. (2009). Modeling data with excess zeros and measurement error: application to evaluating relationships between episodically consumed foods and health outcomes. *Biometrics* 65(4):1003-10. — **INDIVINT method.**

3. Tooze JA, et al. (2010). A mixed-effects model approach for estimating the distribution of usual intake of nutrients: the NCI method. *Stat Med* 29(27):2857-68. — **Full statistical development.**

4. Herrick KA, et al. (2018). Estimating usual dietary intake from NHANES data using the NCI method. *Vital Health Stat* Series 2, No. 178. — **Definitive NHANES application guide. Covers BRR weight construction, all macro parameters, worked examples.**

5. Luo H, et al. (2021). Introduction to the SIMPLE Macro, a tool to increase the accessibility of 24-hour dietary recall analysis and modeling. *J Nutr* 151(5):1329-40. — **SIMPLE macro wrapping MIXTRAN + DISTRIB + BRR.**

6. Luo H, et al. (2022). Advanced dietary analysis and modeling: a deep dive into the NCI method. *Curr Dev Nutr* 6(10):nzac135. — **Most detailed published guide to advanced NCI macro usage. Essential reading for implementation.**

7. Arnold CD, et al. (2019). A new statistical method for estimating usual intakes of nearly-daily consumed foods and nutrients through use of only one 24-hour dietary recall. *J Nutr* 149(10):1667-73. — **NCI 1-day method (TRAN1 macro). Could be added as an extension.**

8. NCI SAS macro documentation: https://epi.grants.cancer.gov/diet/usualintakes/macros_single.html — **User's Guide (PDF), example programs, example datasets.**

9. Wang L, et al. (2023). Evaluation of health and economic effects of United States school meal standards. *AJCN* 118:605-13. — **The CRA analysis this package supports. Uses NCI method for baseline intake estimation.**
