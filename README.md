# nciusual

An R package implementing the National Cancer Institute (NCI) method for estimating usual dietary intake distributions from 24-hour dietary recall data. It is a port of the NCI **MIXTRAN, DISTRIB and INDIVINT SAS macros, version 2.1**, and is validated against NCI's published SAS output for the macro example programs (see [Validation](#validation)).

## Overview

A single 24-hour dietary recall captures only one day of intake and conflates two sources of variation:

- **Between-person variation**: real differences in habitual (usual) diets
- **Within-person variation**: day-to-day fluctuation within an individual

The NCI method uses repeated recalls to separate the two and estimate the distribution of *usual* (long-run average) intake. This package does that entirely in R, without a SAS license.

## Installation

```r
remotes::install_github("comcgovern/NCI-intakes-inR")
```

The package has no compiled code and depends only on base R packages. Optional: `survey` (BRR weight generation), `future.apply` (parallel BRR), `haven` (reading NCI's example data).

## Quick Start

```r
library(nciusual)

# Step 1 - MIXTRAN: fit the model (survey weights enter the likelihood,
# as with MIXTRAN's replicate_var; lambda is estimated)
fit <- mixtran(
  data        = nhanes_data,
  intake_var  = "sodium_mg",
  subject_var = "SEQN",
  repeat_var  = "day_num",
  weekend_var = "weekend",
  weight_var  = "WTDR2D",
  covariates  = c("age_19_30", "age_31_50", "female"),
  model_type  = "amount"    # or "uncorr" / "corr" for episodic foods
)

# Step 2 - DISTRIB: usual intake distribution
result <- distrib(fit, n_sims = 100, seed = 42, cutpoints = 2300)
print(result)

# Step 3 - BRR standard errors (base weight in mixtran_args)
brr_wts <- generate_nhanes_brr_weights(nhanes_data, weight_var = "WTDR2D",
                                       subject_var = "SEQN")
se_result <- brr_usual_intake(
  nhanes_data, brr_wts, subject_var = "SEQN",
  mixtran_args = list(intake_var = "sodium_mg", subject_var = "SEQN",
                      repeat_var = "day_num", weekend_var = "weekend",
                      weight_var = "WTDR2D"))
```

## Core Functions

### Model Fitting: `mixtran()`

```r
mixtran(
  data, intake_var, subject_var, repeat_var,
  model_type      = "amount",    # "amount" | "uncorr" | "corr"  (SAS: amount | nocorr | corr)
  covariates      = NULL,        # amount-part covariates (SAS covars_amt)
  covariates_prob = covariates,  # probability-part covariates (SAS covars_prob)
  weekend_var     = NULL,        # 0/1 Fri-Sun indicator (SAS weekend=)
  weight_var      = NULL,        # survey weight (SAS replicate_var)
  include_seq     = TRUE,        # sequence indicator for recalls after the first (SAS seq=)
  lambda          = NULL,        # NULL = estimated by ML; a value fixes it
  subgroup_var    = NULL,        # carried to the predictions for distrib()
  start           = NULL,        # prior fit for starting values (SAS start_val1-3)
  verbose         = TRUE
)
```

The models are fitted the way MIXTRAN fits them with PROC NLMIXED:

- **Weights** multiply each subject's marginal log-likelihood (the NLMIXED `REPLICATE` statement). Multiplying all weights by a constant changes nothing.
- **Box-Cox λ** is a free parameter of the amount model, bounded below by 0.01 (`A_LAMBDA > 0.01`). Its starting value comes from a TRANSREG-style grid search with the amount covariates.
- **Uncorrelated model**: the probability part (logistic, random intercept) and the amount part (Box-Cox normal on consumption days) are fitted separately.
- **Correlated model**: the joint likelihood is maximised starting from the uncorrelated estimates. The random effects are integrated out: the amount random effect in closed form, the probability random effect by adaptive Gauss–Hermite quadrature (`n_quad`, default 21).
- Sequence and weekend effects enter both parts of the model.

### Distribution Estimation: `distrib()`

```r
result <- distrib(
  mixtran_obj,
  subgroup_var  = NULL,          # in fit$predicted, or supply subgroup_data (SAS add_da)
  cutpoints     = NULL,          # proportion below/above each cutpoint
  percentiles   = c(5, 10, 25, 50, 75, 90, 95),
  n_sims        = 100,           # SAS nsim_mc
  seed          = 12345,
  wkend_prop    = 3/7            # SAS wkend_prop (weekday weight 4/7)
)
distrib_to_df(result)
```

As in the DISTRIB macro, each subject's usual intake is computed from simulated random effects as:

- the consumption probability, times
- the consumption-day amount *integrated over within-person error*: DISTRIB's 9-point rule, or exact `exp(μ + σ²_e/2)` when λ = 0, bounded below by half the minimum amount.

The sequence effect is removed. With a weekend indicator, weekday and weekend usual intake are averaged 4/7 : 3/7. Percentiles and cutpoint proportions use the DISTRIB PC-macro interpolation.

### Individual Predictions: `indivint()`

```r
preds <- indivint(fit, lamt = NULL)   # lamt: optional Box-Cox scale for the prediction
# Returns: subject, indusint, predicted_usual_orig (+ prob_usual, amt_usual_orig)
```

As in the INDIVINT macro (Kipnis et al., 2009), this returns each subject's conditional expected usual intake given their recalls, `E[T | recalls]`, computed by adaptive quadrature. It is intended for regression calibration.

### Bivariate Usual Intake: `distrib_bivariate()`

```r
biv <- distrib_bivariate(
  data = nhanes_data, intake_var1 = "sodium_mg", intake_var2 = "energy_kcal",
  subject_var = "SEQN", repeat_var = "day_num",
  fn_ratio = function(u1, u2) (u1 / u2) * 1000, ratio_label = "sodium_density")
```

This is **not** a port of NCI's multivariate macros. It fits each component separately, approximates the cross-component correlation from the empirical Bayes estimates, and simulates both components jointly with DISTRIB's back-transformation.

### BRR Standard Errors: `brr_usual_intake()`

This uses Balanced Repeated Replication with Fay's modification (F = 0.3). The full-sample run uses `mixtran_args$weight_var` (the base weight). Each replicate refits MIXTRAN and DISTRIB with its replicate weight, starting from the full-sample estimates:

```r
brr_wts <- generate_nhanes_brr_weights(nhanes_data, strata_var = "SDMVSTRA",
                                       psu_var = "SDMVPSU", weight_var = "WTDR2D",
                                       subject_var = "SEQN")
se_result <- brr_usual_intake(nhanes_data, brr_wts, subject_var = "SEQN",
                              mixtran_args = list(..., weight_var = "WTDR2D"),
                              parallel = TRUE)
brr_to_df(se_result)
```

### Multi-Cycle NHANES: `combine_nhanes_cycles()`

```r
pooled <- combine_nhanes_cycles(list(nhanes_1516, nhanes_p1720),
                                weight_var = "WTDR2D", subject_var = "SEQN",
                                cycle_years = c(2, 3.2))   # unequal cycle lengths
```

Without `cycle_years`, weights are divided by the number of cycles, which is correct only for cycles of equal length. Check current CDC NHANES weighting guidance for the cycles you pool.

### Box-Cox and Quadrature Utilities

```r
boxcox_transform(y, lambda); boxcox_inverse(z, lambda)
find_optimal_lambda(y, weights, lambda_grid, covariates)   # grid search (starting value)
gh_nodes_weights(n); gh_nodes_normal(n)
gh_nodes_bivariate(n, rho, sigma_v1, sigma_v2); gh_nodes_adaptive(...)
```

## Statistical Background

1. **MIXTRAN** fits a Box-Cox transformed mixed model. Amount-only: `T(Y_ij) = X'β + u_i + ε_ij`. Two-part: `logit P(Y_ij > 0) = Z'α + u1_i` together with `T(Y_ij | Y_ij > 0) = X'β + u2_i + ε_ij`, where `(u1, u2)` is bivariate normal (correlated model) or independent (uncorrelated model).
2. **DISTRIB** draws the random effects by Monte Carlo. Usual intake is `P(u1) × E_ε[T⁻¹(X'β + u2 + ε)]`: the within-person error is *integrated over* in the back-transformation, not drawn. That is what removes day-to-day variation while keeping usual intake on the scale of expected daily intake.
3. **INDIVINT** computes `E[usual intake | person's recalls]`.

## Validation

### Parity with NCI's published SAS output

NCI publishes the v2.1 macros, the analytic dataset, and the SAS listing output for its example programs ([macros](https://epi.grants.cancer.gov/diet/usualintakes/macros_single.html), [data](https://epi.grants.cancer.gov/diet/usualintakes/dataset.html)). The package was run on that dataset with the same model specifications:

| | Result |
|---|---|
| **Example 1** (added sugars, males 9+, amount model) | λ 0.3939 vs SAS 0.3941; σ²_u 3.514 vs 3.518; σ²_e 4.072 vs 4.077; fixed effects within 0.02; −2 log L equal to SAS's printed value. DISTRIB mean, P1–P99 and all 10 cutpoint proportions, overall and by age group, within about 1–2% (simulation noise). |
| **Example 2** (total fruit, children 1–8, correlated model) | Uncorrelated stage matches SAS to 3–4 decimals (λ 0.4086 vs 0.4086). Correlated stage: λ 0.4108 vs 0.4109, ρ 0.544 vs 0.556, σ²_u1 3.12 vs 3.00. DISTRIB mean and P10–P99 within about 1%; P1/P5 2–3% lower, within the simulation error of SAS's single 100-draw run; cutpoint proportions within 0.005. |

Where parameters differ, nciusual's estimates have a *higher* likelihood than SAS's. SAS's NLMIXED stops at a relative-gradient criterion (its reported gradients are non-zero), and the correlated model's likelihood is very flat in the (σ²_u1, ρ) direction. Reproduce with:

```sh
Rscript data-raw/sas_parity_examples.R                       # side-by-side report
NCIUSUAL_SAS_PARITY=true Rscript -e 'testthat::test_file("tests/testthat/test-sas-examples.R")'
```

SAS's Example 4 (INDIVINT) dataset is distributed in SAS CPORT format, which R cannot read. INDIVINT is instead checked against brute-force numerical integration of the same conditional expectation; it agrees to 6 decimals.

### Other checks

- The marginal likelihoods agree with brute-force numerical integration (amount, probability and correlated models), and the probability part agrees with `lme4::glmer(nAGQ = 25)`.
- Analytic ground-truth fixtures (`data-raw/generate_reference.R`, `tests/testthat/fixtures/`): datasets simulated from known parameters, with the true usual-intake distribution under the NCI definition. All three model types recover it within 4–7% across P5–P95.

```r
testthat::test_dir("tests/testthat")
```

## SAS Macro Parity

| SAS macro feature | R implementation | Status |
|---|---|---|
| MIXTRAN amount-only model (NLMIXED, REPLICATE weights, λ estimated) | `mixtran(model_type = "amount")` | Matches Example 1 |
| MIXTRAN uncorrelated two-part (`nocorr`) | `mixtran(model_type = "uncorr")` | Matches Example 2 nocorr stage |
| MIXTRAN correlated two-part (`corr`) | `mixtran(model_type = "corr")` | Matches Example 2 (see Validation) |
| `covars_prob` / `covars_amt`, `seq`, `weekend`, `subgroup` | `covariates_prob`, `covariates`, `include_seq`, `weekend_var`, `subgroup_var` | Implemented |
| `vargroup` (separate residual variance by group) | not implemented | Gap |
| `start_val1-3` | `start` | Implemented |
| DISTRIB Monte Carlo, 9-point back-transform, min-amount floor | `distrib()` | Matches Examples 1–2 |
| DISTRIB `wkend_prop`, subgroups, cutpoints, percentiles | `distrib()` | Implemented |
| DISTRIB `recamt` (proportion meeting a per-person recommendation) | not implemented | Gap |
| INDIVINT (adaptive-quadrature conditional expectation, `lamt`) | `indivint()` | Verified numerically |
| BRR with Fay's modification | `brr_usual_intake()` | Implemented |
| NCI multivariate (bivariate/ratio) macros | `distrib_bivariate()` | Approximation, not a port |

## Changelog

### v0.5.0 (2026-10-14)

**SAS parity rewrite.** An audit against the NCI v2.1 macro source and NCI's published example output found material differences; see [`docs/AUDIT_2026-10.md`](docs/AUDIT_2026-10.md). Results from earlier versions are not comparable to SAS: DISTRIB means were about 10% low and lower-tail percentiles 25–40% low in Example 1.

- **MIXTRAN**: new weighted maximum-likelihood engine mirroring PROC NLMIXED. It applies subject-level `REPLICATE` weights, estimates the Box-Cox λ in the likelihood, fits the uncorrelated parts separately and the correlated model jointly, and supports separate probability/amount covariates (`covariates_prob`). This replaces `nlme::lme` with `varFixed` weights (which mis-scaled σ²_w and failed with NHANES-scale weights), PQL/glmer, profile-ρ and the non-adaptive GHQ engine. `prob_engine`, `corr_engine` and `ghq_n_nodes` are deprecated and ignored. `nlme` and `MASS` are no longer dependencies.
- **DISTRIB**: integrates the back-transformation over within-person error (9-point rule) and floors amounts at half the minimum amount. It averages weekday/weekend intake (`wkend_prop`, default 3/7), predicts the amount part from each never-consumer's own covariates, and uses the PC-macro percentile and cutpoint interpolation. It also handles factor subgroups, accepts `subgroup_data`, and no longer changes the caller's RNG state.
- **INDIVINT**: returns the conditional expectation `E[T | recalls]` by adaptive quadrature, including correlation, within-person error and the `lamt` option. The meaningless two-part `predicted_usual` column (a logit × Box-Cox product) is gone; `predicted_usual` is now an alias of `indusint`.
- **BRR**: the full-sample run uses the base weight (previously unweighted while replicates were weighted, inflating SEs). Replicates warm-start from the full-sample fit, and statistics are aligned by name.
- `combine_nhanes_cycles()` gains `cycle_years` for cycles of unequal length.
- `gh_nodes_adaptive()` weights now include the prior/standard-normal density ratio. The 7-point Gauss–Hermite node table is corrected.
- Data sufficiency follows MIXTRAN: stop when fewer than 2 subjects have repeat (positive) recalls, warn for 2–10. `skip_if_empty = FALSE` now errors instead of returning NULL.
- Tests: new SAS-example parity tests, brute-force likelihood/INDIVINT checks, and a `glmer` cross-check. Ground-truth fixtures are regenerated with the NCI definition of usual intake.

### v0.4.1 (2026-06-15)

**Correctness fixes (code audit)**

- **Random-effect extraction bug (major)** — `nlme::ranef(fit)[[1]]` was
  collapsing the per-subject random-effect vector to a single unnamed value,
  so every subject lookup missed. This silently broke three features:
  - the correlated two-part model (`corr_engine = "profile_rho"`) always fell
    back to `rho = 0`, regardless of the true correlation;
  - two-part `indivint()` returned the population-level prediction for *every*
    subject (the BLUP step was inert);
  - bivariate `cross_rho` was always estimated as 0.

  All three now recover the intended quantities. Added regression tests that
  assert non-degenerate `rho`, per-subject variation in `indivint()`
  predictions, and positive `cross_rho` for correlated components.
- **Weighted quantile convention** — `weighted_quantiles()` switched from a
  right-edge cumulative-weight convention (which biased DISTRIB percentiles
  downward) to the midpoint/Hazen convention, so equal-weight percentiles now
  match base R `quantile(..., type = 5)`.
- **`plot.distrib_result()`** — replaced the invalid colour name `"teal"`
  (which errored for exactly 7 subgroups) with `"darkcyan"`.
- **Gauss-Hermite tables** — added the previously missing precomputed node/
  weight tables for `n = 5` (the default) and `n = 9`, avoiding a repeated
  eigendecomposition inside the GHQ likelihood loop; removed dead code in
  `gh_nodes_adaptive()`.

**Accuracy improvement**

- **GHQ correlated engine now bias-corrects the intercepts** — the amount
  sub-model is fit on consumed days only, so under positive ρ those days
  over-represent high-amount subjects and the consumers-only intercept is
  biased upward; the previous engines left it uncorrected, inflating the
  lower-tail percentiles of the correlated two-part usual-intake distribution
  by 20–40% against analytic ground truth. The `corr_engine = "ghq"` path now
  estimates the probability/amount intercept shifts jointly with the variance
  components and ρ (NLMIXED-style), recovering the *full* distribution —
  including the lower tail — to within ~3%. The faster `profile_rho` default
  remains approximate in the tail; use `"ghq"` when tail accuracy matters.

### v0.4.0 (2026-03-21)

**Robustness and compatibility improvements**

- **Pure R implementation** — removed Rcpp/C++ dependency; the GHQ inner loop is now implemented entirely in R, eliminating the need for Rtools/a C compiler on any platform.
- **`skip_if_empty` parameter** — `mixtran()` gains `skip_if_empty` (default `TRUE`) for graceful handling of dietary components with all-zero or missing intake values, returning `NULL` instead of stopping.
- **nlme compatibility fixes** — resolved a series of `base_lme`/namespace resolution errors when calling `nlme::lme()` across different nlme versions; the package now uses a wrapper approach that is robust to `match.call()` internals.
- **Degenerate model handling** — `mixtran()` now guards against degenerate two-part (corr/uncorr) fits (singular variance, constant covariates, MEEM errors) with informative warnings rather than hard stops.
- **`subgroup_var` carry-through** — fixed a bug in `distrib()` where the subgroup variable was not correctly merged into predicted usual intakes for subgroup-stratified analyses.
- **Empty random-effect vectors** — fixed a crash when `prob_re` vectors are empty or subject names are misaligned in two-part model predictions.
- **Profile-ρ boundary fix** — corrected profile likelihood degeneracy that produced ρ estimates at the grid boundary (±1) for near-degenerate correlated models.
- **New test coverage** — added `tests/testthat/test-nhanes-helpers.R` covering NHANES BRR weight generation and multi-cycle pooling utilities.

### v0.3.0

- GHQ engine for correlated two-part model (`corr_engine = "ghq"`)
- Bivariate usual intake and ratio estimation (`distrib_bivariate()`)
- Adaptive Gauss-Hermite quadrature utilities (`gh_nodes_adaptive()`)
- Vectorized DISTRIB Monte Carlo simulation

### v0.2.0

- `glmer` probability sub-model engine (`prob_engine = "glmer"`)
- Individual BLUP predictions (`indivint()`)
- Warm-starting via `start` parameter in `mixtran()`
- Multi-cycle NHANES weight pooling (`combine_nhanes_cycles()`)

### v0.1.0

- Initial release: `mixtran()`, `distrib()`, `brr_usual_intake()` for amount-only and uncorrelated two-part models

## References

1. Tooze JA et al. (2006). A new statistical method for estimating the usual intake of episodically consumed foods. *JADA* 106(10):1575–87.
2. Tooze JA et al. (2010). A mixed-effects model approach for estimating the distribution of usual intake of nutrients. *Stat Med* 29(27):2857–68.
3. Kipnis V et al. (2009). Modeling data with excess zeros and measurement error. *Biometrics* 65(4):1003–10.
4. Herrick KA et al. (2018). Estimating usual dietary intake from NHANES data using the NCI method. *Vital Health Stat* Series 2, No. 178.
5. NCI SAS macro documentation: <https://epi.grants.cancer.gov/diet/usualintakes/macros.html>

## License

MIT © 2026 Conor McGovern
