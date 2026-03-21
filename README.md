# nciusual

An R package implementing the National Cancer Institute (NCI) method for estimating usual dietary intake distributions from 24-hour dietary recall data. R port of the NCI MIXTRAN v2.21, DISTRIB v2.2, and INDIVINT SAS macros.

## Overview

A single 24-hour dietary recall captures only one day of intake and conflates two sources of variation:

- **Between-person variation** — real differences in habitual (usual) diets
- **Within-person variation** — day-to-day fluctuation within an individual

The NCI method uses repeated recalls to estimate and remove within-person variation, recovering the narrower distribution of *usual* (long-run average) intake. This package implements that methodology entirely in R, making it accessible without a SAS license.

## Installation

```r
remotes::install_github("comcgovern/NCI-intakes-inR")
```

Core dependencies (`nlme`, `MASS`, `stats`, `survey`) are installed automatically. Optional enhancements: `lme4` (glmer engine), `future.apply` (parallel BRR).

## Quick Start

```r
library(nciusual)

# Step 1 — Fit MIXTRAN model
fit <- mixtran(
  data        = nhanes_data,
  intake_var  = "sodium_mg",
  subject_var = "SEQN",
  repeat_var  = "day_num",
  weekend_var = "weekend",
  weight_var  = "WTDR2D",
  model_type  = "amount"    # or "uncorr" / "corr" for episodic foods
)

# Step 2 — Estimate usual intake distribution (DISTRIB)
result <- distrib(fit, n_sims = 100, seed = 42)
print(result)

# Step 3 — BRR standard errors
brr_wts <- generate_nhanes_brr_weights(nhanes_data, weight_var = "WTDR2D",
                                        subject_var = "SEQN")
se_result <- brr_usual_intake(nhanes_data, brr_wts, subject_var = "SEQN",
                               mixtran_args = list(intake_var = "sodium_mg", ...))
```

## Core Functions

### Model Fitting — `mixtran()`

```r
mixtran(
  data, intake_var, subject_var, repeat_var,
  model_type  = "amount",       # "amount" | "uncorr" | "corr"
  covariates  = NULL,           # demographic covariates (retained in usual intake)
  weekend_var = NULL,           # weekend indicator (removed from usual intake)
  weight_var  = NULL,           # survey weights (used in DISTRIB, not likelihood)
  lambda      = NULL,           # Box-Cox lambda; NULL = auto-estimated
  prob_engine = "glmmPQL",      # "glmmPQL" | "glmer" (two-part models)
  corr_engine = "profile_rho",  # "profile_rho" | "ghq" (corr model only)
  ghq_n_nodes = 5L,             # GH nodes per dimension (3/5/7/9) for GHQ engine
  start       = NULL,           # prior mixtran_fit for warm starting
  verbose     = TRUE
)
```

**Model types:**

| `model_type` | Use for | Notes |
|---|---|---|
| `"amount"` | Nutrients consumed every day (sodium, energy, protein) | `nlme::lme` |
| `"uncorr"` | Episodic foods, independent random effects | `glmmPQL`/`glmer` + `lme` |
| `"corr"` | Episodic foods, correlated random effects (full NCI parity) | Profile-ρ or GHQ |

**Correlated model engines:**

| `corr_engine` | Method | Speed | Accuracy |
|---|---|---|---|
| `"profile_rho"` | Profile likelihood over ρ grid (default) | Fast | Approximate |
| `"ghq"` | Gauss-Hermite quadrature, joint optimisation of σ²_v1/v2/e and ρ | Moderate | More accurate |

### Distribution Estimation — `distrib()`

```r
result <- distrib(
  mixtran_obj,
  subgroup_var = NULL,           # optional stratification variable
  cutpoints    = NULL,           # intake thresholds for % below/above
  percentiles  = c(5,10,25,50,75,90,95),
  n_sims       = 100,
  seed         = 12345
)

distrib_to_df(result)            # tidy data frame of results
```

### Individual Predictions — `indivint()`

BLUP-based individual usual intake predictions for regression calibration:

```r
preds <- indivint(fit, zero_seq = TRUE, zero_weekend = TRUE)
# Returns: subject, predicted_usual, predicted_usual_orig
#          (+ prob_usual, amt_usual_orig for two-part models)
```

### Bivariate Usual Intake — `distrib_bivariate()`

Joint estimation for two dietary components with cross-component correlation, enabling ratio and derived-quantity distributions:

```r
biv <- distrib_bivariate(
  data        = nhanes_data,
  intake_var1 = "sodium_mg",
  intake_var2 = "energy_kcal",
  subject_var = "SEQN",
  repeat_var  = "day_num",
  fn_ratio    = function(u1, u2) (u1 / u2) * 1000,  # sodium density (mg/1000 kcal)
  ratio_label = "sodium_density",
  n_sims      = 100
)
print(biv)
bivariate_to_df(biv)
```

### BRR Standard Errors — `brr_usual_intake()`

Balanced Repeated Replication with Fay's modification, following NHANES convention:

```r
# Generate BRR weights from NHANES design variables
brr_wts <- generate_nhanes_brr_weights(
  data        = nhanes_data,
  strata_var  = "SDMVSTRA",
  psu_var     = "SDMVPSU",
  weight_var  = "WTDR2D",
  subject_var = "SEQN"
)

# Run BRR pipeline (parallel-capable)
se_result <- brr_usual_intake(
  data            = nhanes_data,
  replicate_weights = brr_wts,
  subject_var     = "SEQN",
  mixtran_args    = list(intake_var = "sodium_mg", ...),
  parallel        = TRUE
)
brr_to_df(se_result)   # estimate, SE, 95% CI for every statistic
```

### Multi-Cycle NHANES — `combine_nhanes_cycles()`

```r
pooled <- combine_nhanes_cycles(
  cycle_list  = list(nhanes_1718, nhanes_1920, nhanes_2122),
  weight_var  = "WTDR2D",
  subject_var = "SEQN"
)
```

### Box-Cox Utilities

```r
boxcox_transform(y, lambda)          # T(y; λ)
boxcox_inverse(z, lambda)            # T⁻¹(z; λ)
find_optimal_lambda(y, weights, ...)  # profile likelihood grid search
```

### Gauss-Hermite Quadrature — `R/quadrature.R`

```r
gh_nodes_weights(n)                          # physicist convention (Σw = √π)
gh_nodes_normal(n)                           # probabilist (Σw = 1)
gh_nodes_bivariate(n, rho, sigma_v1, sigma_v2)   # tensor-product BVN nodes
gh_nodes_adaptive(n, eta_i, mu_i, ...)       # adaptive nodes (mode + Hessian scale)
```

## Statistical Background

The NCI method is a two-step procedure:

1. **MIXTRAN** — Fits a Box-Cox transformed mixed model, decomposing observed intake variation into between-person (σ²_b) and within-person (σ²_w) components.

2. **DISTRIB** — Monte Carlo simulation: for each person, draws `n_sims` values of the between-person random effect, back-transforms to the original scale, and summarises the resulting usual intake distribution. Within-person error is intentionally omitted — this is what narrows the distribution.

For episodically consumed foods, a two-part model separates the probability of consumption from the amount consumed on consumption days. The correlation ρ between their subject-level random effects is estimated via profile likelihood (fast) or Gauss-Hermite quadrature (more accurate).

For bivariate analyses, the cross-component correlation between random effects is estimated from empirical Bayes predictions and propagated through the joint Monte Carlo simulation.

See [`PACKAGE_SPEC.md`](PACKAGE_SPEC.md) for the full statistical specification, parity matrix with SAS macros, and implementation notes.

## Repository Structure

```
nciusual/
├── DESCRIPTION
├── NAMESPACE
├── PACKAGE_SPEC.md          # Full statistical and API specification
├── README.md
├── R/
│   ├── boxcox.R             # Box-Cox transform, inverse, lambda search
│   ├── mixtran.R            # MIXTRAN: amount, uncorr, corr, GHQ engine
│   ├── distrib.R            # DISTRIB: Monte Carlo usual intake distributions
│   ├── indivint.R           # INDIVINT: individual BLUP predictions
│   ├── brr.R                # BRR: standard errors + NHANES weight helpers
│   ├── bivariate.R          # Bivariate usual intake and ratio estimation
│   ├── quadrature.R         # GH quadrature utilities (nodes/weights/adaptive)
│   └── print_methods.R      # print/summary/plot methods
└── tests/testthat/
    ├── test-amount-model.R
    ├── test-boxcox.R
    ├── test-brr.R
    ├── test-distrib.R
    ├── test-ghq-engine.R
    ├── test-glmer-engine.R
    ├── test-indivint.R
    ├── test-bivariate.R
    ├── test-nhanes-helpers.R
    ├── test-twopart-corr.R
    └── test-twopart-uncorr.R
```

## Validation

Key results from synthetic-data recovery tests:

| Check | Target | Result |
|---|---|---|
| Between-person variance (σ²_b) recovery | < 15% error | 5.3% |
| Within-person variance (σ²_w) recovery | < 15% error | 1.0% |
| Mean usual intake recovery | < 10% error | 1.4% |
| Usual intake SD narrower than observed | SD ratio < 1 | 0.79 |
| Uncorrelated model: σ²_v1/v2/e recovery | < 25% error | < 6% |

Run the full test suite:

```r
testthat::test_dir("tests/testthat")
```

## SAS Macro Parity

| SAS Macro Feature | R Implementation | Status |
|---|---|---|
| MIXTRAN amount-only model | `nlme::lme` | Full |
| MIXTRAN uncorrelated two-part | `glmmPQL`/`glmer` + `lme` | Approximate (PQL) / Full (glmer) |
| MIXTRAN correlated two-part | Profile-ρ or GHQ | Approximate / Improved |
| MIXTRAN Box-Cox lambda search | Profile likelihood grid | Full |
| DISTRIB Monte Carlo simulation | Vectorized R | Full |
| DISTRIB weighted percentiles/means | Custom weighted quantile | Full |
| INDIVINT BLUP predictions | Empirical Bayes from lme/glmer | Full |
| BRR standard errors | Parallel replicate pipeline | Full |
| NHANES BRR weight generation | `survey` package wrapper | Full |
| Multi-cycle weight pooling | `combine_nhanes_cycles()` | Full |
| Bivariate usual intake | `distrib_bivariate()` | Implemented |

## Changelog

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
