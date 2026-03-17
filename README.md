# nciusual

An R package implementing the National Cancer Institute (NCI) method for estimating usual dietary intake distributions from 24-hour dietary recall data. This is an R port of the NCI's MIXTRAN v2.21 and DISTRIB v2.2 SAS macros.

## Overview

A single 24-hour dietary recall captures only one day of intake and conflates two sources of variation:

- **Between-person variation** — real differences in habitual (usual) diets
- **Within-person variation** — day-to-day fluctuation within an individual

The NCI method uses repeated recalls to estimate and remove within-person variation, recovering the narrower distribution of *usual* (long-run average) intake. This package implements that methodology entirely in R, making it accessible without a SAS license.

## Installation

```r
# Install from source
install.packages("remotes")
remotes::install_github("comcgovern/NCI-intakes-inR")
```

Dependencies (`nlme`, `MASS`, `stats`, `survey`) are installed automatically.

## Usage

### 1. Fit a MIXTRAN model

`mixtran()` fits a mixed-effects model to long-format dietary recall data and separates between- and within-person variance components.

```r
library(nciusual)

# Amount-only model (ubiquitously consumed nutrients: sodium, energy, protein)
fit <- mixtran(
  data       = nhanes_data,
  intake_var = "sodium_mg",
  id_var     = "SEQN",
  recall_var = "day_num",
  covariates = c("age", "sex", "weekend"),
  model_type = "amount"
)

# Two-part model (episodically consumed foods: whole grains, alcohol)
fit <- mixtran(
  data       = nhanes_data,
  intake_var = "wholegrains_g",
  id_var     = "SEQN",
  recall_var = "day_num",
  covariates = c("age", "sex", "weekend"),
  model_type = "uncorr"   # or "corr" for correlated random effects
)
```

**Model types:**

| `model_type` | Use for | Random effects |
|---|---|---|
| `"amount"` | Nutrients consumed every day (sodium, energy) | Correlated intercepts |
| `"uncorr"` | Episodically consumed foods, simpler structure | Uncorrelated |
| `"corr"` | Episodically consumed foods, full NCI parity | Correlated |

### 2. Estimate the usual intake distribution

`distrib()` takes a fitted `mixtran_fit` object and performs Monte Carlo simulation to estimate the population-level usual intake distribution.

```r
result <- distrib(
  fit       = fit,
  n_sims    = 100,
  seed      = 42,
  subgroup  = "sex"   # optional stratification
)

print(result)
#> Usual Intake Distribution
#> Model type: amount
#> N subjects: 2847
#>
#> Percentiles of usual intake:
#>    P5   P10   P25   P50   P75   P90   P95
#>  1203  1487  1894  2541  3312  4089  4621
#>
#> Mean usual intake: 2612 mg/day

# Convert to tidy data frame
df <- distrib_to_df(result)
```

### 3. Box-Cox transformation utilities

```r
# Apply Box-Cox transformation
z <- boxcox_transform(y = intake_values, lambda = 0.30)

# Back-transform
y_hat <- boxcox_inverse(z = z, lambda = 0.30)

# Find optimal lambda (NCI convention: use log if lambda < 0.15)
lambda_opt <- find_optimal_lambda(
  y          = positive_intakes,
  weights    = survey_weights,
  covariates = covariate_matrix
)
```

## Statistical Background

The NCI method is a two-step procedure:

1. **MIXTRAN** — Fits a Box-Cox transformed mixed model to observed recall data, estimating mean intake on a transformed scale along with between-person (σ²_b) and within-person (σ²_w) variance components.

2. **DISTRIB** — Draws Monte Carlo samples from the fitted model, back-transforms to the original scale, and summarises the resulting usual intake distribution (percentiles, means, prevalence above/below cutpoints).

For episodically consumed foods, a two-part model is used:

- **Part 1 (probability):** Logistic GLMM estimating P(consumption on a given day)
- **Part 2 (amount | consumed):** Linear mixed model on Box-Cox transformed positive intakes

The correlation between Part 1 and Part 2 random effects is optionally modelled (`model_type = "corr"`).

For more detail, see [`PACKAGE_SPEC.md`](PACKAGE_SPEC.md) and the references below.

## Validation

`tests/test_validate.R` contains synthetic-data recovery tests. Key results from prototype validation:

| Check | Target | Result |
|---|---|---|
| Between-person variance (σ²_b) recovery | < 15% error | 5.3% |
| Within-person variance (σ²_w) recovery | < 15% error | 1.0% |
| Mean usual intake recovery | < 10% error | 1.4% |
| Usual intake SD narrower than observed | ratio ≈ 0.79 | confirmed |

Run the tests with:

```r
source("tests/test_validate.R")
```

## Repository Structure

```
NCI-intakes-inR/
├── DESCRIPTION          # R package metadata
├── NAMESPACE            # Exported functions
├── LICENSE              # MIT license
├── README.md            # This file
├── PACKAGE_SPEC.md      # Full statistical and API specification
├── R/
│   ├── boxcox.R         # Box-Cox transformation utilities
│   ├── mixtran.R        # MIXTRAN mixed-model fitting
│   └── distrib.R        # DISTRIB Monte Carlo simulation
└── tests/
    └── test_validate.R  # Synthetic-data validation tests
```

## References

- Tooze JA et al. (2006). A new statistical method for estimating the usual intake of episodically consumed foods with application to their distribution. *J Am Diet Assoc*, 106(10):1575–87.
- Tooze JA et al. (2010). Advanced methodological issues in the dietary assessment of episodically consumed foods. *Stat Methods Med Res*, 19(1):42–59.
- Kipnis V et al. (2009). Modeling data with excess zeros and measurement error. *Biometrics*, 65(4):1003–10.
- NCI MIXTRAN/DISTRIB macro documentation: <https://prevention.cancer.gov/research-groups/biometry/measurement-error-impact/software-measurement-error/usual-dietary-intake>

## License

MIT © 2026 Conor McGovern
