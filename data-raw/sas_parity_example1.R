# =============================================================================
# SAS parity check against NCI's published Example 1 (MIXTRAN/DISTRIB v2.1)
# =============================================================================
#
# Reproduces the comparison reported in docs/AUDIT_2026-10.md.
#
# NCI publishes the v2.1 macros, the analytic dataset, and the SAS listing
# (.lst) output for "Example 1" (added sugars, males 9+, amount-only model,
# covariates age4-age8 race1 race3, seq = day2, weekend, replicate_var = rndw1,
# DISTRIB nsim_mc = 100, wkend_prop = 3/7, subgroup = agegrp):
#   https://epi.grants.cancer.gov/diet/usualintakes/macros_single.html
#   https://epi.grants.cancer.gov/diet/usualintakes/dataset.html
#
# This script:
#   1. downloads the NCI analytic dataset (indiv_ims.sas7bdat);
#   2. fits the amount model the way MIXTRAN does (subject-level weighted ML
#      via the NLMIXED REPLICATE statement, lambda estimated jointly) and shows
#      it reproduces the SAS parameter estimates;
#   3. takes SAS's OWN published parameter estimates and runs
#        (a) a DISTRIB written from the v2.1 macro source, and
#        (b) the simulation logic used by nciusual::distrib(),
#      comparing both to the SAS .lst percentiles. This isolates the
#      DISTRIB step from any model-fitting differences.
#
# Requires: haven, nciusual. Run from the package root:
#   Rscript data-raw/sas_parity_example1.R
# =============================================================================

suppressMessages(library(nciusual))

# --- 1. Data -----------------------------------------------------------------
dl_dir <- file.path(tempdir(), "nci_example")
dir.create(dl_dir, showWarnings = FALSE)
zip_path <- file.path(dl_dir, "indiv_ims.zip")
if (!file.exists(zip_path)) {
  utils::download.file(
    "https://epi.grants.cancer.gov/diet/usualintakes/indiv_ims.zip",
    zip_path, mode = "wb")
}
utils::unzip(zip_path, exdir = dl_dir)
d <- haven::read_sas(file.path(dl_dir, "indiv_ims.sas7bdat"))
names(d) <- tolower(names(d))
m <- as.data.frame(d[d$stratum == 2,
  c("seqn", "day", "day2", "add_sug", "age4", "age5", "age6", "age7", "age8",
    "weekend", "rndw1", "agegrp", "race")])
m$race1 <- as.numeric(m$race == 1)
m$race3 <- as.numeric(m$race == 3)
m <- as.data.frame(lapply(m, as.numeric))
covs <- c("age4", "age5", "age6", "age7", "age8", "race1", "race3")

# --- Helpers written from the v2.1 macro source ------------------------------
bc <- function(y, l) if (l == 0) log(y) else (y^l - 1) / l

# MIXTRAN amount model: T(y) = Xb + u + e with the REPLICATE weight applied to
# each subject's log-likelihood, lambda a free parameter (bounds A_LAMBDA>0.01),
# Jacobian (lambda-1)*log(y) included, zeros -> half the minimum positive amount.
fit_amount_sas <- function(d, resp, subj, X_names, wt, lambda_start = 0.4) {
  d <- d[order(d[[subj]], d$day), ]
  minpos <- min(d[[resp]][d[[resp]] > 0])
  d[[resp]][d[[resp]] == 0] <- 0.5 * minpos
  y <- d[[resp]]; X <- cbind(1, as.matrix(d[, X_names]))
  id <- match(d[[subj]], unique(d[[subj]])); n_i <- tabulate(id)
  w <- d[[wt]][!duplicated(id)]; ly <- log(y); k <- ncol(X)
  nll <- function(p) {
    b <- p[1:k]; se2 <- exp(2 * p[k + 1]); su2 <- exp(2 * p[k + 2]); l <- p[k + 3]
    r  <- as.numeric(bc(y, l) - X %*% b)
    S1 <- rowsum(r, id)[, 1]; S2 <- rowsum(r^2, id)[, 1]
    J  <- rowsum((l - 1) * ly, id)[, 1]
    ll <- -0.5 * n_i * log(2 * pi) -
      0.5 * ((n_i - 1) * log(se2) + log(se2 + n_i * su2)) -
      0.5 / se2 * (S2 - su2 * S1^2 / (se2 + n_i * su2)) + J
    -sum(w * ll)
  }
  z0 <- bc(y, lambda_start); b0 <- qr.coef(qr(X), z0)
  o <- stats::optim(c(b0, log(stats::sd(z0)), log(stats::sd(z0)), lambda_start),
    nll, method = "L-BFGS-B",
    lower = c(rep(-Inf, k), -10, -10, 0.01), upper = c(rep(Inf, k), 10, 10, 2),
    control = list(maxit = 5000, factr = 1e3))
  list(beta = stats::setNames(o$par[1:k], c("(Intercept)", X_names)),
       lambda = o$par[k + 3], var_u2 = exp(2 * o$par[k + 2]),
       var_e = exp(2 * o$par[k + 1]), convergence = o$convergence)
}

# DISTRIB 9-point back-transformation over within-person error, floored at
# half the minimum amount (distrib_macro_v2.1.sas, macro MC).
cj <- c(-2.1, -1.3, -0.8, -0.5, 0, 0.5, 0.8, 1.3, 2.1)
wj <- c(0.063345, 0.080255, 0.070458, 0.159698, 0.252489,
        0.159698, 0.070458, 0.080255, 0.063345)
bt9 <- function(mu, l, ve, min_amt) {
  out <- if (l == 0) exp(mu + ve / 2) else
    Reduce(`+`, lapply(1:9, function(j)
      wj[j] * pmax((mu + cj[j] * sqrt(ve)) * l + 1, 0)^(1 / l)))
  pmax(out, 0.5 * min_amt)
}

# DISTRIB PC macro percentiles: midpoint plotting positions + interpolation
sas_pctl <- function(x, w, probs) {
  o <- order(x); x <- x[o]; w <- w[o] / sum(w)
  ux <- unique(x); wt <- rowsum(w, match(x, ux))[, 1]
  stats::approx(cumsum(wt) - wt / 2, ux, xout = probs, rule = 2,
                ties = "ordered")$y
}

# --- 2. MIXTRAN parameter parity ---------------------------------------------
sas_par <- c("(Intercept)" = 6.2277, age4 = 0.6813, age5 = 0.2975,
             age6 = -0.2031, age7 = -1.5547, age8 = -2.0490, race1 = -0.1547,
             race3 = 0.04750, weekend = 0.1367, day2 = -0.5416,
             lambda = 0.3941, var_u2 = 3.5177, var_e = 4.0771)  # example1_amount.lst

ref <- fit_amount_sas(m, "add_sug", "seqn", c(covs, "weekend", "day2"), "rndw1")
ref_par <- c(ref$beta, lambda = ref$lambda, var_u2 = ref$var_u2, var_e = ref$var_e)

pkg_fit <- mixtran(m, "add_sug", "seqn", "day", model_type = "amount",
                   covariates = covs, weekend_var = "weekend", verbose = FALSE)
pkg_par <- c(pkg_fit$beta[c("(Intercept)", covs, "weekend")],
             day2 = unname(pkg_fit$beta["seq_num"]), lambda = pkg_fit$lambda,
             var_u2 = pkg_fit$sigma2_b, var_e = pkg_fit$sigma2_w)

cat("\n== MIXTRAN parameters: SAS vs SAS-faithful R vs nciusual (unweighted) ==\n")
print(round(cbind(SAS = sas_par, faithful_R = ref_par[names(sas_par)],
                  nciusual = pkg_par[names(sas_par)]), 4))

pkg_w <- tryCatch(
  mixtran(m, "add_sug", "seqn", "day", model_type = "amount", covariates = covs,
          weekend_var = "weekend", weight_var = "rndw1", verbose = FALSE),
  error = function(e) conditionMessage(e))
cat("\nnciusual::mixtran(..., weight_var = 'rndw1'):\n  ",
    if (is.character(pkg_w)) paste("ERROR:", pkg_w) else "fitted", "\n")

# --- 3. DISTRIB parity using SAS's own parameters ----------------------------
b <- sas_par
min_amt <- min(m$add_sug[m$add_sug > 0])
s <- m[!duplicated(m$seqn), ]; s <- s[order(s$seqn), ]
xb0 <- as.numeric(b["(Intercept)"] + as.matrix(s[, covs]) %*% b[covs])  # seq = 0, weekday
xb1 <- xb0 + b[["weekend"]]                                             # weekend day

sas_out <- data.frame(
  grp  = c("Overall", 3:8),
  mean = c(25.4684, 28.1862, 33.1728, 30.1548, 26.6445, 18.4897, 15.9576),
  p1   = c(3.78757, 5.72069, 7.75554, 6.56467, 5.29683, 2.70647, 1.97508),
  p5   = c(7.2260, 9.7743, 12.5922, 10.9022, 9.0729, 5.1724, 4.1127),
  p10  = c(9.7658, 12.5933, 15.7673, 13.7782, 11.6484, 6.9446, 5.6476),
  p25  = c(15.2218, 18.0439, 22.2254, 19.6931, 17.0356, 10.7992, 9.0366),
  p50  = c(23.1779, 26.0912, 31.0057, 28.0696, 24.6893, 16.6851, 14.2047),
  p75  = c(33.2551, 36.0068, 41.8198, 38.3139, 34.0706, 24.1118, 20.9476),
  p90  = c(44.0691, 46.6128, 53.1966, 49.2491, 44.1306, 32.6036, 28.5617),
  p95  = c(51.4586, 53.9442, 61.1982, 56.7357, 51.0571, 38.0720, 33.8313),
  p99  = c(67.1537, 69.1836, 77.2138, 72.7663, 66.1512, 50.6311, 45.4636),
  below10 = c(0.10535, 0.05282, 0.02387, 0.03859, 0.06523, 0.21546, 0.29901))

nsim <- 100; set.seed(5454768)
u <- matrix(stats::rnorm(nrow(s) * nsim, 0, sqrt(b[["var_u2"]])), nrow(s))
lam <- b[["lambda"]]; ve <- b[["var_e"]]

ui_faithful <- (4 / 7) * bt9(xb0 + u, lam, ve, min_amt) +
               (3 / 7) * bt9(xb1 + u, lam, ve, min_amt)
ui_nciusual <- pmax(boxcox_inverse(xb0 + u, lam), 0)  # simulate_amount_model()

probs <- c(1, 5, 10, 25, 50, 75, 90, 95, 99) / 100
summarise <- function(UI) {
  rows <- lapply(sas_out$grp, function(g) {
    idx <- if (g == "Overall") seq_len(nrow(s)) else which(s$agegrp == as.numeric(g))
    x <- as.numeric(UI[idx, ]); w <- rep(s$rndw1[idx], nsim)
    c(mean = sum(x * w) / sum(w), stats::setNames(sas_pctl(x, w, probs), paste0("p", probs * 100)),
      below10 = sum(w[x < 10]) / sum(w))
  })
  as.data.frame(do.call(rbind, rows))
}
pct_diff <- function(R) cbind(grp = sas_out$grp,
  round(100 * (as.matrix(R) / as.matrix(sas_out[, -1]) - 1), 1))

cat("\n== (a) SAS-faithful DISTRIB vs SAS .lst, % difference ==\n")
print(pct_diff(summarise(ui_faithful)), quote = FALSE)
cat("\n== (b) nciusual DISTRIB logic with the SAME parameters vs SAS .lst, % difference ==\n")
print(pct_diff(summarise(ui_nciusual)), quote = FALSE)
