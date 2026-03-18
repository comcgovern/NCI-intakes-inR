#' Summary and Plot Methods for nciusual Objects
#'
#' summary() and plot() methods for `mixtran_fit` and `distrib_result` objects.
#'
#' @name print_methods
NULL


# ============================================================
# mixtran_fit
# ============================================================

#' Summarise a fitted MIXTRAN model
#'
#' @param object A `mixtran_fit` object from mixtran()
#' @param ... Unused
#' @return Invisibly returns a list with formatted summary components
#' @export
summary.mixtran_fit <- function(object, ...) {
  x <- object

  cat("═══════════════════════════════════════════════════════\n")
  cat("  NCI Usual Intake Model — MIXTRAN\n")
  cat("═══════════════════════════════════════════════════════\n")
  cat(sprintf("  Intake variable : %s\n", x$intake_var))
  cat(sprintf("  Model type      : %s\n", x$model_type))
  cat(sprintf("  Converged       : %s\n", x$converged))
  cat(sprintf("  Box-Cox lambda  : %.4f\n", x$lambda))
  cat(sprintf("  N subjects      : %d\n", x$n_subjects))
  cat(sprintf("  N recalls       : %d\n", x$n_recalls))
  cat(sprintf("  N positive      : %d  (%.1f%% of recalls)\n",
              x$n_positive, x$pct_consuming))

  if (!is.null(x$covariates)) {
    cat(sprintf("  Covariates      : %s\n", paste(x$covariates, collapse = ", ")))
  }
  if (!is.null(x$weekend_var)) {
    cat(sprintf("  Weekend var     : %s\n", x$weekend_var))
  }

  if (x$model_type == "amount") {
    cat("\n  Variance Components\n")
    cat(sprintf("    Between-person (σ²_b) : %10.4f\n", x$sigma2_b))
    cat(sprintf("    Within-person  (σ²_w) : %10.4f\n", x$sigma2_w))
    cat(sprintf("    Ratio (w/b)           : %10.4f\n", x$var_ratio))
    cat("\n  Fixed Effects (amount model)\n")
    for (nm in names(x$beta)) {
      cat(sprintf("    %-22s: %10.4f\n", nm, x$beta[nm]))
    }
  } else {
    cat("\n  Variance Components\n")
    cat(sprintf("    σ²_v1 (between, prob)  : %10.4f\n", x$sigma2_v1))
    cat(sprintf("    σ²_v2 (between, amount): %10.4f\n", x$sigma2_v2))
    cat(sprintf("    σ²_e  (within, amount) : %10.4f\n", x$sigma2_e))
    cat(sprintf("    ρ (v1-v2 correlation)  : %10.4f\n", x$rho))
    cat("\n  Fixed Effects (probability model)\n")
    for (nm in names(x$alpha)) {
      cat(sprintf("    %-22s: %10.4f\n", nm, x$alpha[nm]))
    }
    cat("\n  Fixed Effects (amount model)\n")
    for (nm in names(x$beta)) {
      cat(sprintf("    %-22s: %10.4f\n", nm, x$beta[nm]))
    }
  }
  cat("═══════════════════════════════════════════════════════\n")
  invisible(x)
}


#' Plot diagnostics for a fitted MIXTRAN model
#'
#' Produces up to four diagnostic plots:
#'   1. Histogram of transformed intakes with fitted normal overlay
#'   2. Q-Q plot of transformed residuals
#'   3. Profile likelihood of Box-Cox lambda (if estimated from data)
#'   4. Histogram of raw intakes for context
#'
#' @param x A `mixtran_fit` object from mixtran()
#' @param which Integer vector: which plots to draw (1–4). Default all.
#' @param ... Passed to plot()
#' @export
plot.mixtran_fit <- function(x, which = 1:4, ...) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))

  n_plots <- length(which)
  nc <- min(n_plots, 2L)
  nr <- ceiling(n_plots / nc)
  if (n_plots > 1) graphics::par(mfrow = c(nr, nc))

  pred <- x$predicted

  # --- Plot 1: Histogram of transformed intake ---
  if (1 %in% which) {
    t_col <- if ("t_intake" %in% names(pred)) "t_intake" else NULL
    if (!is.null(t_col)) {
      tz <- pred[[t_col]][!is.na(pred[[t_col]])]
      graphics::hist(tz, breaks = 40, freq = FALSE,
                     main = "Transformed intake (Box-Cox)",
                     xlab = paste0("T(intake), lambda=", round(x$lambda, 3)),
                     col = "steelblue", border = "white")
      xg <- seq(min(tz), max(tz), length.out = 200)
      if (x$model_type == "amount") {
        mu    <- mean(x$beta)
        sigma <- sqrt(x$sigma2_b + x$sigma2_w)
      } else {
        mu    <- mean(x$beta)
        sigma <- sqrt(x$sigma2_v2 + x$sigma2_e)
      }
      graphics::lines(xg, stats::dnorm(xg, mu, sigma),
                      col = "tomato", lwd = 2)
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Transformed data\nnot stored in predicted",
                     cex = 0.9)
    }
  }

  # --- Plot 2: Q-Q plot of residuals ---
  if (2 %in% which) {
    if (x$model_type == "amount" && !is.null(x$model_fit)) {
      res <- stats::residuals(x$model_fit, type = "normalized")
      stats::qqnorm(res, main = "Q-Q: Normalised residuals", pch = 16,
                    cex = 0.4, col = grDevices::rgb(0.2, 0.4, 0.7, 0.5))
      stats::qqline(res, col = "tomato", lwd = 2)
    } else if (!is.null(x$amt_fit)) {
      res <- stats::residuals(x$amt_fit, type = "normalized")
      stats::qqnorm(res, main = "Q-Q: Amount model residuals", pch = 16,
                    cex = 0.4, col = grDevices::rgb(0.2, 0.4, 0.7, 0.5))
      stats::qqline(res, col = "tomato", lwd = 2)
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Residuals not available", cex = 0.9)
    }
  }

  # --- Plot 3: Profile likelihood of rho (corr model) ---
  if (3 %in% which) {
    if (x$model_type == "corr" && !is.null(x$rho_profile)) {
      rp <- x$rho_profile
      graphics::plot(rp$rho, rp$loglik, type = "l",
                     main = expression("Profile likelihood: " * rho),
                     xlab = expression(rho), ylab = "Log-likelihood",
                     lwd = 2, col = "steelblue")
      graphics::abline(v = x$rho, col = "tomato", lty = 2, lwd = 2)
      graphics::legend("topright",
                       legend = sprintf("hat(rho) = %.3f", x$rho),
                       col = "tomato", lty = 2, bty = "n")
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "rho profile only available\nfor corr model", cex = 0.9)
    }
  }

  # --- Plot 4: Histogram of raw (back-transformed) intake ---
  if (4 %in% which) {
    raw_col <- x$intake_var
    if (!is.null(pred) && raw_col %in% names(pred)) {
      ry <- pred[[raw_col]]
      ry <- ry[!is.na(ry) & ry > 0]
      graphics::hist(ry, breaks = 50, freq = FALSE,
                     main = paste("Raw intake:", raw_col),
                     xlab = raw_col,
                     col = "mediumseagreen", border = "white")
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Raw intake data not stored\nin predicted", cex = 0.9)
    }
  }

  invisible(x)
}


# ============================================================
# distrib_result
# ============================================================

#' Summarise a DISTRIB usual intake result
#'
#' @param object A `distrib_result` object from distrib()
#' @param ... Unused
#' @return Invisibly returns the object
#' @export
summary.distrib_result <- function(object, ...) {
  x <- object
  cat("═══════════════════════════════════════════════════════\n")
  cat("  NCI Usual Intake Distribution — DISTRIB\n")
  cat("═══════════════════════════════════════════════════════\n")
  cat(sprintf("  Model type  : %s\n", x$mixtran_obj$model_type))
  cat(sprintf("  Intake var  : %s\n", x$mixtran_obj$intake_var))
  cat(sprintf("  Lambda      : %.4f\n", x$mixtran_obj$lambda))
  cat(sprintf("  N subjects  : %d\n", x$mixtran_obj$n_subjects))
  cat(sprintf("  N sims/subj : %d\n", x$n_sims))
  cat(sprintf("  Seed        : %d\n", x$seed))
  if (!is.null(x$subgroup_var)) {
    cat(sprintf("  Subgroup    : %s\n", x$subgroup_var))
  }
  cat("\n")

  for (grp in names(x$results)) {
    r <- x$results[[grp]]
    cat(sprintf("  ── %s ──\n", grp))
    cat(sprintf("    Mean (SD): %.2f  (%.2f)\n", r$mean, r$sd))
    cat("    Percentiles:\n")
    pct_vals <- r$percentiles
    cat(sprintf("      %s\n",
                paste(sprintf("%s=%.2f", names(pct_vals), pct_vals), collapse = "  ")))
    if (!is.null(r$cutpoint_below)) {
      cat("    Cutpoints (% below):\n")
      for (nm in names(r$cutpoint_below)) {
        cat(sprintf("      %s: %.1f%%\n", nm, r$cutpoint_below[[nm]] * 100))
      }
    }
    cat("\n")
  }
  cat("═══════════════════════════════════════════════════════\n")
  invisible(x)
}


#' Plot the estimated usual intake distribution
#'
#' Produces up to three plots:
#'   1. Density of simulated usual intakes (one curve per subgroup)
#'   2. Percentile bar chart comparing subgroups
#'   3. (If cutpoints provided) Bar chart of prevalence below each cutpoint
#'
#' @param x A `distrib_result` object from distrib()
#' @param which Integer vector: which plots to draw (1–3). Default 1:2.
#' @param ... Passed to plot()
#' @export
plot.distrib_result <- function(x, which = 1:2, ...) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))

  n_plots <- length(which)
  nc <- min(n_plots, 2L)
  nr <- ceiling(n_plots / nc)
  if (n_plots > 1) graphics::par(mfrow = c(nr, nc))

  grps <- names(x$results)
  n_grps <- length(grps)

  # Colour palette
  pal <- if (n_grps <= 8) {
    c("steelblue", "tomato", "mediumseagreen", "goldenrod",
      "mediumpurple", "coral", "teal", "sienna")[seq_len(n_grps)]
  } else {
    grDevices::rainbow(n_grps)
  }

  # --- Plot 1: Density of simulated usual intakes ---
  if (1 %in% which) {
    # Compute densities from the simulated values
    sim <- x$simulated
    n_subj <- sim$n_subjects
    n_s    <- sim$n_sims

    # For overall density we use all simulated values
    all_ui <- sim$usual_intake
    d_all  <- stats::density(all_ui[all_ui >= 0], na.rm = TRUE)

    graphics::plot(d_all,
                   main = paste("Usual intake distribution\n(", x$mixtran_obj$intake_var, ")"),
                   xlab = x$mixtran_obj$intake_var,
                   ylab = "Density",
                   col = pal[1], lwd = 2,
                   xlim = range(d_all$x))
    graphics::legend("topright", legend = grps[1], col = pal[1],
                     lwd = 2, bty = "n")
  }

  # --- Plot 2: Percentile bar chart ---
  if (2 %in% which) {
    pct_names <- names(x$results[[1]]$percentiles)
    pct_matrix <- do.call(cbind, lapply(grps, function(g) {
      x$results[[g]]$percentiles
    }))
    colnames(pct_matrix) <- grps

    graphics::barplot(
      t(pct_matrix),
      beside   = TRUE,
      names.arg = pct_names,
      col      = pal,
      border   = NA,
      main     = paste("Percentiles of usual intake\n(", x$mixtran_obj$intake_var, ")"),
      xlab     = "Percentile",
      ylab     = x$mixtran_obj$intake_var,
      legend.text = grps,
      args.legend = list(bty = "n", x = "topleft")
    )
  }

  # --- Plot 3: Cutpoint prevalences ---
  if (3 %in% which) {
    first_cp <- x$results[[grps[1]]]$cutpoint_below
    if (is.null(first_cp)) {
      graphics::plot.new()
      graphics::text(0.5, 0.5,
                     "No cutpoints specified.\nUse cutpoints= argument in distrib().",
                     cex = 0.9)
    } else {
      cp_names <- names(first_cp)
      cp_matrix <- do.call(cbind, lapply(grps, function(g) {
        x$results[[g]]$cutpoint_below * 100
      }))
      colnames(cp_matrix) <- grps

      graphics::barplot(
        t(cp_matrix),
        beside    = TRUE,
        names.arg = gsub("pct_below_", "", cp_names),
        col       = pal,
        border    = NA,
        main      = "Prevalence below cutpoints (%)",
        xlab      = "Cutpoint",
        ylab      = "% of population",
        legend.text = grps,
        args.legend = list(bty = "n", x = "topright")
      )
    }
  }

  invisible(x)
}
