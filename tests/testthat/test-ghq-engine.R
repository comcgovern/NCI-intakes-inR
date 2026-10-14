# Tests for the GHQ (Gauss-Hermite Quadrature) engine of the correlated
# two-part model, and the quadrature utility functions.

# ---------------------------------------------------------------------------
# Quadrature utility tests
# ---------------------------------------------------------------------------

test_that("gh_nodes_weights returns correct structure", {
  for (n in c(3, 5, 7, 9)) {
    gh <- gh_nodes_weights(n)
    expect_length(gh$nodes, n)
    expect_length(gh$weights, n)
    expect_true(all(gh$weights > 0))
    # Weights should sum to sqrt(pi) for physicist convention
    expect_equal(sum(gh$weights), sqrt(pi), tolerance = 1e-3)
    # Nodes should be symmetric about zero
    expect_equal(gh$nodes, -rev(gh$nodes), tolerance = 1e-10)
  }
})

test_that("gh_nodes_normal weights sum to 1", {
  for (n in c(3, 5, 7)) {
    ghn <- gh_nodes_normal(n)
    expect_equal(sum(ghn$weights), 1, tolerance = 1e-10)
    expect_true(all(ghn$weights > 0))
  }
})

test_that("gh_nodes_normal integrates x^2 against N(0,1) ≈ 1", {
  # E[X^2] = 1 for X ~ N(0,1)
  ghn <- gh_nodes_normal(7)
  integral <- sum(ghn$weights * ghn$nodes^2)
  expect_equal(integral, 1, tolerance = 0.01)
})

test_that("gh_nodes_bivariate has correct dimensions", {
  biv <- gh_nodes_bivariate(5, rho = 0.3, sigma_v1 = 1, sigma_v2 = 1.2)
  expect_equal(nrow(biv), 25)
  expect_true(all(c("v1", "v2", "log_weight") %in% names(biv)))
})

test_that("gh_nodes_bivariate marginal variances match inputs", {
  # With rho=0, v1 and v2 should be independent
  # E[v1^2] = sigma_v1^2, E[v2^2] = sigma_v2^2
  biv <- gh_nodes_bivariate(9, rho = 0, sigma_v1 = 2, sigma_v2 = 3)
  w   <- exp(biv$log_weight)
  ev1_sq <- sum(w * biv$v1^2)
  ev2_sq <- sum(w * biv$v2^2)
  expect_equal(ev1_sq, 4, tolerance = 0.05)
  expect_equal(ev2_sq, 9, tolerance = 0.05)
})

test_that("gh_nodes_bivariate covariance matches rho * sigma_v1 * sigma_v2", {
  rho <- 0.5; sv1 <- 1.5; sv2 <- 2.0
  biv <- gh_nodes_bivariate(9, rho = rho, sigma_v1 = sv1, sigma_v2 = sv2)
  w   <- exp(biv$log_weight)
  cov_v1v2 <- sum(w * biv$v1 * biv$v2)
  expect_equal(cov_v1v2, rho * sv1 * sv2, tolerance = 0.05)
})

test_that("gh_nodes_weights works for non-table n via eigenvalue method", {
  gh <- gh_nodes_weights(4)
  expect_length(gh$nodes, 4)
  expect_equal(sum(gh$weights), sqrt(pi), tolerance = 1e-3)
})


test_that("tabulated GH nodes match the Golub-Welsch computation", {
  gw <- function(n) {
    off <- sqrt(seq_len(n - 1) / 2); J <- diag(0, n)
    J[cbind(1:(n - 1), 2:n)] <- off; J[cbind(2:n, 1:(n - 1))] <- off
    e <- eigen(J, symmetric = TRUE); o <- order(e$values)
    list(nodes = e$values[o], weights = (e$vectors[1, ]^2 * sqrt(pi))[o])
  }
  for (n in c(3, 5, 7, 9)) {
    expect_equal(gh_nodes_weights(n)$nodes, gw(n)$nodes, tolerance = 1e-10)
    expect_equal(gh_nodes_weights(n)$weights, gw(n)$weights, tolerance = 1e-10)
  }
})

test_that("gh_nodes_adaptive weights give prior expectations", {
  # With no data the posterior is the prior: weights must integrate 1 and v1^2
  a <- gh_nodes_adaptive(9, eta_i = numeric(0), mu_i = numeric(0),
                         consumed_i = numeric(0), t_y_i = numeric(0),
                         sigma_v1 = 2, sigma_v2 = 0.5, sigma_e = 1, rho = 0.3,
                         lambda = 0.5)
  expect_equal(sum(exp(a$log_weight)), 1, tolerance = 1e-6)
  expect_equal(sum(exp(a$log_weight) * a$v1^2), 4, tolerance = 1e-4)
  expect_equal(sum(exp(a$log_weight) * a$v1 * a$v2), 0.3 * 2 * 0.5, tolerance = 1e-4)
})

test_that("adaptive quadrature in the engine matches integrate() for the probability part", {
  set.seed(8)
  n <- 40
  dat <- data.frame(id = rep(1:n, each = 3), day = rep(1:3, n))
  v <- rep(rnorm(n, 0, 1.5), each = 3)
  dat$y <- rbinom(3 * n, 1, plogis(0.4 + v)) * exp(rnorm(3 * n, 1, 0.5))
  fit <- suppressMessages(mixtran(dat, "y", "id", "day", model_type = "uncorr",
                                  lambda = 0, include_seq = FALSE, verbose = FALSE))
  d <- nciusual:::.nci_design(fit$data, fit$prob_covs_used, fit$amt_covs_used, 0, 21, has_prob = TRUE)
  eng <- nciusual:::.ll_subject(fit$theta, d, "prob")
  sd1 <- sqrt(fit$sigma2_v1); a0 <- fit$alpha[["(Intercept)"]]
  bf <- vapply(1:n, function(i) {
    yn <- fit$data$consumed[fit$data$subject == i]
    log(integrate(function(u) vapply(u, function(uu)
      prod(dbinom(yn, 1, plogis(a0 + uu))) * dnorm(uu, 0, sd1), 1), -Inf, Inf,
      rel.tol = 1e-12)$value)
  }, 1)
  expect_equal(unname(eng), bf, tolerance = 1e-7)
})

test_that("distrib() works after a correlated fit", {
  set.seed(42)
  n <- 400
  nr <- sample(1:2, n, replace = TRUE, prob = c(0.35, 0.65))
  sid <- rep(seq_len(n), nr); day <- unlist(lapply(nr, seq_len))
  z1 <- rnorm(n); z2 <- rnorm(n)
  v1 <- rep(z1, nr); v2 <- rep(0.9 * (0.5 * z1 + sqrt(0.75) * z2), nr)
  y <- rbinom(length(sid), 1, plogis(0.2 + v1)) *
    pmax(boxcox_inverse(2.8 + v2 + rnorm(length(sid), 0, sqrt(0.7)), 0.3), 0.1)
  fit <- suppressMessages(mixtran(data.frame(id = sid, day = day, y = y), "y", "id", "day",
                                  model_type = "corr", lambda = 0.3, verbose = FALSE))
  res <- distrib(fit, n_sims = 50, seed = 1)
  expect_s3_class(res, "distrib_result")
  expect_true(is.finite(res$results$Overall$mean))
  expect_true(all(diff(res$results$Overall$percentiles) >= 0))
})
