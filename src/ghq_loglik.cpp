// ghq_loglik.cpp
//
// Rcpp-accelerated Gauss-Hermite quadrature log-likelihood for the NCI
// two-part correlated model.  This is the C++ equivalent of the inner loop
// in fit_twopart_corr_ghq() (R/mixtran.R).
//
// Compilation: add `Rcpp` to `LinkingTo` and `Imports` in DESCRIPTION,
// then call `Rcpp::compileAttributes()` and rebuild the package.
//
// Usage (once compiled):
//   ghq_loglik_cpp(
//     log_sv1, log_sv2, log_se, atanh_rho,
//     eta_pop, mu_pop, consumed, t_y_pos, log_y_pos, lambda,
//     person_id, person_id_pos, n_subj,
//     gh_nodes, gh_log_weights
//   )
//
// Speedup over pure-R: ~50-100× for large n_subj × n_q loops due to
// elimination of tapply() overhead and vectorised C++ inner loop.
//
// NOTE: This file is a Phase 2 stub.  The R fallback in R/mixtran.R is fully
// functional.  Activate this implementation by:
//   1. Adding `Rcpp` to Imports and LinkingTo in DESCRIPTION
//   2. Uncommenting the [[Rcpp::export]] attribute below
//   3. Running Rcpp::compileAttributes() and rebuilding

// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// Log-sum-exp for a vector (numerically stable)
static inline double log_sum_exp(const std::vector<double>& x) {
  double m = *std::max_element(x.begin(), x.end());
  double s = 0.0;
  for (double xi : x) s += std::exp(xi - m);
  return m + std::log(s);
}

// plogis(x) = 1 / (1 + exp(-x)), log version via log1p
static inline double log_plogis(double x) {
  return -std::log1p(std::exp(-x));
}
static inline double log_plogis_lower(double x) {
  return -std::log1p(std::exp(x));
}

//' Gauss-Hermite quadrature log-likelihood (C++ implementation)
//'
//' Computes the marginal log-likelihood for the NCI two-part correlated model
//' by integrating out (v1_i, v2_i) ~ BVN(0, Sigma) via tensor-product GHQ.
//'
//' @param log_sv1 log(sigma_v1)
//' @param log_sv2 log(sigma_v2)
//' @param log_se  log(sigma_e)
//' @param atanh_rho atanh(rho)
//' @param eta_pop Population-level probability linear predictor (all obs)
//' @param mu_pop Population-level amount linear predictor (positive obs)
//' @param consumed Integer vector: 1 = consumed, 0 = not (all obs)
//' @param t_y_pos Box-Cox transformed intake (positive obs)
//' @param log_y_pos log(y) for Jacobian (positive obs)
//' @param lambda Box-Cox lambda
//' @param person_id Integer subject index 1..n_subj (all obs)
//' @param person_id_pos Integer subject index 1..n_subj (positive obs)
//' @param n_subj Number of unique subjects
//' @param gh_nodes GH nodes on the N(0,1) scale (length n_nodes)
//' @param gh_log_weights Log of GH weights (sum-to-one, length n_nodes)
//' @return Scalar: total marginal log-likelihood
//'
//' @keywords internal
// [[Rcpp::export]]
double ghq_loglik_cpp(double log_sv1, double log_sv2, double log_se,
                      double atanh_rho,
                      NumericVector eta_pop,
                      NumericVector mu_pop,
                      IntegerVector consumed,
                      NumericVector t_y_pos,
                      NumericVector log_y_pos,
                      double lambda,
                      IntegerVector person_id,
                      IntegerVector person_id_pos,
                      int n_subj,
                      NumericVector gh_nodes,
                      NumericVector gh_log_weights) {

  double sv1 = std::exp(log_sv1);
  double sv2 = std::exp(log_sv2);
  double se  = std::exp(log_se);
  double rho = std::tanh(atanh_rho);
  double sqrt_1r2 = std::sqrt(std::max(0.0, 1.0 - rho * rho));

  int n_nodes = gh_nodes.size();
  int n_obs   = eta_pop.size();
  int n_pos   = t_y_pos.size();
  int n_q     = n_nodes * n_nodes;

  // ll_mat[subj][q]: log-likelihood contribution for subject subj at quad point q
  std::vector<std::vector<double>> ll_mat(n_subj,
                                          std::vector<double>(n_q, 0.0));
  std::vector<double> log_wt(n_q);

  int q = 0;
  for (int q1 = 0; q1 < n_nodes; ++q1) {
    for (int q2 = 0; q2 < n_nodes; ++q2, ++q) {
      double z1 = gh_nodes[q1];
      double z2 = gh_nodes[q2];
      double v1 = sv1 * z1;
      double v2 = sv2 * (rho * z1 + sqrt_1r2 * z2);
      log_wt[q]  = gh_log_weights[q1] + gh_log_weights[q2];

      // Initialise per-subject accumulator for this quadrature point
      std::vector<double> ll_subj(n_subj, 0.0);

      // Probability contributions (all observations)
      for (int n = 0; n < n_obs; ++n) {
        int subj = person_id[n] - 1;  // 0-based
        double eta_q = eta_pop[n] + v1;
        ll_subj[subj] += (consumed[n] == 1) ?
          log_plogis(eta_q) :
          log_plogis_lower(eta_q);
      }

      // Amount contributions (positive observations)
      for (int n = 0; n < n_pos; ++n) {
        int subj = person_id_pos[n] - 1;  // 0-based
        double mu_q  = mu_pop[n] + v2;
        double resid = t_y_pos[n] - mu_q;
        double ll_norm = -0.5 * std::log(2.0 * M_PI) - std::log(se)
                         - 0.5 * (resid * resid) / (se * se);
        double ll_jac  = (lambda - 1.0) * log_y_pos[n];
        ll_subj[subj] += ll_norm + ll_jac;
      }

      for (int i = 0; i < n_subj; ++i) {
        ll_mat[i][q] = ll_subj[i];
      }
    }
  }

  // Marginal log-likelihood: sum_i log(sum_q w_q exp(ll_i_q))
  double total_ll = 0.0;
  for (int i = 0; i < n_subj; ++i) {
    std::vector<double> ll_plus_lw(n_q);
    for (int qi = 0; qi < n_q; ++qi) {
      ll_plus_lw[qi] = ll_mat[i][qi] + log_wt[qi];
    }
    total_ll += log_sum_exp(ll_plus_lw);
  }

  return total_ll;
}
