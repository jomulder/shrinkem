// shrinkem_samplers.cpp
// RcppArmadillo ports of normal.horseshoe / normal.lasso / normal.ridge.
// Drop into src/. Requires DESCRIPTION: LinkingTo: Rcpp, RcppArmadillo.
//
// Uses R's RNG stream throughout (R::rgamma / R::norm_rand / R::unif_rand),
// so set.seed() works -- but draws will NOT be bit-identical to the R chain
// (different RNG call order). Validate by posterior summaries, not by seed.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// X ~ InvGamma(shape, scale): density propto x^{-shape-1} exp(-scale/x)
// (matches extraDistr::rinvgamma's (alpha, beta) = (shape, scale))
static inline double rinvgamma1(double shape, double scale) {
  return 1.0 / R::rgamma(shape, 1.0 / scale);   // rgamma takes (shape, scale)
}

// X ~ InverseGaussian(mu, lambda); Michael, Schucany & Haas (1976)
static inline double rinvgauss1(double mu, double lambda) {
  double nu = R::norm_rand();
  double y  = nu * nu;
  double x  = mu + (mu * mu * y) / (2.0 * lambda)
                 - (mu / (2.0 * lambda)) *
                   std::sqrt(4.0 * mu * lambda * y + mu * mu * y * y);
  double u  = R::unif_rand();
  return (u <= mu / (mu + x)) ? x : (mu * mu) / x;
}

// beta ~ N(Q^{-1} b, Q^{-1}),  Q = covInv + diag(dinv),  b = covInv * estimate
static inline vec sample_beta(const mat& covInv, const vec& dinv, const vec& b) {
  uword P = b.n_elem;
  mat Q = covInv;
  Q.diag() += dinv;
  mat R  = chol(Q);                       // upper, R.t() * R = Q
  vec y  = solve(trimatl(R.t()), b);      // R.t() y = b
  vec mu = solve(trimatu(R),   y);        // R mu  = y   ->  mu = Q^{-1} b
  vec z(P);
  for (uword i = 0; i < P; ++i) z[i] = R::norm_rand();
  return mu + solve(trimatu(R), z);       // Cov(R^{-1} z) = Q^{-1}
}

// build groups: members[g] = column indices with group_idx == g
static std::vector<std::vector<uword>> build_groups(const ivec& group_idx, int numGroup) {
  std::vector<std::vector<uword>> members(numGroup);
  for (uword j = 0; j < group_idx.n_elem; ++j)
    members[group_idx[j]].push_back(j);
  return members;
}

// ---------------------------------------------------------------- HORSESHOE
// [[Rcpp::export]]
List normal_horseshoe_cpp(const arma::vec& estimate, const arma::mat& covmatrix,
                          const arma::ivec& group_idx, int numGroup,
                          int iterations, int burnin, int store,
                          double a1, double a2, double b1,
                          double a3, double a4, double b2,
                          bool lambda2_fixed, const arma::vec& lambda2_input,
                          double nugget) {
  uword P = estimate.n_elem;
  mat covInv = inv_sympd(covmatrix);
  vec b = covInv * estimate;
  auto members = build_groups(group_idx, numGroup);
  vec lenGroup(numGroup, fill::zeros);
  for (int g = 0; g < numGroup; ++g) lenGroup[g] = (double) members[g].size();

  vec lambda2 = lambda2_fixed ? lambda2_input : vec(numGroup, fill::ones);
  vec gamma2  = lambda2_fixed ? lambda2_input : vec(numGroup, fill::ones);
  vec tau2(P, fill::ones), psi2(P, fill::ones), beta = estimate;
  vec lambda2vec(P);
  for (uword j = 0; j < P; ++j) lambda2vec[j] = lambda2[group_idx[j]];

  int nstore = iterations / store, sc = 0, total = burnin + iterations;
  mat beta_S(nstore, P), tau2_S(nstore, P), psi2_S(nstore, P);
  mat lam_S(nstore, numGroup), gam_S(nstore, numGroup);

  for (int t = 0; t < total; ++t) {
    beta = sample_beta(covInv, 1.0 / (tau2 % lambda2vec), b);

    for (uword j = 0; j < P; ++j) {
      double scale = psi2[j] + beta[j] * beta[j] / (2.0 * lambda2vec[j]);
      tau2[j] = rinvgamma1(a3 + 0.5, scale) + nugget;
    }
    for (uword j = 0; j < P; ++j)
      psi2[j] = R::rgamma(a3 + a4, 1.0 / (1.0 / b2 + 1.0 / tau2[j]));

    if (!lambda2_fixed) {
      for (int g = 0; g < numGroup; ++g) {
        double ss = gamma2[g];
        for (uword k : members[g]) ss += beta[k] * beta[k] / (2.0 * tau2[k]);
        lambda2[g] = rinvgamma1(a1 + lenGroup[g] / 2.0, ss);
        gamma2[g]  = R::rgamma(a1 + a2, 1.0 / (1.0 / b1 + 1.0 / lambda2[g]));
      }
    }
    for (uword j = 0; j < P; ++j) lambda2vec[j] = lambda2[group_idx[j]] + nugget;

    if (t >= burnin && ((t - burnin + 1) % store == 0)) {
      beta_S.row(sc) = beta.t(); tau2_S.row(sc) = tau2.t();
      psi2_S.row(sc) = psi2.t(); lam_S.row(sc) = lambda2.t();
      gam_S.row(sc) = gamma2.t(); ++sc;
    }
  }
  return List::create(_["beta"] = beta_S, _["tau2"] = tau2_S,
                      _["gamma2"] = gam_S, _["psi2"] = psi2_S,
                      _["lambda2"] = lam_S);
}

// -------------------------------------------------------------------- LASSO
// [[Rcpp::export]]
List normal_lasso_cpp(const arma::vec& estimate, const arma::mat& covmatrix,
                      const arma::ivec& group_idx, int numGroup,
                      int iterations, int burnin, int store,
                      double a1, double a2, double b1,
                      bool lambda2_fixed, const arma::vec& lambda2_input,
                      double nugget) {
  uword P = estimate.n_elem;
  mat covInv = inv_sympd(covmatrix);
  vec b = covInv * estimate;
  auto members = build_groups(group_idx, numGroup);
  vec lenGroup(numGroup, fill::zeros);
  for (int g = 0; g < numGroup; ++g) lenGroup[g] = (double) members[g].size();

  vec lambda2 = lambda2_fixed ? lambda2_input : vec(numGroup, fill::ones);
  vec gamma2  = lambda2_fixed ? lambda2_input : vec(numGroup, fill::ones);
  vec tau2(P, fill::ones), beta = estimate;
  vec lambda2vec(P);
  for (uword j = 0; j < P; ++j) lambda2vec[j] = lambda2[group_idx[j]];

  int nstore = iterations / store, sc = 0, total = burnin + iterations;
  mat beta_S(nstore, P), tau2_S(nstore, P);
  mat lam_S(nstore, numGroup), gam_S(nstore, numGroup);

  for (int t = 0; t < total; ++t) {
    beta = sample_beta(covInv, 1.0 / (tau2 % lambda2vec), b);

    for (uword j = 0; j < P; ++j) {
      double mu_tau = std::sqrt(lambda2vec[j] / (beta[j] * beta[j]));
      tau2[j] = 1.0 / rinvgauss1(mu_tau, 1.0) + nugget;
    }
    if (!lambda2_fixed) {
      for (int g = 0; g < numGroup; ++g) {
        double ss = gamma2[g];
        for (uword k : members[g]) ss += beta[k] * beta[k] / (2.0 * tau2[k]);
        lambda2[g] = rinvgamma1(a1 + lenGroup[g] / 2.0, ss);
        gamma2[g]  = R::rgamma(a1 + a2, 1.0 / (1.0 / b1 + 1.0 / lambda2[g]));
      }
    }
    for (uword j = 0; j < P; ++j) lambda2vec[j] = lambda2[group_idx[j]] + nugget;

    if (t >= burnin && ((t - burnin + 1) % store == 0)) {
      beta_S.row(sc) = beta.t(); tau2_S.row(sc) = tau2.t();
      lam_S.row(sc) = lambda2.t(); gam_S.row(sc) = gamma2.t(); ++sc;
    }
  }
  return List::create(_["beta"] = beta_S, _["lambda2"] = lam_S,
                      _["gamma2"] = gam_S, _["tau2"] = tau2_S);
}

// -------------------------------------------------------------------- RIDGE
// [[Rcpp::export]]
List normal_ridge_cpp(const arma::vec& estimate, const arma::mat& covmatrix,
                      const arma::ivec& group_idx, int numGroup,
                      int iterations, int burnin, int store,
                      double a1, double a2, double b1,
                      bool lambda2_fixed, const arma::vec& lambda2_input,
                      double nugget) {
  uword P = estimate.n_elem;
  mat covInv = inv_sympd(covmatrix);
  vec b = covInv * estimate;
  auto members = build_groups(group_idx, numGroup);
  vec lenGroup(numGroup, fill::zeros);
  for (int g = 0; g < numGroup; ++g) lenGroup[g] = (double) members[g].size();

  // fixed lambda2 -> posterior is exactly Gaussian; draw i.i.d. (no chain)
  if (lambda2_fixed) {
    vec lambda2vec(P);
    for (uword j = 0; j < P; ++j) lambda2vec[j] = lambda2_input[group_idx[j]];
    vec dinv = 1.0 / lambda2vec;
    mat beta_S(iterations, P);
    for (int t = 0; t < iterations; ++t) beta_S.row(t) = sample_beta(covInv, dinv, b).t();
    return List::create(_["beta"] = beta_S, _["lambda2"] = R_NilValue,
                        _["gamma2"] = R_NilValue);
  }

  vec lambda2(numGroup, fill::ones), gamma2(numGroup, fill::ones), beta = estimate;
  vec lambda2vec(P);
  for (uword j = 0; j < P; ++j) lambda2vec[j] = lambda2[group_idx[j]];

  int nstore = iterations / store, sc = 0, total = burnin + iterations;
  mat beta_S(nstore, P), lam_S(nstore, numGroup), gam_S(nstore, numGroup);

  for (int t = 0; t < total; ++t) {
    beta = sample_beta(covInv, 1.0 / lambda2vec, b);
    for (int g = 0; g < numGroup; ++g) {
      double ss = gamma2[g];
      for (uword k : members[g]) ss += beta[k] * beta[k] / 2.0;
      lambda2[g] = rinvgamma1(a1 + lenGroup[g] / 2.0, ss);
      gamma2[g]  = R::rgamma(a1 + a2, 1.0 / (1.0 / b1 + 1.0 / lambda2[g]));
    }
    for (uword j = 0; j < P; ++j) lambda2vec[j] = lambda2[group_idx[j]] + nugget;

    if (t >= burnin && ((t - burnin + 1) % store == 0)) {
      beta_S.row(sc) = beta.t(); lam_S.row(sc) = lambda2.t();
      gam_S.row(sc) = gamma2.t(); ++sc;
    }
  }
  return List::create(_["beta"] = beta_S, _["lambda2"] = lam_S, _["gamma2"] = gam_S);
}
