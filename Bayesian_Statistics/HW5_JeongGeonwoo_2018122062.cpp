#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Function to sample from a multivariate normal distribution
// [[Rcpp::export]]
arma::vec rmvnorm(const arma::vec& mean, const arma::mat& sigma) {
  int ncols = sigma.n_cols;
  arma::vec Y = arma::randn(ncols);
  return mean + arma::chol(sigma) * Y;
}

// Function to sample from an inverse gamma distribution
// [[Rcpp::export]]
double rinvgamma(double shape, double scale) {
  double d = R::rgamma(shape, 1.0 / scale);
  return 1.0 / d;
}

// The main Gibbs sampler function
// [[Rcpp::export]]
List RcppGibbs(int B, const arma::mat& X_design, const arma::vec& Y, double m_beta, double v_beta, double a_s2, double b_s2) {
  int p = X_design.n_cols;
  arma::mat beta_samps(p, B);
  arma::vec s2_samps(B);

  // Initial values
  beta_samps.col(0).fill(m_beta);
  s2_samps(0) = 1;

  // Gibbs sampling loop
  for(int i = 1; i < B; ++i) {
    // Sample beta given s2
    arma::mat V_beta = inv(X_design.t() * X_design / s2_samps(i - 1) + diagmat(arma::vec(p, arma::fill::value(v_beta))));
    arma::vec m_beta = V_beta * X_design.t() * Y / s2_samps(i - 1);
    beta_samps.col(i) = rmvnorm(m_beta, V_beta);

    // Sample s2 given beta
    double a_s2_new = a_s2 + Y.n_elem / 2.0;
    double residuals = as_scalar((Y - X_design * beta_samps.col(i)).t() * (Y - X_design * beta_samps.col(i)));
    double b_s2_new = b_s2 + residuals / 2.0;
    s2_samps(i) = rinvgamma(a_s2_new, b_s2_new);
  }

  return List::create(Named("beta_samps") = beta_samps, Named("s2_samps") = s2_samps);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List RcppGibbsparallel(int B, const arma::mat& X_design, const arma::vec& Y, double m_beta, double v_beta, double a_s2, double b_s2, int core) {
  int p = X_design.n_cols;
  cube beta_samps(p, B, core);
  mat s2_samps(B, core);

  // Set up parallel environment
  // omp_set_num_threads(core);

  // Parallel Gibbs sampling loop
  #pragma omp parallel for
  for (int j = 0; j < core; j++) {
    // Initial values for each chain
    beta_samps.slice(j).col(0).fill(m_beta);
    s2_samps(0, j) = 1;

    // Gibbs sampling loop for each chain
    for (int i = 1; i < B; ++i) {
      // Sample beta given s2
      arma::mat V_beta = inv(X_design.t() * X_design / s2_samps(i - 1, j) + diagmat(arma::vec(p, arma::fill::value(v_beta))));
      arma::vec m_beta = V_beta * X_design.t() * Y / s2_samps(i - 1, j);
      beta_samps.slice(j).col(i) = rmvnorm(m_beta, V_beta);

      // Sample s2 given beta
      double a_s2_new = a_s2 + Y.n_elem / 2.0;
      double residuals = as_scalar((Y - X_design * beta_samps.slice(j).col(i)).t() * (Y - X_design * beta_samps.slice(j).col(i)));
      double b_s2_new = b_s2 + residuals / 2.0;
      s2_samps(i, j) = rinvgamma(a_s2_new, b_s2_new);
    }
  }

  return List::create(Named("beta_samps") = beta_samps, Named("s2_samps") = s2_samps);
}










