// Moment-based approximations for the Wright-Fisher model of population dynamics under natural selection at two linked loci
// Zhangyi He, Wenyang Lyu, Mark Beaumont, Feng Yu

// version 1.0.1

// C functions

#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]
#include <math.h>

using namespace Rcpp;
using namespace std;
// using namespace arma;

/********** WFM **********/
// Calculate the fitness matrix for the Wright-Fisher model
// [[Rcpp::export]]
arma::dmat calculateFitnessMat_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the fitness matrix
  arma::dmat fts_mat = arma::ones<arma::dmat>(4, 4);
  // fts_mat(0, 0) = 1;
  fts_mat(1, 0) = (1 - dom_par_B * sel_cof_B);
  fts_mat(2, 0) = (1 - dom_par_A * sel_cof_A);
  fts_mat(3, 0) = (1 - dom_par_A * sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts_mat(0, 1) = fts_mat(1, 0);
  fts_mat(1, 1) = (1 - sel_cof_B);
  fts_mat(2, 1) = (1 - dom_par_A * sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts_mat(3, 1) = (1 - dom_par_A * sel_cof_A) * (1 - sel_cof_B);
  fts_mat(0, 2) = fts_mat(2, 0);
  fts_mat(1, 2) = fts_mat(2, 1);
  fts_mat(2, 2) = (1 - sel_cof_A);
  fts_mat(3, 2) = (1 - sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts_mat(0, 3) = fts_mat(3, 0);
  fts_mat(1, 3) = fts_mat(3, 1);
  fts_mat(2, 3) = fts_mat(3, 2);
  fts_mat(3, 3) = (1 - sel_cof_A) * (1 - sel_cof_B);

  // return the fitness matrix for the Wright-Fisher model
  return fts_mat;
}

// Simulate the haplotype frequency trajectories according to the Wright-Fisher model
// [[Rcpp::export]]
arma::dmat simulateWFM_arma(const arma::dmat& fts_mat, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the haplotype frequency trajectories
  arma::dmat frq_pth(4, arma::uword(lst_gen - int_gen) + 1);

  // initialise the haplotype frequencies in generation 0
  frq_pth.col(0) = int_frq;

  // declare eta
  arma::dcolvec eta = {-1.0, 1.0, 1.0, -1.0};

  // simulate the haplotype frequency trajectories
  arma::dcolvec hap_frq = int_frq;
  for (arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    // calculate the sampling probabilities
    arma::dcolvec prob = hap_frq;
    prob = hap_frq % (fts_mat * hap_frq) / arma::as_scalar(hap_frq.t() * fts_mat * hap_frq);
    prob = prob + eta * rec_rat * (prob(0) * prob(3) - prob(1) * prob(2));

    // proceed the Wright-Fisher sampling
    IntegerVector hap_cnt(4);
    R::rmultinom(2 * pop_siz, prob.begin(), 4, hap_cnt.begin());
    hap_frq = as<arma::dcolvec>(hap_cnt) / 2 / pop_siz;
    frq_pth.col(k) = hap_frq;
  }

  // return the haplotype frequency trajectories under the Wright-Fisher model
  return frq_pth;
}
/*************************/


/***** DiffusApprox ******/
// Simulate the haplotype frequency trajectories according to the Wright-Fisher diffusion using the Euler-Maruyama method
// [[Rcpp::export]]
arma::dmat simulateDiffusApprox_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // rescale the selection coefficients and the recombination rate
  double scl_sel_cof_A = 2 * pop_siz * sel_cof_A;
  double scl_sel_cof_B = 2 * pop_siz * sel_cof_B;
  double scl_rec_rat = 4 * pop_siz * rec_rat;

  // calculate delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  // generate delta W
  arma::dmat dW = pow(dt, 0.5) * arma::randn<arma::dmat>(6, arma::uword(lst_gen - int_gen) * ptn_num);

  // declare the haplotype frequency trajectories
  arma::dmat frq_pth(4, arma::uword(lst_gen - int_gen) * ptn_num + 1);

  // initialise the haplotype frequencies in generation 0
  frq_pth.col(0) = int_frq;

  // simulate the haplotype frequency trajectories
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) * ptn_num + 1; t++) {
    // calculate the drift coefficient vector
    arma::dcolvec mu = arma::zeros<arma::dcolvec>(4);
    mu(0) =  scl_sel_cof_A * frq_pth(0, t - 1) * (frq_pth(2, t - 1) + frq_pth(3, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(1, t - 1)) * dom_par_A + (frq_pth(2, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_A)) +
      scl_sel_cof_B * frq_pth(0, t - 1) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par_B + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_B)) -
      0.5 * scl_rec_rat * (frq_pth(0, t - 1) * frq_pth(3, t - 1) - frq_pth(1, t - 1) * frq_pth(2, t - 1));
    mu(1) =  scl_sel_cof_A * frq_pth(1, t - 1) * (frq_pth(2, t - 1) + frq_pth(3, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(1, t - 1)) * dom_par_A + (frq_pth(2, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_A)) -
      scl_sel_cof_B * frq_pth(1, t - 1) * (frq_pth(0, t - 1) + frq_pth(2, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par_B + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_B)) +
      0.5 * scl_rec_rat * (frq_pth(0, t - 1) * frq_pth(3, t - 1) - frq_pth(1, t - 1) * frq_pth(2, t - 1));
    mu(2) = -scl_sel_cof_A * frq_pth(2, t - 1) * (frq_pth(0, t - 1) + frq_pth(1, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(1, t - 1)) * dom_par_A + (frq_pth(2, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_A)) +
      scl_sel_cof_B * frq_pth(2, t - 1) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par_B + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_B)) +
      0.5 * scl_rec_rat * (frq_pth(0, t - 1) * frq_pth(3, t - 1) - frq_pth(1, t - 1) * frq_pth(2, t - 1));
    mu(3) = -scl_sel_cof_A * frq_pth(3, t - 1) * (frq_pth(0, t - 1) + frq_pth(1, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(1, t - 1)) * dom_par_A + (frq_pth(2, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_A)) -
      scl_sel_cof_B * frq_pth(3, t - 1) * (frq_pth(0, t - 1) + frq_pth(2, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par_B + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_B)) -
      0.5 * scl_rec_rat * (frq_pth(0, t - 1) * frq_pth(3, t - 1) - frq_pth(1, t - 1) * frq_pth(2, t - 1));

    // calculate the diffusion coefficient matrix
    arma::dmat sigma = arma::zeros<arma::dmat>(4, 6);
    sigma(0, 0) = pow(frq_pth(0, t - 1) * frq_pth(1, t - 1), 0.5);
    sigma(0, 1) = pow(frq_pth(0, t - 1) * frq_pth(2, t - 1), 0.5);
    sigma(0, 2) = pow(frq_pth(0, t - 1) * frq_pth(3, t - 1), 0.5);
    // sigma(0, 3) = 0;
    // sigma(0, 4) = 0;
    // sigma(0, 5) = 0;
    sigma(1, 0) = -pow(frq_pth(1, t - 1) * frq_pth(0, t - 1), 0.5);
    // sigma(1, 1) = 0;
    // sigma(1, 2) = 0;
    sigma(1, 3) = pow(frq_pth(1, t - 1) * frq_pth(2, t - 1), 0.5);
    sigma(1, 4) = pow(frq_pth(1, t - 1) * frq_pth(3, t - 1), 0.5);
    // sigma(1, 5) = 0;
    // sigma(2, 0) = 0;
    sigma(2, 1) = -pow(frq_pth(2, t - 1) * frq_pth(0, t - 1), 0.5);
    // sigma(2, 2) = 0;
    sigma(2, 3) = -pow(frq_pth(2, t - 1) * frq_pth(1, t - 1), 0.5);
    // sigma(2, 4) = 0;
    sigma(2, 5) = pow(frq_pth(2, t - 1) * frq_pth(3, t - 1), 0.5);
    // sigma(3, 0) = 0;
    // sigma(3, 1) = 0;
    sigma(3, 2) = -pow(frq_pth(3, t - 1) * frq_pth(0, t - 1), 0.5);
    // sigma(3, 3) = 0;
    sigma(3, 4) = -pow(frq_pth(3, t - 1) * frq_pth(1, t - 1), 0.5);
    sigma(3, 5) = -pow(frq_pth(3, t - 1) * frq_pth(2, t - 1), 0.5);

    // proceed the Euler-Maruyama scheme
    frq_pth.col(t) = frq_pth.col(t - 1) + mu * dt + sigma * dW.col(t - 1);

    // remove the noise from the numerical techniques
    for (arma::uword i = 0; i < 4; i++) {
      if (frq_pth(i, t) < 0) {
        frq_pth(i, t) = 0;
      }
      if (frq_pth(i, t) > 1) {
        frq_pth(i, t) = 1;
      }
    }
    frq_pth.col(t) = frq_pth.col(t) / sum(frq_pth.col(t));
  }

  // return the haplotype frequency trajectories under the Wright-Fisher diffusion
  return frq_pth;
}
/*************************/


/***** MomentApprox ******/
// Approximate the first two moments of the Wright-Fisher model using Monte Carlo simulation
// [[Rcpp::export]]
List approximateMoment_MonteCarlo_arma(const arma::dmat& fts_mat, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // simulate the haplotype frequency trajectories under the Wright-Fisher model
  arma::dcube frq_pth(4, arma::uword(lst_gen - int_gen) + 1, sim_num);
  for (arma::uword i = 0; i < sim_num; i++) {
    // cout << "iteration: " << i + 1 << endl;
    frq_pth.slice(i) = simulateWFM_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen);
  }

  // declare the mean vector and variance matrix
  arma::dmat mu = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen));
  arma::dcube sigma = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen));

  // calculate the mean vector and variance matrix
  for (arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    arma::dmat frq_smp = frq_pth.col(k + 1);
    mu.col(k) = arma::mean(frq_smp, 1);
    sigma.slice(k) = arma::cov(frq_smp.t());
  }

  // return the Monte Carlo approximations for the mean vector and variance matrix of the Wright-Fisher model
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}

// Calculate the mean vector for the Wright-Fisher model
// [[Rcpp::export]]
arma::dcolvec calculateMean_arma(const arma::dcolvec& hap_frq, const arma::dmat& fts_mat, const double& rec_rat) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // initialise eta
  arma::dcolvec eta = {-1.0, 1.0, 1.0, -1.0};

  // calculate the mean vector
  arma::dcolvec mu = hap_frq;
  mu = hap_frq % (fts_mat * hap_frq) / arma::as_scalar(hap_frq.t() * fts_mat * hap_frq);
  mu = mu + eta * rec_rat * (mu(0) * mu(3) - mu(1) * mu(2));

  // return the mean vector of the Wright-Fisher model
  return mu;
}

// Calculate the Jacobian matrix over the mean vector for the Wright-Fisher model
// [[Rcpp::export]]
arma::dmat calculateJacobianMean_arma(const arma::dcolvec& hap_frq, const arma::dmat& fts_mat, const double& rec_rat) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the marginal and mean fitness
  arma::dcolvec fts_vec = fts_mat * hap_frq;
  double fts_val = arma::as_scalar(hap_frq.t() * fts_mat * hap_frq);

  // initialise eta
  arma::dcolvec eta = {-1.0, 1.0, 1.0, -1.0};

  //
  arma::dcolvec q = (hap_frq % fts_vec) / fts_val;
  arma::dmat Q = arma::zeros<arma::dmat>(4, 4);
  for (arma::uword i = 0; i < 4; i++) {
    Q.col(i) = q + pow(-1, i) * eta * q(3 - i);
  }

  // calculate the Jacobian matrix over the mean vector
  arma::dmat Jacobian_mu = (arma::eye<arma::dmat>(4, 4) + rec_rat * Q) * (arma::diagmat(hap_frq) * (fts_mat / fts_val - 2 * (fts_vec / fts_val) * (fts_vec.t() / fts_val)) + arma::diagmat(fts_vec / fts_val));

  // return the Jacobian matrix over the mean vector
  return Jacobian_mu;
}

// Calculate the Hessian matrix over the mean vector for the Wright-Fisher model
// [[Rcpp::export]]
arma::cube calculateHessianMean_arma(const arma::dcolvec& hap_frq, const arma::dmat& fts_mat, const double& rec_rat) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the marginal and mean fitness
  arma::dcolvec fts_vec = fts_mat * hap_frq;
  double fts_val = arma::as_scalar(hap_frq.t() * fts_mat * hap_frq);

  //
  arma::dcolvec q = (hap_frq % fts_vec) / fts_val;

  // calculate the Hessian matrix over the mean vector
  arma::dmat Jacobian_q = arma::zeros<arma::dmat>(4, 4);
  for (arma::uword i = 0; i < 4; i++) {
    for (arma::uword j = 0; j < 4; j++) {
      Jacobian_q(i, j) = (fts_mat(i, j) * hap_frq(i) + fts_vec(i) * ((i == j) ? 1.0 : 0.0)) / fts_val - 2 * fts_vec(i) * fts_vec(j) * hap_frq(i) / fts_val / fts_val;
    }
  }
  arma::cube Hessian_q = arma::zeros<arma::dcube>(4, 4, 4);
  for (arma::uword i = 0; i < 4; i++) {
    for (arma::uword j = 0; j < 4; j++) {
      for (arma::uword k = 0; k < 4; k++) {
        Hessian_q(i, j, k) = fts_mat(k, j) * (((k == i) ? 1.0 : 0.0) + ((i == j) ? 1.0 : 0.0)) / fts_val - 2 * (fts_mat(i, j) * fts_vec(k) * hap_frq(i) + fts_mat(i, k) * fts_vec(j) * hap_frq(i) + fts_mat(j, k) * fts_vec(i) * hap_frq(i) + fts_vec(j) * fts_vec(k) * (((j == i) ? 1.0 : 0.0) + ((i == k) ? 1.0 : 0.0))) / fts_val / fts_val + 8 * (fts_vec(i) * fts_vec(j) * fts_vec(k) * hap_frq(i)) / fts_val / fts_val / fts_val;
      }
    }
  }
  arma::cube Hessian_mu = arma::zeros<arma::dcube>(4, 4, 4);
  for (arma::uword i = 0; i < 4; i++) {
    for (arma::uword j = 0; j < 4; j++) {
      Hessian_mu(0, i, j) = (1 - rec_rat) * Hessian_q(0, i, j) + rec_rat * ((Hessian_q(0, i, j) + Hessian_q(1, i, j)) * (q(0) + q(2)) + (Jacobian_q(0, j) + Jacobian_q(1, j)) * (Jacobian_q(0, i) + Jacobian_q(2, i)) + (Jacobian_q(0, i) + Jacobian_q(1, i)) * (Jacobian_q(0, j) + Jacobian_q(2, j)) + (q(0) + q(1)) * (Hessian_q(0, i, j) + Hessian_q(2, i, j)));
      Hessian_mu(1, i, j) = (1 - rec_rat) * Hessian_q(1, i, j) + rec_rat * ((Hessian_q(0, i, j) + Hessian_q(1, i, j)) * (q(1) + q(3)) + (Jacobian_q(0, j) + Jacobian_q(1, j)) * (Jacobian_q(1, i) + Jacobian_q(3, i)) + (Jacobian_q(0, i) + Jacobian_q(1, i)) * (Jacobian_q(1, j) + Jacobian_q(3, j)) + (q(0) + q(1)) * (Hessian_q(1, i, j) + Hessian_q(3, i, j)));
      Hessian_mu(2, i, j) = (1 - rec_rat) * Hessian_q(2, i, j) + rec_rat * ((Hessian_q(2, i, j) + Hessian_q(3, i, j)) * (q(0) + q(2)) + (Jacobian_q(2, j) + Jacobian_q(3, j)) * (Jacobian_q(0, i) + Jacobian_q(2, i)) + (Jacobian_q(2, i) + Jacobian_q(3, i)) * (Jacobian_q(0, j) + Jacobian_q(2, j)) + (q(2) + q(3)) * (Hessian_q(0, i, j) + Hessian_q(2, i, j)));
      Hessian_mu(3, i, j) = (1 - rec_rat) * Hessian_q(3, i, j) + rec_rat * ((Hessian_q(2, i, j) + Hessian_q(3, i, j)) * (q(1) + q(3)) + (Jacobian_q(2, j) + Jacobian_q(3, j)) * (Jacobian_q(1, i) + Jacobian_q(3, i)) + (Jacobian_q(2, i) + Jacobian_q(3, i)) * (Jacobian_q(1, j) + Jacobian_q(3, j)) + (q(2) + q(3)) * (Hessian_q(1, i, j) + Hessian_q(3, i, j)));
    }
  }

  // return the Hessian matrix over the mean vector
  return Hessian_mu;
}

// Approximate the first two moments of the Wright-Fisher model using the extension of Lacerda & Seoighe (2014)
// [[Rcpp::export]]
List approximateMoment_Lacerda_arma(const arma::dmat& fts_mat, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the mean vector and variance matrix
  arma::dmat mu = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);

  // calculate the mean vector and variance matrix
  mu.col(0) = int_frq;
  // sigma.slice(0) = arma::zeros<arma::dmat>(4, 4);
  for (arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    arma::dmat Jacobian_mu = calculateJacobianMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    mu.col(k) = calculateMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    sigma.slice(k) = (1.0 / 2 / pop_siz) * (arma::diagmat(mu.col(k)) - mu.col(k) * mu.col(k).t()) + Jacobian_mu * sigma.slice(k - 1) * Jacobian_mu.t();
  }
  mu.shed_col(0);
  sigma.shed_slice(0);

  // return the approximations for the mean vector and variance matrix of the Wright-Fisher model
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}

// Approximate the first two moments of the Wright-Fisher model using the extension of Terhorst et al. (2015)
// [[Rcpp::export]]
List approximateMoment_Terhorst_arma(const arma::dmat& fts_mat, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the mean vector, variance matrix and error terms
  arma::dmat x = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dmat mu = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);
  arma::dmat epsilon = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube xi = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);

  // calculate the mean vector and variance matrix
  x.col(0) = int_frq;
  mu.col(0) = int_frq;
  // sigma.slice(0) = arma::zeros<arma::dmat>(4, 4);
  // epsilon.col(0) = arma::zeros<arma::dcolvec>(4);
  // xi.slice(0) = arma::zeros<arma::dmat>(4, 4);
  for (arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    arma::dmat Jacobian_mu = calculateJacobianMean_arma(x.col(k - 1), fts_mat, rec_rat);
    arma::dcube Hessian_mu = calculateHessianMean_arma(x.col(k - 1), fts_mat, rec_rat);
    arma::dcolvec term = arma::zeros<arma::dcolvec>(4);
    for (arma::uword i = 1; i < 4; i++) {
      term(i) = 0.5 * arma::as_scalar(sum(sum(arma::dmat(Hessian_mu.row(i)) % arma::dmat(xi.slice(k - 1)), 0), 1));
    }

    x.col(k) = calculateMean_arma(x.col(k - 1), fts_mat, rec_rat);
    mu.col(k) = x.col(k) + Jacobian_mu * epsilon.col(k - 1) + term;
    sigma.slice(k) = (1.0 / 2 / pop_siz) * (arma::diagmat(mu.col(k)) - (x.col(k) + Jacobian_mu * epsilon.col(k - 1)) * (x.col(k) + Jacobian_mu * epsilon.col(k - 1)).t()) + (1 - 1.0 / 2 / pop_siz) * Jacobian_mu * sigma.slice(k - 1) * Jacobian_mu.t();
    epsilon.col(k) = mu.col(k) - x.col(k);
    xi.slice(k) = sigma.slice(k) + epsilon.col(k) * epsilon.col(k).t();
  }
  mu.shed_col(0);
  sigma.shed_slice(0);

  // return the approximations for the mean vector and variance matrix of the Wright-Fisher model
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}

// Approximate the first two moments of the Wright-Fisher model using the extension of Paris et al. (2019) with the first-order Taylor expansion
// [[Rcpp::export]]
List approximateMoment_Paris1_arma(const arma::dmat& fts_mat, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the mean vector and variance matrix
  arma::dmat mu = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);

  // calculate the mean vector and variance matrix
  mu.col(0) = int_frq;
  // sigma.slice(0) = arma::zeros<arma::dmat>(4, 4);
  for (arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    arma::dmat Jacobian_mu = calculateJacobianMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    mu.col(k) = calculateMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    sigma.slice(k) = (1.0 / 2 / pop_siz) * (arma::diagmat(mu.col(k)) - mu.col(k) * mu.col(k).t()) + (1 - 1.0 / 2 / pop_siz) * Jacobian_mu * sigma.slice(k - 1) * Jacobian_mu.t();
  }
  mu.shed_col(0);
  sigma.shed_slice(0);

  // return the approximations for the mean vector and variance matrix of the Wright-Fisher model
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}

// Approximate the first two moments of the Wright-Fisher model using the extension of Paris et al. (2019) with the second-order Taylor expansion
// [[Rcpp::export]]
List approximateMoment_Paris2_arma(const arma::dmat& fts_mat, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the mean vector and variance matrix
  arma::dmat mu = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);

  // calculate the mean vector and variance matrix
  mu.col(0) = int_frq;
  // sigma.slice(0) = arma::zeros<arma::dmat>(4, 4);
  for (arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    arma::dmat Jacobian_mu = calculateJacobianMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    arma::dcube Hessian_mu = calculateHessianMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    arma::dcolvec term = arma::zeros<arma::dcolvec>(4);
    for (arma::uword i = 1; i < 4; i++) {
      term(i) = 0.5 * arma::as_scalar(sum(sum(arma::dmat(Hessian_mu.row(i)) % arma::dmat(sigma.slice(k - 1)), 0), 1));
    }

    mu.col(k) = calculateMean_arma(mu.col(k - 1), fts_mat, rec_rat) + term;
    sigma.slice(k) = (1.0 / 2 / pop_siz) * (arma::diagmat(mu.col(k)) - mu.col(k) * (mu.col(k)).t()) + (1 - 1.0 / 2 / pop_siz) * Jacobian_mu * sigma.slice(k - 1) * Jacobian_mu.t();
  }
  mu.shed_col(0);
  sigma.shed_slice(0);

  // return the approximations for the mean vector and variance matrix of the Wright-Fisher model
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}
/*************************/


/****** ParamApprox ******/
// Approximate the Wright-Fisher model using the normal distribution
// [[Rcpp::export]]
List calculateParam_Norm_arma(const arma::dmat& mean, const arma::dcube& variance, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the parameters of the normal approximation
  arma::dmat m = mean.rows(0, 2);
  arma::dcube V = variance.tube(0, 0, 2, 2);

  // return the parameters of the normal approximation
  return List::create(Named("mean", m),
                      Named("variance", V));
}

// Generate the samples under the normal approximation of the Wright-Fisher model
// [[Rcpp::export]]
arma::dcube generateSample_Norm_arma(const arma::dmat& mean, const arma::dcube& variance, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // simulate the haplotype frequency trajectories under the normal approximation
  arma::dcube frq_pth = arma::zeros<arma::dcube>(4, arma::uword(lst_gen - int_gen), sim_num);
  for (arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    arma::dcolvec eigval;
    arma::dmat eigvec;
    arma::eig_sym(eigval, eigvec, variance.slice(k));
    eigval.elem(arma::find(eigval < 0)).zeros();

    arma::dmat hap_frq = arma::randn<arma::dmat>(3, sim_num);
    hap_frq = pow(eigval, 0.5) % hap_frq.each_col();
    hap_frq = eigvec * hap_frq;
    hap_frq = mean.col(k) + hap_frq.each_col();
    hap_frq.insert_rows(3, 1);
    hap_frq.row(3) = 1 - hap_frq.row(0) - hap_frq.row(1) - hap_frq.row(2);
    frq_pth.col(k) = hap_frq;
  }

  // return the haplotype frequency trajectories under the normal approximation
  return frq_pth;
}

// Approximate the moments of the normal approximation of the Wright-Fisher model
// [[Rcpp::export]]
List approximateMoment_Norm_arma(const arma::dmat& mean, const arma::dcube& variance, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // simulate the haplotype frequency trajectories
  arma::dcube frq_pth = generateSample_Norm_arma(mean, variance, int_gen, lst_gen, sim_num);

  // calculate the mean vector and variance matrix
  arma::dmat m = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen));
  arma::dcube V = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen));
  for (arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    arma::dmat frq_smp = frq_pth.col(k);
    for (arma::uword i = 0; i < 4; i++) {
      arma::drowvec frq = frq_smp.row(i);
      frq_smp.shed_cols(arma::find(frq < 0));
    }
    m.col(k) = arma::mean(frq_smp, 1);
    V.slice(k) = arma::cov(frq_smp.t());
  }

  // return the approximations for the mean vector and variance matrix of the normal approximation
  return List::create(Named("mean", m),
                      Named("variance", V));
}

// Calculate the location vector for the logistic normal approximation
// [[Rcpp::export]]
arma::dmat calculateLocation_arma(const arma::dmat& hap_frq) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the location vector
  arma::dmat phi_frq = arma::zeros<arma::dmat>(3, hap_frq.n_cols);
  phi_frq.row(0) = log(hap_frq.row(0) / hap_frq.row(3));
  phi_frq.row(1) = log(hap_frq.row(1) / hap_frq.row(3));
  phi_frq.row(2) = log(hap_frq.row(2) / hap_frq.row(3));

  // return the location vector
  return phi_frq;
}

// Calculate the Jacobian matrix over the location vector for the logistic normal approximation
// [[Rcpp::export]]
arma::dcube calculateJacobianLocation_arma(const arma::dmat& hap_frq) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the Jacobian matrix over the location vector
  arma::dcube Jacobian_phi = arma::zeros<arma::dcube>(3, 4, hap_frq.n_cols);
  Jacobian_phi.tube(0, 0) = 1.0 / hap_frq.row(0);
  Jacobian_phi.tube(1, 1) = 1.0 / hap_frq.row(1);
  Jacobian_phi.tube(2, 2) = 1.0 / hap_frq.row(2);
  Jacobian_phi.tube(0, 3) = -1.0 / hap_frq.row(3);
  Jacobian_phi.tube(1, 3) = -1.0 / hap_frq.row(3);
  Jacobian_phi.tube(2, 3) = -1.0 / hap_frq.row(3);

  // return the Jacobian matrix over the location vector
  return Jacobian_phi;
}

// Approximate the Wright-Fisher model using the logistic normal distribution
// [[Rcpp::export]]
List calculateParam_LogisticNorm_arma(const arma::dmat& mean, const arma::dcube& variance, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the parameters of the logistic normal approximation
  arma::dmat location = calculateLocation_arma(mean);
  arma::dcube scalesq = arma::zeros<arma::dcube>(3, 3, arma::uword(lst_gen - int_gen));
  arma::dcube Jacobian_location = calculateJacobianLocation_arma(mean);
  for (arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    scalesq.slice(k) = Jacobian_location.slice(k) * variance.slice(k) * Jacobian_location.slice(k).t();
  }

  // return the parameters of the logistic normal approximation
  return List::create(Named("location", location),
                      Named("scalesq", scalesq));
}

// Calculate the parameters for the logistic normal approximation of the Wright-Fisher model
// [[Rcpp::export]]
arma::dmat calculateALT_arma(const arma::dmat& phi_frq) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the additive logistic transformation
  arma::dmat hap_frq = arma::zeros<arma::dmat>(4, phi_frq.n_cols);
  hap_frq.row(0) = exp(phi_frq.row(0));
  hap_frq.row(1) = exp(phi_frq.row(1));
  hap_frq.row(2) = exp(phi_frq.row(2));
  hap_frq.row(3) = arma::ones<arma::drowvec>(phi_frq.n_cols);
  hap_frq = arma::normalise(hap_frq, 1, 0);

  // return the additive logistic transformation
  return hap_frq;
}

// Generate the samples under the logistic normal approximation of the Wright-Fisher model
// [[Rcpp::export]]
arma::dcube generateSample_LogisticNorm_arma(const arma::dmat& location, const arma::dcube& scalesq, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // simulate the haplotype frequency trajectories under the logistic normal approximation
  arma::dcube frq_pth = arma::zeros<arma::dcube>(4, arma::uword(lst_gen - int_gen), sim_num);
  for (arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    cout << "iteration: " << k << endl;
    arma::dcolvec eigval;
    arma::dmat eigvec;
    arma::eig_sym(eigval, eigvec, scalesq.slice(k));
    eigval.elem(arma::find(eigval < 0)).zeros();

    arma::dmat phi_frq = arma::randn<arma::dmat>(3, sim_num);
    phi_frq = pow(eigval, 0.5) % phi_frq.each_col();
    phi_frq = eigvec * phi_frq;
    phi_frq = location.col(k) + phi_frq.each_col();
    frq_pth.col(k) = calculateALT_arma(phi_frq);

    // arma::dmat phi_frq = arma::mvnrnd(location.col(k), scalesq.slice(k), sim_num);
    // frq_pth.col(k) = calculateALT_arma(phi_frq);
  }

  // return the haplotype frequency trajectories under the logistic normal approximation
  return frq_pth;
}

// Approximate the moments of the logistic normal approximation of the Wright-Fisher model
// [[Rcpp::export]]
List approximateMoment_LogisticNorm_arma(const arma::dmat& location, const arma::dcube& scalesq, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // simulate the haplotype frequency trajectories
  arma::dcube frq_pth = generateSample_LogisticNorm_arma(location, scalesq, int_gen, lst_gen, sim_num);

  // calculate the mean vector and variance matrix
  arma::dmat mean = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen));
  arma::dcube variance = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen));
  for (arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    arma::dmat frq_smp = frq_pth.col(k);
    mean.col(k) = arma::mean(frq_smp, 1);
    variance.slice(k) = arma::cov(frq_smp.t());
  }

  // return the approximations for the mean vector and variance matrix of the logistic normal approximation
  return List::create(Named("mean", mean),
                      Named("variance", variance));
}

// Calculate the permutation matrix
// [[Rcpp::export]]
arma::dmat calculatePermutationMatrix_arma() {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat perm = arma::zeros<arma::dmat>(24, 4);
  perm(arma::span(0, 5), 0) = 0 * arma::ones<arma::dcolvec>(6);
  perm(arma::span(0, 1), 1) = 1 * arma::ones<arma::dcolvec>(2);
  perm(arma::span(2, 3), 1) = 2 * arma::ones<arma::dcolvec>(2);
  perm(arma::span(4, 5), 1) = 3 * arma::ones<arma::dcolvec>(2);
  perm(arma::span(0, 5), 2) = {2, 3, 1, 3, 1, 2};
  perm(arma::span(0, 5), 3) = {3, 2, 3, 1, 2, 1};
  perm(arma::span(6, 11), arma::span(0, 2)) = perm(arma::span(0, 5), arma::span(1, 3));
  perm(arma::span(6, 11), 3) = perm(arma::span(0, 5), 0);
  perm(arma::span(12, 17), arma::span(0, 2)) = perm(arma::span(6, 11), arma::span(1, 3));
  perm(arma::span(12, 17), 3) = perm(arma::span(6, 11), 0);
  perm(arma::span(18, 23), arma::span(0, 2)) = perm(arma::span(12, 17), arma::span(1, 3));
  perm(arma::span(18, 23), 3) = perm(arma::span(12, 17), 0);

  return perm;
}

// Calculate the parameters for the hierarchical beta approximation of the Wright-Fisher model (a single generation)
// [[Rcpp::export]]
arma::dmat calculateParam_HierarchicalBeta_OneStep_arma(const arma::dcolvec& mean, const arma::dmat& variance) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the mean of the beta distributions
  arma::dcolvec m = arma::zeros<arma::dcolvec>(3);
  m(0) = mean(0) + mean(1);
  m(1) = mean(0) / (mean(0) + mean(1));
  m(2) = mean(2) / (mean(2) + mean(3));

  // calculate the variance of the beta distributions
  arma::dcolvec V = arma::zeros<arma::dcolvec>(3);
  V(0) = (variance(0, 0) - variance(1, 1)) / (m(1) * m(1) - (1 - m(1)) * (1 - m(1)));
  V(1) = (variance(0, 0) + variance(1, 1) - (m(1) * m(1) + (1 - m(1)) * (1 - m(1))) * V(0)) / 2 / (m(0) * m(0) + V(0));
  // V(2) = (variance(2, 2) + variance(3, 3) - (m(2) * m(2) + (1 - m(2)) * (1 - m(2))) * V(0)) / 2 / ((1 - m(0)) * (1 - m(0)) + V(0));
  V(2) = (variance(2, 2) - m(2) * m(2) * V(0)) / ((1 - m(0)) * (1 - m(0)) + V(0));

  // calculate the parameters of the hierarchical beta approximation
  arma::dmat param = arma::zeros<arma::dmat>(3, 2);
  param.col(0) = (m % (1 - m) / V - 1) % m;
  param.col(1) = (m % (1 - m) / V - 1) % (1 - m);

  // return the parameters of the hierarchical beta approximation
  return param;
}

// Calculate the parameters for the hierarchical beta approximation of the Wright-Fisher model
// [[Rcpp::export]]
List calculateParam_HierarchicalBeta_arma(const arma::dmat& mean, const arma::dcube& variance, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the parameters of the hierarchical beta approximation
  arma::dmat alpha = arma::zeros<arma::dmat>(3, arma::uword(lst_gen - int_gen));
  arma::dmat beta = arma::zeros<arma::dmat>(3, arma::uword(lst_gen - int_gen));
  arma::dmat perm = calculatePermutationMatrix_arma();
  arma::dmat index = arma::zeros<arma::dmat>(2, arma::uword(lst_gen - int_gen));
  arma::dcube trans = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen));

  for (arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    // calculate the parameters of the hierarchical beta approximation
    arma::dmat param = calculateParam_HierarchicalBeta_OneStep_arma(mean.col(k), variance.slice(k));
    if (param.min() < 0) {
      for (arma::uword i = 1; i < 24; i++) {
        arma::dmat perm_mat = arma::zeros<arma::dmat>(4, 4);
        perm_mat(0, perm(i, 0)) = 1;
        perm_mat(1, perm(i, 1)) = 1;
        perm_mat(2, perm(i, 2)) = 1;
        perm_mat(3, perm(i, 3)) = 1;

        arma::dcolvec mean_perm = perm_mat * mean.col(k);
        arma::dmat variance_perm = perm_mat * variance.slice(k) * perm_mat;
        arma::dmat param_perm = calculateParam_HierarchicalBeta_OneStep_arma(mean_perm, variance_perm);
        if (param_perm.min() > 0 && param_perm.is_finite() == true) {
          index(1, k) = 1;
          trans.slice(k) = perm_mat;
          alpha.col(k) = param_perm.col(0);
          beta.col(k) = param_perm.col(1);

          break;
        }
      }

      double delta = (2 * 5000 - 1); // an ad-hoc procedure for avoiding non-positive shape parameters
      if (alpha(0, k) == 0 || beta(0, k) == 0) {
        arma::uvec order_vec = arma::sort_index(mean.col(k));
        arma::dmat perm_mat = arma::zeros<arma::dmat>(4, 4);
        perm_mat(0, order_vec(0)) = 1;
        perm_mat(1, order_vec(1)) = 1;
        perm_mat(2, order_vec(2)) = 1;
        perm_mat(3, order_vec(3)) = 1;

        trans.slice(k) = perm_mat;
        arma::dcolvec mean_perm = trans.slice(k) * mean.col(k);
        arma::dmat variance_perm = (trans.slice(k)) * variance.slice(k) * (trans.slice(k));
        arma::dmat param_perm = calculateParam_HierarchicalBeta_OneStep_arma(mean_perm, variance_perm);
        if (param_perm(0, 0) > 0 && param_perm.row(0).is_finite() == true) {
          alpha(0, k) = param_perm(0, 0);
          beta(0, k) = param_perm(0, 1);
        } else {
          alpha(0, k) = delta * (mean_perm(0) + mean_perm(1));
          beta(0, k) = delta * (1.0 - mean_perm(0) - mean_perm(1));
        }
        if (param_perm(1, 0) > 0 && param_perm.row(1).is_finite() == true) {
          alpha(1, k) = param_perm(1, 0);
          beta(1, k) = param_perm(1, 1);
        } else {
          alpha(1, k) = delta * (mean_perm(0) / (mean_perm(0) + mean_perm(1)));
          beta(1, k) = delta * (1.0 - mean_perm(0) / (mean_perm(0) + mean_perm(1)));
        }
        if (param_perm(2, 0) > 0 && param_perm.row(2).is_finite() == true) {
          alpha(2, k) = param_perm(2, 0);
          beta(2, k) = param_perm(2, 1);
        } else {
          alpha(2, k) = delta * (mean_perm(2) / (mean_perm(2) + mean_perm(3)));
          beta(2, k) = delta * (1.0 - mean_perm(2) / (mean_perm(2) + mean_perm(3)));
        }
      }
    } else{
      alpha.col(k) = param.col(0);
      beta.col(k) = param.col(1);
      index(0, k) = 1;
      trans.slice(k) = arma::eye(4, 4);
    }
  }

  // return the parameters of the hierarchical beta approximation
  return List::create(Named("alpha", alpha),
                      Named("beta", beta),
                      Named("index", index),
                      Named("trans", trans));
}

// Generate the samples under the hierarchical beta approximation of the Wright-Fisher model
// [[Rcpp::export]]
arma::dcube generateSample_HierarchicalBeta_arma(const arma::dmat& alpha, const arma::dmat& beta, const arma::dcube& trans, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // simulate the haplotype frequency trajectories
  arma::dcube frq_pth(4, arma::uword(lst_gen - int_gen), sim_num);
  for (arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    arma::dmat psi_frq = arma::zeros<arma::dmat>(3, sim_num);
    psi_frq.row(0) = as<arma::drowvec>(Rcpp::rbeta(sim_num, alpha(0, k), beta(0, k)));
    psi_frq.row(1) = as<arma::drowvec>(Rcpp::rbeta(sim_num, alpha(1, k), beta(1, k)));
    psi_frq.row(2) = as<arma::drowvec>(Rcpp::rbeta(sim_num, alpha(2, k), beta(2, k)));

    frq_pth.tube(0, k) = psi_frq.row(0) % psi_frq.row(1);
    frq_pth.tube(1, k) = psi_frq.row(0) % (1 - psi_frq.row(1));
    frq_pth.tube(2, k) = (1 - psi_frq.row(0)) % psi_frq.row(2);
    frq_pth.tube(3, k) = (1 - psi_frq.row(0)) % (1 - psi_frq.row(2));

    frq_pth.col(k) = trans.slice(k) * arma::dmat(frq_pth.col(k));
  }

  // return the haplotype frequency trajectories under the hierarchical beta approximation
  return frq_pth;
}

// Approximate the moments of the hierarchical beta approximation of the Wright-Fisher model
// [[Rcpp::export]]
List approximateMoment_HierarchicalBeta_arma(const arma::dmat& alpha, const arma::dmat& beta, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the mean vector and variance matrix
  arma::dmat mean = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen));
  arma::dcube variance = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen));
  for (arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    // calculate the mean and variance of the beta distributions
    arma::dcolvec m = alpha.col(k) / (alpha.col(k) + beta.col(k));
    arma::dcolvec V = alpha.col(k) % beta.col(k) / (alpha.col(k) + beta.col(k)) / (alpha.col(k) + beta.col(k)) / (alpha.col(k) + beta.col(k) + 1);

    // calculate the mean vector
    mean(0, k) = m(0) * m(1);
    mean(1, k) = m(0) * (1 - m(1));
    mean(2, k) = (1 - m(0)) * m(2);
    mean(3, k) = (1 - m(0)) * (1 - m(2));

    // calculate the variance matrix
    variance(0, 0, k) = V(0) * m(1) * m(1) + (V(0) + m(0) * m(0)) * V(1);
    variance(0, 1, k) = -(V(0) + m(0) * m(0)) * V(1) + V(0) * m(1) * (1 - m(1));
    variance(0, 2, k) = -V(0) * m(1) * m(2);
    variance(0, 3, k) = -V(0) * m(1) * (1 - m(2));
    variance(1, 0, k) = variance(0, 1, k);
    variance(1, 1, k) = V(0) * (1 - m(1)) * (1 - m(1)) + (V(0) + m(0) * m(0)) * V(1);
    variance(1, 2, k) = -V(0) * (1 - m(1)) * m(2);
    variance(1, 3, k) = -V(0) * (1 - m(1)) * (1 - m(2));
    variance(2, 0, k) = variance(0, 2, k);
    variance(2, 1, k) = variance(1, 2, k);
    variance(2, 2, k) = V(0) * m(2) * m(2) + (V(0) + (1 - m(0)) * (1 - m(0))) * V(2);
    variance(2, 3, k) = -(V(0) + (1 - m(0)) * (1 - m(0))) * V(2) + V(0) * m(2) * (1 - m(2));
    variance(3, 0, k) = variance(0, 3, k);
    variance(3, 1, k) = variance(1, 3, k);
    variance(3, 2, k) = variance(2, 3, k);
    variance(3, 3, k) = V(0) * (1 - m(2)) * (1 - m(2)) + (V(0) + (1 - m(0)) * (1 - m(0))) * V(2);
  }

  // return the approximations for the mean vector and variance matrix of the hierarchical beta approximation
  return List::create(Named("mean", mean),
                      Named("variance", variance));
}

// Calculate the parameters for the pyramidal hierarchical beta approximation of the Wright-Fisher model
// [[Rcpp::export]]
List calculateParam_PyramidalHierarchicalBeta_arma(const arma::dmat& mean, const arma::dcube& variance, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the parameters of the pyramidal hierarchical beta approximation
  arma::dmat alpha = arma::zeros<arma::dmat>(6, arma::uword(lst_gen - int_gen));
  arma::dmat beta = arma::zeros<arma::dmat>(6, arma::uword(lst_gen - int_gen));
  for(arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    // calculate the mean of the beta distributions
    arma::dcolvec m = arma::zeros<arma::dcolvec>(6);
    m(0) = 0.5;
    m(1) = 0.5;
    m(2) = 0.5;
    m(3) = 4 * mean(0, k);
    m(4) = (m(3) + 4 * mean(1, k) - 1) / 2;
    m(5) = 1 - 4 * mean(3, k);

    // calculate the variance of the beta distributions
    arma::dcolvec V = arma::zeros<arma::dcolvec>(6);
    V(0) = - 4 * variance(0, 3, k) / (m(3) * (1 - m(5)));
    V(2) = 0.25 - (variance(1, 3, k) / m(4) / (1 - m(5)) + 0.25 * (1 - m(5)) * (V(0) * (1 - m(3) + m(4)) + 0.25 * m(4))) / (V(0) + 0.25);
    V(5) = ((4 * variance(3, 3, k) - V(0) * (1 - m(5)) * (1 - m(5))) / (4 * V(0) + 1) - V(2) * (1 - m(5)) * (1 - m(5))) / (V(2) + 0.25);
    V(1) = - (4 * variance(0, 2, k) + V(0) * m(3) * (m(3) - m(4) + m(5))) / ((4 * V(0) + 1) * m(3) * (1 - m(3)));
    // V(1) = (variance(2, 2, k) - variance(1, 1, k) + variance(0, 0, k) - 0.25 * V(0) * m(3) * m(3) - (V(0) + 0.25) * (V(5) * (0.25 + V(2)) + V(2) * (1 - 2 * m(4) + m(5) * m(5))) - V(0) * (- 2 * (1 - m(4)) * m(5) * V(2) + 0.25 * m(5) * m(5) - 0.25 * (1 - m(3)) * (1 - m(3))) + 0.5 * (1 - m(4)) * m(5) * V(2)) / (- 2 * m(4) * (1 - m(3)) * V(0) + 0.5 * (1 - m(3)) * m(4) + 2 * (V(0) + 0.25) * (m(3) - m(4)));
    V(3) = ((4 * variance(0, 0, k) - V(0) * m(3) * m(3)) / (4 * V(0) + 1) - V(1) * m(3) * m(3)) / (V(1) + 0.25);
    V(4) = (variance(1, 1, k) - (V(0) + 0.25) * (V(3) * (0.25 + V(1)) + V(1) * ((1 - m(3)) * (1 - m(3)) + m(4) * m(4)) + V(2) * m(4) * m(4)) - V(0) * (2 * m(4) * (1 - m(3)) * V(1) + 0.25 * (1 - m(3)) * (1 - m(3))) + 0.5 * (1 - m(3)) * m(4) * V(1)) / ((V(0) + 0.25) * (0.5 + V(1) + V(2)) - V(0) / 2 + 0.125);

    // calculate the parameters of the pyramidal hierarchical beta approximation
    alpha.col(k) = (m % (1 - m) / V - 1) % m;
    beta.col(k) = (m % (1 - m) / V - 1) % (1 - m);
  }

  // return the parameters of the pyramidal hierarchical beta approximation
  return List::create(Named("alpha", alpha),
                      Named("beta", beta));
}

// Generate the samples under the pyramidal hierarchical beta approximation of the Wright-Fisher model
// [[Rcpp::export]]
arma::dcube generateSample_PyramidalHierarchicalBeta_arma(const arma::dmat& alpha, const arma::dmat& beta, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // simulate the haplotype frequency trajectories
  arma::dcube frq_pth(4, arma::uword(lst_gen - int_gen), sim_num);
  for(arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    arma::dmat psi_frq = arma::zeros<arma::dmat>(6, sim_num);
    psi_frq.row(0) = as<arma::drowvec>(Rcpp::rbeta(sim_num, alpha(0, k), beta(0, k)));
    psi_frq.row(1) = as<arma::drowvec>(Rcpp::rbeta(sim_num, alpha(1, k), beta(1, k)));
    psi_frq.row(2) = as<arma::drowvec>(Rcpp::rbeta(sim_num, alpha(2, k), beta(2, k)));
    psi_frq.row(3) = as<arma::drowvec>(Rcpp::rbeta(sim_num, alpha(3, k), beta(3, k)));
    psi_frq.row(4) = as<arma::drowvec>(Rcpp::rbeta(sim_num, alpha(4, k), beta(4, k)));
    psi_frq.row(5) = as<arma::drowvec>(Rcpp::rbeta(sim_num, alpha(5, k), beta(5, k)));

    frq_pth.tube(0, k) = psi_frq.row(0) % psi_frq.row(1) % psi_frq.row(3);
    frq_pth.tube(1, k) = psi_frq.row(0) % (psi_frq.row(1) % (1 - psi_frq.row(3)) + (1 - psi_frq.row(1)) % psi_frq.row(4)) + (1 - psi_frq.row(0)) % psi_frq.row(2) % psi_frq.row(4);
    frq_pth.tube(2, k) = psi_frq.row(0) % (1 - psi_frq.row(1)) % (1 - psi_frq.row(4)) + (1 - psi_frq.row(0)) % (psi_frq.row(2) % (1 - psi_frq.row(4)) + (1 - psi_frq.row(2)) % psi_frq.row(5));
    frq_pth.tube(3, k) = (1 - psi_frq.row(0)) % (1 - psi_frq.row(2)) % (1 - psi_frq.row(5));
  }

  // return the haplotype frequency trajectories under the pyramidal hierarchical beta approximation
  return frq_pth;
}

// Approximate the moments of the pyramidal hierarchical beta approximation of the Wright-Fisher model
// [[Rcpp::export]]
List approximateMoment_PyramidalHierarchicalBeta_arma(const arma::dmat& alpha, const arma::dmat& beta, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the mean vector and variance matrix
  arma::dmat mean = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen));
  arma::dcube variance = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen));
  for(arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    // calculate the mean and variance of the beta distributions
    arma::dcolvec m = alpha.col(k) / (alpha.col(k) + beta.col(k));
    arma::dcolvec V = alpha.col(k) % beta.col(k) / (alpha.col(k) + beta.col(k)) / (alpha.col(k) + beta.col(k)) / (alpha.col(k) + beta.col(k) + 1);

    // calculate the mean vector
    mean(0, k) = m(0) * m(1) * m(3);
    // mean(1, k) = m(0) * (m(1) * (1 - m(3)) + (1 - m(1)) * m(4)) + (1 - m(0)) * m(2) * m(4);
    mean(1, k) = (1 - m(3) + 2 * m(4)) / 4;
    // mean(2, k) = m(0) * (1 - m(1)) * (1 - m(4)) + (1 - m(0)) * (m(2) * (1 - m(4)) + (1 - m(2)) * m(5));
    mean(2, k) = (2 - 2 * m(4) + m(5)) / 4;
    mean(3, k) = (1 - m(0)) * (1 - m(2)) * (1 - m(5));

    // calculate the variance matrix
    variance(0, 0, k) = (V(0) + 0.25) * (V(3) * (0.25 + V(1)) + V(1) * m(3) * m(3))  + V(0) * m(3) * m(3) / 4;
    variance(0, 1, k) = (V(0) + 0.25) * (V(1) * (m(3) * (1 - m(3)) - m(3) * m(4)) - V(3) * (0.25 + V(1))) + V(0) * m(3) / 4 * (1 - m(3));
    variance(0, 2, k) = - (V(0) + 0.25) * V(1) * m(3) * (1 - m(3)) + V(0) * m(3) / 4 * (- m(3) + m(4) - m(5));
    variance(0, 3, k) = - V(0) * m(3) * (1 - m(5)) / 4;
    variance(1, 0, k) = variance(0, 1, k);
    variance(1, 1, k) = (V(0) + 0.25) * (V(3) * (0.25 + V(1)) + V(4) * (0.5 + V(1) + V(2)) + V(2) * m(4) * m(4) + V(1) * ((1 - m(3)) * (1 - m(3)) + m(4) * m(4))) + V(0) * (2 * m(4) * (1 - m(3)) * V(1) + 0.25 * (1 - m(3)) * (1 - m(3)) - 0.5 * V(4)) + 0.125 * V(4) - 0.5 * (1 - m(3)) * m(4) * V(1);
    variance(1, 2, k) = (V(0) + 0.25) * (V(1) * (m(4) * (1 - m(4)) - (1 - m(3)) * (1 - m(4))) + V(2) * (m(4) * (1 - m(4)) - m(4) * m(5)) - V(4) * (0.75 + V(1) + V(2))) + V(0) * (0.5 * m(4) * m(5) - 0.25 * (1 - m(3)) * m(5)) - 0.5 * (0.25 - V(0)) * V(4);
    variance(1, 3, k) = m(4) * (1 - m(5)) * ((V(0) + 0.25) * (0.25 - V(2)) - 0.5 * (1 - m(5)) * (V(0) / 2 * (1 - m(3) + m(4)) + 0.125 * m(4)));
    variance(2, 0, k) = variance(0, 2, k);
    variance(2, 1, k) = variance(1, 2, k);
    variance(2, 2, k) = (V(0) + 0.25) * (V(1) * (1 - m(4)) * (1 - m(4)) + V(4) * (0.5 + V(1) + V(2)) + V(5) * (0.25 + V(2)) + V(2) * ((1 - m(4)) * (1 - m(4)) + m(5) * m(5))) + V(0) * (- (1 - m(4)) * 2 * m(5) * V(2) + 0.25 * m(5) * m(5) - 0.5 * V(4)) + 0.125 * V(4) - 0.5 * (1 - m(4)) * m(5) * V(2);
    variance(2, 3, k) = - (V(0) + 0.25) * (V(2) + 0.25) * V(5) * (1 - m(4)) * (1 - m(5)) * (- V(2) * (V(0) + 0.25)) + m(5) * (1 - m(5)) * (V(0) * (V(2) + 0.25) + 0.25 * V(2));
    variance(3, 0, k) = variance(0, 3, k);
    variance(3, 1, k) = variance(1, 3, k);
    variance(3, 2, k) = variance(2, 3, k);
    variance(3, 3, k) = (V(0) + 0.25) * (V(5) * (0.25 + V(2)) + V(2) * (1 - m(5)) * (1 - m(5))) + V(0) / 4 * (1 - m(5)) * (1 - m(5));
  }

  // return the approximations for the mean vector and variance matrix of the pyramidal hierarchical beta approximation
  return List::create(Named("mean", mean),
                      Named("variance", variance));
}
/*************************/


/********* ECDF *********/
// Generate the grids for empirical cumulative distribution function
// [[Rcpp::export]]
arma::dmat generateFixGrid_2L_arma(const arma::uword& grd_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat frq_grd = arma::zeros<arma::dmat>(1, 4);
  for (arma::uword i = 1; i < grd_num; i++) {
    for (arma::uword j = 1; j < grd_num - i; j++) {
      for (arma::uword k = 1; k < grd_num - i - j; k++) {
        arma::drowvec hap_frq = {double(i), double(j), double(k), double(grd_num - i - j - k)};
        frq_grd.insert_rows(0, 1);
        frq_grd.row(0) = hap_frq / grd_num;
      }
    }
  }

  return frq_grd;
}

// Generate the grids for empirical cumulative distribution function
// [[Rcpp::export]]
arma::dmat generateRndGrid_2L_arma(const arma::uword& grd_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  NumericMatrix frq_grd(grd_num, 4);
  for(int j = 0; j < 4; j++) {
    frq_grd(_, j) = rgamma(grd_num, 1.0, 1.0);
  }
  for(int i = 0; i < grd_num; i++) {
    frq_grd(i, _) = frq_grd(i, _) / sum(frq_grd(i, _));
  }

  return as<arma::dmat>(frq_grd);
}

// calculate the empirical cumulative distribution function
// [[Rcpp::export]]
arma::dmat calculateECDF_2L_arma(const arma::dcube& frq_pth, const arma::uword& sim_num, const arma::dmat& frq_grd) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat cdf = arma::zeros<arma::dmat>(frq_grd.n_rows, frq_pth.n_cols);
  for(arma::uword k = 0; k < frq_pth.n_cols; k++) {
    arma::dmat frq_smp = frq_pth.col(k);
    frq_smp = frq_smp.t();

    arma::dcolvec frq;
    arma::ucolvec idx;
    for(arma::uword i = 0; i < frq_grd.n_rows; i++) {
      frq = frq_smp.col(0);
      idx = frq < frq_grd(i, 0);
      frq = frq_smp.col(1);
      idx = arma::intersect(idx, frq < frq_grd(i, 1));
      frq = frq_smp.col(2);
      idx = arma::intersect(idx, frq < frq_grd(i, 2));

      cdf(i, k) = idx.n_elem;
    }
  }
  cdf = cdf / sim_num;

  return cdf;
}

// Calculate the root mean square deviation between the two empirical cumulative distribution functions
// [[Rcpp::export]]
arma::colvec calculateRMSD_arma(const arma::dcube& smp_mod, const arma::dcube& smp_apx, const arma::uword& sim_num, const arma::uword& grd_num, const arma::dmat& frq_grd) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat cdf_mod = calculateECDF_2L_arma(smp_mod, sim_num, frq_grd);
  arma::dmat cdf_apx = calculateECDF_2L_arma(smp_apx, sim_num, frq_grd);
  arma::colvec dist = arma::zeros<arma::dcolvec>(smp_mod.n_cols);
  for(arma::uword k = 0; k < smp_mod.n_cols; k++) {
    cout << "iteration: " << k << endl;
    dist(k) = pow(sum((cdf_mod.col(k) - cdf_apx.col(k)) % (cdf_mod.col(k) - cdf_apx.col(k))) / grd_num, 0.5);
  }
  return dist;
}

// Generate initial haplotype frequencies (uniform generation from the flat Dirichlet distribution)
// [[Rcpp::export]]
arma::dmat generateInitFreq_arma(const arma::uword& smp_siz) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // arma::dmat hap_frq = arma::zeros<arma::dmat>(4, pcl_num);
  // arma::dmat mut_frq = arma::randu<arma::dmat>(2, pcl_num);
  // arma::drowvec LD = arma::randu<arma::drowvec>(pcl_num);
  // for (arma::uword i = 0; i < pcl_num; i++) {
  //   double a = -mut_frq(0, i) * mut_frq(1, i);
  //   a = (a >= -(1 - mut_frq(0, i)) * (1 - mut_frq(1, i)))? a : -(1 - mut_frq(0, i)) * (1 - mut_frq(1, i));
  //   double b = mut_frq(0, i) * (1 - mut_frq(1, i));
  //   b = (b <= (1 - mut_frq(0, i)) * mut_frq(1, i))? b : (1 - mut_frq(0, i)) * mut_frq(1, i);
  //   LD(i) = a + (b - a) * LD(i);
  // }
  // hap_frq.row(0) = (1 - mut_frq.row(0)) % (1 - mut_frq.row(1)) + LD;
  // hap_frq.row(1) = (1 - mut_frq.row(0)) % mut_frq.row(1) - LD;
  // hap_frq.row(2) = mut_frq.row(0) % (1 - mut_frq.row(1)) - LD;
  // hap_frq.row(3) = mut_frq.row(0) % mut_frq.row(1) + LD;
  //
  // return hap_frq;

  NumericMatrix hap_frq(4, smp_siz);
  for (int i = 0; i < 4; i++) {
    hap_frq(i, _) = rgamma(smp_siz, 1.0, 1.0);
  }
  for (int j = 0; j < smp_siz; j++) {
    hap_frq(_, j) = hap_frq(_, j) / sum(hap_frq(_, j));
  }

  return as<arma::dmat>(hap_frq);
}
/*************************/
