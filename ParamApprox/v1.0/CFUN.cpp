// Moment-based approximations for the Wright-Fisher model of population dynamics under natural selection at two linked loci
// Zhangyi He, Wenyang Lyu, Mark Beaumont and Feng Yu

// version 1.0

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

  // declare and calculate the fitness matrix
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

// Calculate the sampling probabilities for the Wright-Fisher model
// [[Rcpp::export]]
arma::dcolvec calculateSamplingProb_arma(const arma::dcolvec& hap_frq, const arma::dmat& fts_mat, const double& rec_rat) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare eta
  arma::dcolvec eta = arma::ones<arma::dcolvec>(4);
  eta(0) = -1;
  // eta(1) = 1;
  // eta(2) = 1;
  eta(3) = -1;

  // declare and calculate the sampling probabilities
  arma::dcolvec prob = hap_frq;
  prob = hap_frq % (fts_mat * hap_frq) / arma::as_scalar(hap_frq.t() * fts_mat * hap_frq);
  prob = prob + eta * rec_rat * (prob(0) * prob(3) - prob(1) * prob(2));

  // return the sampling probabilities for the Wright-Fisher model
  return prob;
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

  // declare and initialise the haplotype frequencies during a single generation of the life cycle
  arma::dcolvec hap_frq = int_frq;

  // simulate the haplotype frequency trajectories
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    // calculate the haplotype frequencies after genetic recombination
    arma::dcolvec prob = calculateSamplingProb_arma(hap_frq, fts_mat, rec_rat)

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

  // declare delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  // declare delta W
  arma::dmat dW = pow(dt, 0.5) * arma::randn<arma::dmat>(6, arma::uword(lst_gen - int_gen) * ptn_num);

  // declare the haplotype frequency trajectories
  arma::dmat frq_pth(4, arma::uword(lst_gen - int_gen) * ptn_num + 1);

  // initialise the haplotype frequencies in generation 0
  frq_pth.col(0) = int_frq;

  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) * ptn_num + 1; t++) {
    // calculate the drift coefficient vector
    arma::dcolvec mu = arma::zeros<arma::dcolvec>(4);
    mu(0) =  scl_sel_cof_A * frq_pth(0, t - 1) * (frq_pth(2, t - 1) + frq_pth(3, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(1, t - 1)) * dom_par_A + (frq_pth(2, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_A)) + scl_sel_cof_B * frq_pth(0, t - 1) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par_B + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_B)) - 0.5 * scl_rec_rat * (frq_pth(0, t - 1) * frq_pth(3, t - 1) - frq_pth(1, t - 1) * frq_pth(2, t - 1));
    mu(1) =  scl_sel_cof_A * frq_pth(1, t - 1) * (frq_pth(2, t - 1) + frq_pth(3, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(1, t - 1)) * dom_par_A + (frq_pth(2, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_A)) - scl_sel_cof_B * frq_pth(1, t - 1) * (frq_pth(0, t - 1) + frq_pth(2, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par_B + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_B)) + 0.5 * scl_rec_rat * (frq_pth(0, t - 1) * frq_pth(3, t - 1) - frq_pth(1, t - 1) * frq_pth(2, t - 1));
    mu(2) = -scl_sel_cof_A * frq_pth(2, t - 1) * (frq_pth(0, t - 1) + frq_pth(1, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(1, t - 1)) * dom_par_A + (frq_pth(2, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_A)) + scl_sel_cof_B * frq_pth(2, t - 1) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par_B + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_B)) + 0.5 * scl_rec_rat * (frq_pth(0, t - 1) * frq_pth(3, t - 1) - frq_pth(1, t - 1) * frq_pth(2, t - 1));
    mu(3) = -scl_sel_cof_A * frq_pth(3, t - 1) * (frq_pth(0, t - 1) + frq_pth(1, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(1, t - 1)) * dom_par_A + (frq_pth(2, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_A)) - scl_sel_cof_B * frq_pth(3, t - 1) * (frq_pth(0, t - 1) + frq_pth(2, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par_B + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par_B)) - 0.5 * scl_rec_rat * (frq_pth(0, t - 1) * frq_pth(3, t - 1) - frq_pth(1, t - 1) * frq_pth(2, t - 1));

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
    sigma(2, 4) = 0;
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
    for(arma::uword i = 0; i < 4; i++) {
      if(frq_pth(i, t) < 0) {
        frq_pth(i, t) = 0;
      }
      if(frq_pth(i, t) > 1) {
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
List approximateMoment_MonteCarlo_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the fitness matrix
  arma::dmat fts_mat = calculateFitnessMat_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B);

  // simulate the haplotype frequency trajectories
  arma::dcube frq_pth(4, arma::uword(lst_gen - int_gen) + 1, sim_num);
  for(arma::uword i = 1; i < sim_num; i++) {
    frq_pth.slice(i) = simulateWFM_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen);
  }

  // declare the mean vector and variance matrix
  arma::dmat mu = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);

  // calculate the mean vector and variance matrix
  mu.col(0) = int_frq;
  // sigma.slice(0) = arma::zeros<arma::dmat>(4, 4)
  arma::dmat frq_smp = frq_pth.col(k);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    mu.col(k) = arma::mean(frq_smp, 1);
    sigma.slice(k) = arma::cov(frq_smp.t(), frq_smp.t());
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

  // declare eta
  arma::dcolvec eta = arma::ones<arma::dcolvec>(4);
  eta(0) = -1;
  // eta(1) = 1;
  // eta(2) = 1;
  eta(3) = -1;

  // declare and calculate the mean vector
  arma::dcolvec mu = hap_frq;
  mu = hap_frq % (fts_mat * hap_frq) / arma::as_scalar(hap_frq.t() * fts_mat * hap_frq);
  mu = mu + eta * rec_rat * (mu(0) * mu(3) - mu(1) * mu(2));

  // return the mean vector of the Wright-Fisher model
  return mu;
}

// Calculate Jacobian matrix over the mean vector for the Wright-Fisher model
// [[Rcpp::export]]
arma::dmat calculateJacobianMean_arma(const arma::dcolvec& hap_frq, const arma::dmat& fts_mat, const double& rec_rat){
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat grad_mu(4, 4);
  arma::dmat gen_frq = hap_frq * hap_frq.t();

  // viability selection
  arma::dcolvec q = sum(fts_mat % gen_frq, 1) / arma::as_scalar(sum(sum(fts_mat % gen_frq, 0), 1));

  arma::dmat Q(4, 4);
  // declare eta
  arma::dcolvec eta(4);
  eta(0) = -1;
  eta(1) = 1;
  eta(2) = 1;
  eta(3) = -1;
  for(arma::uword i = 1; i < 4; i++) {
    Q.col(i) = q + eta * q(3-i);
    eta = eta * (-1);
  }

  arma::dmat iden_mat(4,4);
  iden_mat.eye();
  grad_mu = (iden_mat + Q) * (arma::diagmat(hap_frq) * (fts_mat / arma::as_scalar(sum(sum(fts_mat % gen_frq, 0), 1)) - 2 * q / hap_frq * (q / hap_frq).t() ) + arma::diagmat(q / hap_frq) );

  return grad_mu;
}

// Approximate the first two moments of the Wright-Fisher model using the extension of Lacerda & Seoighe (2014)
// [[Rcpp::export]]
List approximateMoment_Lacerda_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::dmat fts_mat) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the mean vector and variance matrix
  arma::dmat mu = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);

  // calculate the mean vector and variance matrix
  mu.col(0) = int_frq;
  // sigma.slice(0) = arma::zeros<arma::dmat>(4, 4);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    mu.col(k) = calculateMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    arma::dmat jacobian_mu = calculateJacobianMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    sigma.slice(k) = (arma::diagmat(mu.col(k)) - mu.col(k) * mu.col(k).t()) / 2 / pop_siz + jacobian_mu * sigma.slice(k - 1) * jacobian_mu.t();
  }

  // return the approximations for the mean vector and variance matrix of the Wright-Fisher model
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}

// Approximate the first two moments of the Wright-Fisher model using the extension of Terhorst et al. (2015)
// [[Rcpp::export]]
List approximateMoment_Terhorst_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::dmat fts_mat) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat x(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dmat mu(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dmat epsilon(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma(4, 4, arma::uword(lst_gen - int_gen) + 1);
  mu.col(0) = int_frq;
  x.col(0) = int_frq;
  arma::dmat int_var = arma::zeros<arma::dmat>(4, 4);
  sigma.slice(0) = int_var;

  for(arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    x.col(t) = calculate_p(rec_rat, fts_mat, x.col(t-1));
    arma::dmat grad_mu = calculate_grad_mu(fts_mat, x.col(t-1));
    epsilon.col(t) = grad_mu * epsilon.col(t-1);
    mu.col(t) = x.col(t) + epsilon.col(t);
    sigma.slice(t) = 1 / 2 / pop_siz * arma::diagmat(mu.col(t)) - 1 / 2 / pop_siz * mu.col(t) * (mu.col(t)).t() + (1 - 1 / 2 / pop_siz) * grad_mu * sigma.slice(t-1) * grad_mu.t();
  }

  // return the approximations for the mean and variance of the Wright-Fisher model at each generation from int_gen to lst_gen
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}

// Approximate the first two moments of the Wright-Fisher model using the extension of Paris et al. (2019)
// [[Rcpp::export]]
List approximateMoment_Paris_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::dmat fts_mat) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the mean vector and variance matrix
  arma::dmat mu = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);

  // calculate the mean vector and variance matrix
  mu.col(0) = int_frq;
  // sigma.slice(0) = arma::zeros<arma::dmat>(4, 4);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    mu.col(k) = calculateMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    arma::dmat jacobian_mu = calculateJacobianMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    sigma.slice(k) = (arma::diagmat(mu.col(k)) - mu.col(k) * mu.col(k).t()) / 2 / pop_siz + (1 - 1.0 / 2 / pop_siz) * jacobian_mu * sigma.slice(k - 1) * jacobian_mu.t();
  }

  // return the approximations for the mean vector and variance matrix of the Wright-Fisher model
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}
/*************************/


/****** ParamApprox ******/
// Approximate the Wright-Fisher model using the normal distribution
// [[Rcpp::export]]
List approximatWFM_norm_arma(const arma::dmat& mean, const arma::dcube& variance) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the parameters of the normal approximation
  arma::dmat mu = mean;
  arma::dcube sigma = variance;

  // return the parameters of the normal approximation
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}

// Simulate the Wright-Fisher model using the normal approximation
// [[Rcpp::export]]
arma::dcube simulateWFM_norm_arma(const arma::dmat& mean, const arma::dcube& variance, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // simulate the haplotype frequency trajectories
  arma::dcube frq_pth(4, arma::uword(lst_gen - int_gen) + 1, sim_num);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    frq_pth.col(k) = arma::mvnrnd(mean.col(k), variance.slice(k), sim_num);
  }

  // return the haplotype frequency trajectories under the normal approximation
  return frq_pth;
}

// Calculate the location vector for the logistic normal approximation
// [[Rcpp::export]]
arma::dmat calculateLocation_arma(const arma::dmat& hap_frq){
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare and calculate the location vector
  arma::dmat phi_frq = arma::zeros<arma::dmat>(3, hap_frq.n_cols);
  phi_frq.row(0) = log(hap_frq.row(0) / hap_frq.row(3));
  phi_frq.row(1) = log(hap_frq.row(1) / hap_frq.row(3));
  phi_frq.row(2) = log(hap_frq.row(2) / hap_frq.row(3));

  // return the location vector
  return phi_frq;
}

// Calculate the Jacobian matrix over the location vector for the logistic normal approximation
// [[Rcpp::export]]
arma::dcube calculateJacobianLocation_arma(const arma::dmat& hap_frq){
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare and calculate the Jacobian matrix over the location vector
  arma::dcube jacobian_mat = arma::zeros<arma::dcube>(3, 4, hap_frq.n_cols);
  jacobian_mat.tube(0, 0) = 1.0 / hap_frq.row(0);
  jacobian_mat.tube(1, 1) = 1.0 / hap_frq.row(1);
  jacobian_mat.tube(2, 2) = 1.0 / hap_frq.row(2);
  jacobian_mat.tube(0, 3) = -1.0 / hap_frq.row(3);
  jacobian_mat.tube(1, 3) = -1.0 / hap_frq.row(3);
  jacobian_mat.tube(2, 3) = -1.0 / hap_frq.row(3);

  // return the Jacobian matrix over the location vector
  return jacobian_mat;
}

// Approximate the Wright-Fisher model using the logistic normal distribution
// [[Rcpp::export]]
List approximatWFM_LogisticNorm_arma(const arma::dmat& mean, const arma::dcube& variance, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare and calculate the parameters of the logistic normal approximation
  arma::dmat location = calculateLocation_arma(mean);
  arma::dcube squared_scale = arma::zeros<arma::dcube>(3, 3, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube jacobian_location = calculateJacobianLocation_arma(mean);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    squared_scale.slice(k) = jacobian_location.slice(k) * variance.slice(k) * jacobian_location.slice(k).t();
  }

  // return the parameters of the logistic normal approximation
  return List::create(Named("location", location),
                      Named("squared_scale", squared_scale));
}

// Calculate the additive logistic transformation
// [[Rcpp::export]]
arma::dmat calculateALT_arma(const arma::dmat& phi_frq){
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare and calculate the additive logistic transformation
  arma::dmat hap_frq = arma::zeros<arma::dmat>(4, phi_frq.n_cols);
  hap_frq.row(0) = exp(phi_frq.row(0));
  hap_frq.row(1) = exp(phi_frq.row(1));
  hap_frq.row(2) = exp(phi_frq.row(2));
  hap_frq.row(3) = arma::ones<arma::drowvec>(phi_frq.n_cols);
  hap_frq = arma::normalise(hap_frq, p = 1, dim = 0);

  // return the additive logistic transformation
  return hap_frq;
}

// Simulate the Wright-Fisher model using the logistic normal approximation
arma::dcube simulateWFM_LogisticNorm_arma(const arma::dmat& location, const arma::dcube& squared_scale, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // simulate the haplotype frequency trajectories under the logistic normal approximation
  arma::dcube frq_pth = arma::zeros<arma::dcube>(4, arma::uword(lst_gen - int_gen) + 1, sim_num);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    arma::dmat phi_frq = arma::mvnrnd(location.col(k), squared_scale.slice(k), sim_num);
    frq_pth.col(k) = calculateALT_arma(phi_frq);
  }

  // return the haplotype frequency trajectories under the logistic normal approximation
  return frq_pth;
}

// Approximate the first two moments of the logistic normal approximation of the Wright-Fisher model
List approximateMoment_LogisticNorm_arma(const arma::dmat& location, const arma::dcube& squared_scale, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // simulate the haplotype frequency trajectories
  arma::dcube frq_pth = simulateWFM_LogisticNorm_arma(location, squared_scale, int_gen, lst_gen, sim_num);

  // calculate the mean vector and variance matrix of the logistic normal approximation
  arma::dmat mean = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube variance = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    mean.col(k) = arma::mean(frq_pth.col(k), 1);
    variance.slice(k) = arma::cov(frq_pth.col(k).t(), frq_pth.col(k).t());
  }

  // return the approximations for the mean vector and variance matrix of the logistic normal approximation
  return List::create(Named("mean", mean),
                      Named("variance", variance));
}

// Approximate the Wright-Fisher model using the hierarchical beta distribution
// [[Rcpp::export]]
List approximatWFM_HierarchicalBeta_arma(const arma::dmat& mean, const arma::dcube& variance, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare and calculate the parameters of the hierarchical beta approximation
  arma::dmat alpha = arma::zeros<arma::dmat>(3, arma::uword(lst_gen - int_gen) + 1);
  arma::dmat beta = arma::zeros<arma::dmat>(3, arma::uword(lst_gen - int_gen) + 1);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    // calculate the mean of the beta distributions 
    arma::dcolvec m = arma::zeros<arma::dcolvec>(3);
    m(0) = mean(0, k) + mean(1, k);
    m(1) = mean(0, k) / (mean(0, k) + mean(1, k));
    m(2) = mean(2, k) / (mean(2, k) + mean(3, k));

    // calculate the variance of the beta distributions 
    arma::dcolvec V = arma::zeros<arma::dcolvec>(3);
    V(0) = (variance(0, 0, k) - variance(1, 1, k)) / (m(1) * m(1) - (1 - m(1)) * (1 - m(1)));
    V(1) = (variance(0, 0, k) + variance(1, 1, k) - (m(1) * m(1) + (1 - m(1)) * (1 - m(1))) * V(0)) / 2 / (m(0) * m(0) + V(0));
    V(2) = (variance(2, 2, k) + variance(3, 3, k) - (m(2) * m(2) + (1 - m(2)) * (1 - m(2))) * V(0)) / 2 / ((1 - m(0)) * (1 - m(0)) + V(0));

    // calculate the parameters of the hierarchical beta approximation
    alpha.col(k) = (m % (1 - m) / V - 1) % m;
    beta.col(k) = (m % (1 - m) / V - 1) % (1 - m);
  }

  // return the parameters of the hierarchical beta approximation
  return List::create(Named("alpha", alpha),
                      Named("beta", beta));
}

// Simulate the Wright-Fisher model using the hierarchical beta approximation
arma::dcube simulateWFM_HierarchicalBeta_arma(const arma::dmat& alpha, const arma::dmat& beta, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // simulate the haplotype frequency trajectories under the hierarchical beta approximation
  arma::dcube frq_pth(4, arma::uword(lst_gen - int_gen) + 1, sim_num);
  for(arma::uword k = 0; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    arma::dmat psi_frq = arma::(3, sim_num);
    psi_frq.row(0) = Rcpp::rbeta(sim_num, alpha(0, k), beta(0, k));
    psi_frq.row(1) = Rcpp::rbeta(sim_num, alpha(1, k), beta(1, k));
    psi_frq.row(2) = Rcpp::rbeta(sim_num, alpha(2, k), beta(2, k));

    frq_pth.tube(0, k) = psi_frq.row(0) % psi_frq.row(1);
    frq_pth.tube(1, k) = psi_frq.row(0) % (1 - psi_frq.row(1));
    frq_pth.tube(2, k) = (1 - psi_frq.row(0)) % psi_frq.row(2);
    frq_pth.tube(3, k) = (1 - psi_frq.row(0)) % (1 - psi_frq.row(2));
  }

  // return the haplotype frequency trajectories under the hierarchical beta approximation
  return frq_pth;
}

// Approximate the first two moments of the hierarchical beta approximation of the Wright-Fisher model
List approximateMoment_HierarchicalBeta_arma(const arma::dmat& alpha, const arma::dmat& beta, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the mean vector and variance matrix of the hierarchical beta approximation
  arma::dmat mean = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube variance = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    //
    arma::dcolvec m = alpha.col(k) / (alpha.col(k) + beta.col(k));
    arma::dcolvec V = alpha.col(k) % beta.col(k) / (alpha.col(k) + beta.col(k)) / (alpha.col(k) + beta.col(k)) / (alpha.col(k) + beta.col(k) + 1);

    //
    mean(0, k) = m(0) * m(1);
    mean(1, k) = m(0) * (1.0 - m(1));
    mean(2, k) = (1.0 - m(0)) * m(2);
    mean(3, k) = (1.0 - m(0)) * (1.0 - m(2));

    //
    variance(0, 0, k) = V(0) * m(1) * m(1) + (V(0) + m(0) * m(0)) * V(1);
    variance(0, 1, k) = -(V(0) + m(0) * m(0)) * V(1) + V(0) * m(1) * (1 - m(1));
    variance(0, 2, k) = -V(0) * m(1) * mt(2);
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

/*************************/
