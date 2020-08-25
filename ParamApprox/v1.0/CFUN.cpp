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

  // declare the fitness
  arma::dmat fts_mat = arma::ones<arma::dmat>(4, 4);
  // fts_mat(0, 0) = 1;
  fts_mat(1, 0) = (1 - dom_par_B * sel_cof_B);
  fts_mat(2, 0) = (1 - dom_par_A * sel_cof_A);
  fts_mat(3, 0) = (1 - dom_par_A * sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts_mat(0, 1) = (1 - dom_par_B * sel_cof_B);
  fts_mat(1, 1) = (1 - sel_cof_B);
  fts_mat(2, 1) = (1 - dom_par_A * sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts_mat(3, 1) = (1 - dom_par_A * sel_cof_A) * (1 - sel_cof_B);
  fts_mat(0, 2) = (1 - dom_par_A * sel_cof_A);
  fts_mat(1, 2) = (1 - dom_par_A * sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts_mat(2, 2) = (1 - sel_cof_A);
  fts_mat(3, 2) = (1 - sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts_mat(0, 3) = (1 - dom_par_A * sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts_mat(1, 3) = (1 - dom_par_A * sel_cof_A) * (1 - sel_cof_B);
  fts_mat(2, 3) = (1 - sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts_mat(3, 3) = (1 - sel_cof_A) * (1 - sel_cof_B);

  return fts_mat;
}

// Calculate the sampling probabilities
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

  // declare and initialise the sampling probabilities for the Wright-Fisher model
  arma::dcolvec prob = hap_frq;

  // calculate the haplotype frequencies after natural selection
  prob = hap_frq % (fts_mat * hap_frq) / arma::as_scalar(hap_frq.t() * fts_mat * hap_frq);

  // calculate the haplotype frequencies after genetic recombination
  prob = prob + eta * rec_rat * (prob(0) * prob(3) - prob(1) * prob(2));

  return prob;
}

// Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
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

  // declare the sampling probabilities
  arma::dcolvec prob = arma::zeros<arma::dcolvec>(4);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    // calculate the haplotype frequencies after genetic recombination
    prob = calculateSamplingProb_arma(hap_frq, fts_mat, rec_rat)

    // proceed the Wright-Fisher sampling
    IntegerVector hap_cnt(4);
    R::rmultinom(2 * pop_siz, prob.begin(), 4, hap_cnt.begin());
    hap_frq = as<arma::dcolvec>(hap_cnt) / 2 / pop_siz;

    frq_pth.col(k) = hap_frq;
  }

  return frq_pth;
}
/*************************/


/***** DiffusApprox ******/
// Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
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

  return frq_pth;
}
/*************************/


/***** MomentApprox ******/
// Approximate the first two moments of the two-locus Wright-Fisher model with selection using Monte Carlo simulation
// [[Rcpp::export]]
List approximateMoment_MonteCarlo_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat fts_mat = calculateFitnessMat_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B);

  // declare the mean vector and variance matrix of the Wright-Fisher model
  arma::dmat mu = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);

  // initialise the mean vector and variance matrix of the Wright-Fisher model
  mu.col(0) = int_frq;
  // sigma.slice(0) = arma::zeros<arma::dmat>(4, 4);

  // simulate the haplotype frequency trajectories according to the Wright-Fisher model
  arma::dcube frq_pth(4, arma::uword(lst_gen - int_gen) + 1, sim_num);
  for(arma::uword i = 1; i < sim_num; i++) {
    frq_pth.slice(i) = simulateWFM_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen);
  }

  // calculate the mean vector and variance matrix of the Wright-Fisher model
  arma::dmat frq_smp = frq_pth.col(k);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    mu.col(k) = arma::mean(frq_smp, 1);
    sigma.slice(k) = arma::cov(frq_smp.t(), frq_smp.t());
  }

  // return the approximations for the mean and variance of the Wright-Fisher model at each generation from int_gen to lst_gen
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}

//
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

// Approximate the first two moments of the two-locus Wright-Fisher model with selection using the extension of Lacerda & Seoighe (2014)
// [[Rcpp::export]]
List approximateMoment_Lacerda_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::dmat fts_mat) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the mean vector and variance matrix of the Wright-Fisher model
  arma::dmat mu = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);

  // initialise the mean vector and variance matrix of the Wright-Fisher model
  mu.col(0) = int_frq;
  // sigma.slice(0) = arma::zeros<arma::dmat>(4, 4);

  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    // approximate the mean vector of the Wright-Fisher model
    mu.col(k) = calculateSamplingProb_arma(mu.col(k - 1), fts_mat, rec_rat);

    // approximate the variance matrix of the Wright-Fisher model
    arma::dmat jacobian_mu = calculateJacobianMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    sigma.slice(k) = (arma::diagmat(mu.col(k)) - mu.col(k) * (mu.col(k)).t()) / 2 / pop_siz + grad_mu * sigma.slice(k - 1) * grad_mu.t();
  }

  // return the approximations for the mean and variance of the Wright-Fisher model at each generation from int_gen to lst_gen
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}


// Approximate the first two moments of the two-locus Wright-Fisher model with selection using the extension of Terhorst et al. (2015)
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


// Approximate the first two moments of the two-locus Wright-Fisher model with selection using the extension of Paris et al. (2019)
// [[Rcpp::export]]
List approximateMoment_Paris_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::dmat fts_mat) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the mean vector and variance matrix of the Wright-Fisher model
  arma::dmat mu = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);

  // initialise the mean vector and variance matrix of the Wright-Fisher model
  mu.col(0) = int_frq;
  // sigma.slice(0) = arma::zeros<arma::dmat>(4, 4);

  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    // approximate the mean vector of the Wright-Fisher model
    mu.col(k) = calculateSamplingProb_arma(mu.col(k - 1), fts_mat, rec_rat);

    // approximate the variance matrix of the Wright-Fisher model
    arma::dmat jacobian_mu = calculateJacobianMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    sigma.slice(k) = (arma::diagmat(mu.col(k)) - mu.col(k) * (mu.col(k)).t()) / 2 / pop_siz + (1 - 1.0 / 2 / pop_siz) * grad_mu * sigma.slice(k - 1) * grad_mu.t();
  }

  // return the approximations for the mean and variance of the Wright-Fisher model at each generation from int_gen to lst_gen
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}
/*************************/


/****** ParamApprox ******/
// Approximate the two-locus Wright-Fisher model with selection using the normal distribution
// [[Rcpp::export]]
List approximatWFM_norm_arma(const arma::dmat& mean, const arma::dcube& variance) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat mu = mean;
  arma::dcube sigma = variance;

  // return the parameters of the normal distribution
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}

// Simulate the two-locus Wright-Fisher model with selection using the normal approximation
// [[Rcpp::export]]
arma::dcube simulateWFM_norm_arma(const arma::dmat& mean, const arma::dcube& variance, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dcube frq_pth(4, arma::uword(lst_gen - int_gen) + 1, sim_num);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    frq_pth.col(k) = arma::mvnrnd(mean.col(k), variance.slice(k), sim_num);
  }

  //
  return frq_pth;
}


// Calculate the inverse additive logistic transformation
// [[Rcpp::export]]
arma::dmat calculateInvALT_arma(const arma::dmat& hap_frq){
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat phi = arma::zeros<arma::dmat>(3, hap_frq.n_cols);
  phi.row(0) = log(hap_frq.row(0) / hap_frq.row(3));
  phi.row(1) = log(hap_frq.row(1) / hap_frq.row(3));
  phi.row(2) = log(hap_frq.row(2) / hap_frq.row(3));

  // return the inverse additive logistic transformation
  return phi;
}

// Calculate the Jacobian matrix over the inverse additive logistic transformation
// [[Rcpp::export]]
arma::dcube calculateJacobianInvALT_arma(const arma::dmat& hap_frq){
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dcube jacobian_mat = arma::zeros<arma::dcube>(3, 4, hap_frq.n_cols);
  jacobian_mat.tube(0, 0) = 1.0 / hap_frq.row(0);
  jacobian_mat.tube(1, 1) = 1.0 / hap_frq.row(1);
  jacobian_mat.tube(2, 2) = 1.0 / hap_frq.row(2);
  jacobian_mat.tube(0, 3) = -1.0 / hap_frq.row(3);
  jacobian_mat.tube(1, 3) = -1.0 / hap_frq.row(3);
  jacobian_mat.tube(2, 3) = -1.0 / hap_frq.row(3);

  // return the Jacobian matrix over the inverse additive logistic transformation
  return jacobian_mat;
}

// Approximate the two-locus Wright-Fisher model with selection using the logistic normal distribution
// [[Rcpp::export]]
List approximatWFM_logisticnorm_arma(const arma::dmat& mean, const arma::dcube& variance, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare and calculate the location and squared scale of the logistic normal distribution
  arma::dmat location = calculateInvALT_arma(mean);

  // declare and calculate the location and squared scale of the logistic normal distribution
  arma::dcube squared_scale = arma::zeros<arma::dcube>(3, 3, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube jacobian_phi = calculateJacobianInvALT_arma(mean);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    squared_scale.slice(k) = jacobian_phi.slice(k) * variance.slice(k) * jacobian_phi.slice(k).t();
  }

  // return the parameters of the logistic normal distribution
  return List::create(Named("location", location),
                      Named("squared_scale", squared_scale));
}

// Calculate the additive logistic transformation
// [[Rcpp::export]]
arma::dmat calculateALT_arma(const arma::dmat& phi_frq){
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat hap_frq = arma::zeros<arma::dmat>(4, phi_frq.n_cols);
  hap_frq.row(0) = exp(phi_frq.row(0));
  hap_frq.row(1) = exp(phi_frq.row(1));
  hap_frq.row(2) = exp(phi_frq.row(2));
  hap_frq.row(3) = arma::ones<arma::drowvec>(phi_frq.n_cols);
  hap_frq = arma::normalise(hap_frq, p = 1, dim = 0);

  // return the additive logistic transformation
  return hap_frq;
}

// Simulate the two-locus Wright-Fisher model with selection using the logistic normal approximation
arma::dcube simulateWFM_logisticnorm_arma(const arma::dmat& location, const arma::dcube& squared_scale, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dcube phi_pth = arma::zeros<arma::dcube>(3, arma::uword(lst_gen - int_gen) + 1, sim_num);
  arma::dcube frq_pth = arma::zeros<arma::dcube>(4, arma::uword(lst_gen - int_gen) + 1, sim_num);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    phi_pth.col(k) = arma::mvnrnd(location.col(k), squared_scale.slice(k), sim_num);
    frq_pth.col(k) = calculateALT_arma(phi.col(k));
  }

  //
  return frq_pth;
}

// Approximate the first two moments of the logistic normal approximation of the two-locus Wright-Fisher model with selection
List approximateMoment_logisticnorm_arma(const arma::dmat& location, const arma::dcube& squared_scale, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dcube frq_pth = simulateWFM_logisticnorm_arma(location, squared_scale, int_gen, lst_gen, sim_num);

  arma::dmat mean = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube variance = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);

  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    arma::dmat frq_smp = frq_pth.col(k);
    mean.col(k) = arma::mean(frq_smp, 1);
    variance.slice(k) = arma::cov(frq_smp.t(), frq_smp.t());
  }

  // return the approximations for the mean and variance of the logistic normal approximation at each generation from int_gen to lst_gen
  return List::create(Named("mean", mean),
                      Named("variance", variance));
}


// Approximate the two-locus Wright-Fisher model with selection using the hierarchical beta distribution
// [[Rcpp::export]]
List approximatWFM_hierarchicalbeta_arma(const arma::dmat& mean, const arma::dcube& variance, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat m(3, arma::uword(lst_gen - int_gen) + 1);
  m.row(0) = mean.row(0) + mean.row(1);
  m.row(1) = mean.row(0) / (mean.row(0) + mean.row(1));
  m.row(2) = mean.row(2) / (mean.row(2) + mean.row(3));

  arma::dcube V(3, 3, arma::uword(lst_gen - int_gen) + 1);

  arma::dmat alpha(3, arma::uword(lst_gen - int_gen) + 1);
  arma::dmat beta(3, arma::uword(lst_gen - int_gen) + 1);

  for(arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    arma::dcolvec col_m(3);
    arma::dcolvec diag_V(3);

    diag_V(0) = (variance(0, 0, t) - variance(1, 1, t)) / (col_m(1) * col_m(1) - (1 - col_m(1)) * (1 - col_m(1)));
    diag_V(1) = (variance(0, 0, t) + variance(1, 1, t) - (col_m(1) * col_m(1) + (1 - col_m(1)) * (1 - col_m(1))) * diag_V(0)) / 2 / (col_m(0) * col_m(0) + diag_V(0));
    diag_V(2) = (variance(2, 2, t) + variance(3, 3, t) - (col_m(2) * col_m(2) + (1 - col_m(2)) * (1 - col_m(2))) * diag_V(0)) / 2 / ((1 - col_m(0)) * (1 - col_m(0)) + diag_V(0));
    m.col(t) = col_m;
    V.slice(t) = arma::diagmat(diag_V);
    alpha.col(t) = (col_m % (1 - col_m) / diag_V - 1 ) % col_m;
    beta.col(t) = (col_m % (1 - col_m) / diag_V - 1 ) % (1 - col_m);
  }


  // return the parameters of the hierarchical beta distribution
  return List::create(Named("alpha", alpha),
                      Named("beta", beta));
}


// Simulate the two-locus Wright-Fisher model with selection using the hierarchical beta approximation
arma::dcube simulateWFM_hierarchicalbeta_arma(const arma::dmat& alpha, const arma::dmat& beta, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dcube Z(3, arma::uword(lst_gen - int_gen) + 1, sim_num);
  for(arma::uword t = 0; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    arma::dcolvec alpha_t = alpha.col(t);
    arma::dcolvec beta_t = beta.col(t);
    arma::dmat Z_t(3, sim_num);
    arma::drowvec Z_0_t = Rcpp::rbeta(sim_num, alpha_t(0), beta_t(0));
    arma::drowvec Z_1_t = Rcpp::rbeta(sim_num, alpha_t(1), beta_t(1));
    arma::drowvec Z_2_t = Rcpp::rbeta(sim_num, alpha_t(2), beta_t(2));
    Z_t.row(0) = Z_0_t;
    Z_t.row(1) = Z_1_t;
    Z_t.row(2) = Z_2_t;
    Z.col(t) = Z_t;
  }
  arma::dcube frq_pth(4, arma::uword(lst_gen - int_gen) + 1, sim_num);
  frq_pth.row(0) = Z.row(0) % Z.row(1);
  frq_pth.row(1) = Z.row(0) - Z.row(0) % Z.row(1);
  frq_pth.row(2) = Z.row(2) - Z.row(0) % Z.row(2);
  arma::dcube ones_cube = arma::ones<arma::dcube>(3, arma::uword(lst_gen - int_gen) + 1, sim_num);
  frq_pth.row(3) = (ones_cube.row(0) - Z.row(0)) % (ones_cube.row(0) - Z.row(2));

  return frq_pth;
}

// Approximate the first two moments of the hierarchical beta approximation of the two-locus Wright-Fisher model with selection
List approximateMoment_hierarchicalbeta_arma(const arma::dmat& alpha, const arma::dmat& beta, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat mean = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube variance = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);

  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    arma::dcolvec m = alpha.col(k) / (alpha.col(k) + beta.col(k));
    arma::dcolvec V = alpha.col(k) % beta.col(k) / (alpha.col(k) + beta.col(k)) / (alpha.col(k) + beta.col(k)) / (alpha.col(k) + beta.col(k) + 1);

    arma::dcolvec mu = arma::zeros<arma::dcolvec>(4);
    mu(0) = m(0) * m(1);
    mu(1) = m(0) * (1.0 - m(1));
    mu(2) = (1.0 - m(0)) * m(2);
    mu(3) = (1.0 - m(0)) * (1.0 - m(2));

    arma::dmat sigma = arma::zeros<arma::dmat>(4, 4);
    sigma(0, 0) = V(0) * m(1) * m(1) + (V(0) + m(0) * m(0)) * V(1);
    sigma(0, 1) = -(V(0) + m(0) * m(0)) * V(1) + V(0) * m(1) * (1 - m(1));
    sigma(0, 2) = -V(0) * m(1) * mt(2);
    sigma(0, 3) = -V(0) * m(1) * (1 - m(2));
    sigma(1, 0) = sigma(0, 1);
    sigma(1, 1) = V(0) * (1 - m(1)) * (1 - m(1)) + (V(0) + m(0) * m(0)) * V(1);
    sigma(1, 2) = -V(0) * (1 - m(1)) * m(2);
    sigma(1, 3) = -V(0) * (1 - m(1)) * (1 - m(2));
    sigma(2, 0) = sigma(0, 2);
    sigma(2, 1) = sigma(1, 2);
    sigma(2, 2) = V(0) * m(2) * m(2) + (V(0) + (1 - m(0)) * (1 - m(0))) * V(2);
    sigma(2, 3) = -(V(0) + (1 - m(0)) * (1 - m(0))) * V(2) + V(0) * m(2) * (1 - m(2));
    sigma(3, 0) = sigma(0, 3);
    sigma(3, 1) = sigma(1, 3);
    sigma(3, 2) = sigma(2, 3);
    sigma(3, 3) = V(0) * (1 - m(2)) * (1 - m(2)) + (V(0) + (1 - m(0)) * (1 - m(0))) * V(2);

    mean.col(k) = mu;
    variance.slice(k) = sigma;
  }

  // return the approximations for the mean and variance of the hierarchical beta approximation at each generation from int_gen to lst_gen
  return List::create(Named("mean", mean),
                      Named("variance", variance));
}

/*************************/
