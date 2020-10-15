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
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
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
List approximateMoment_MonteCarlo_arma(const arma::dmat& fts_mat, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // simulate the haplotype frequency trajectories under the Wright-Fisher model
  arma::dcube frq_pth(4, arma::uword(lst_gen - int_gen) + 1, sim_num);
  for(arma::uword i = 0; i < sim_num; i++) {
    cout << "iteration: " << i + 1 << endl;
    frq_pth.slice(i) = simulateWFM_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen);
  }

  // declare the mean vector and variance matrix
  arma::dmat mu = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen));
  arma::dcube sigma = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen));

  // calculate the mean vector and variance matrix
  for(arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    arma::dmat frq_smp = frq_pth.col(k + 1);
    mu.col(k) = arma::mean(frq_smp, 1);
    sigma.slice(k) = arma::cov(frq_smp.t(), frq_smp.t());
  }

  // return the Monte Carlo approximations for the mean vector and variance matrix of the Wright-Fisher model
  return List::create(Named("mean", mu.shed_col(0)),
                      Named("variance", sigma.shed_slice(0)));
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
arma::dmat calculateJacobianMean_arma(const arma::dcolvec& hap_frq, const arma::dmat& fts_mat, const double& rec_rat){
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
  for(arma::uword i = 0; i < 4; i++) {
    Q.col(i) = q + pow(-1, i) * eta * q(3 - i);
  }

  // calculate the Jacobian matrix over the mean vector
  arma::dmat Jacobian_mu = (arma::eye<arma::dmat>(4, 4) + rec_rat * Q) * (arma::diagmat(hap_frq) * (fts_mat / fts_val - 2 * (fts_vec / fts_val) * (fts_vec.t() / fts_val)) + arma::diagmat(fts_vec / fts_val));

  // return the Jacobian matrix over the mean vector
  return Jacobian_mu;
}

// Calculate the Hessian matrix over the mean vector for the Wright-Fisher model
// [[Rcpp::export]]
arma::cube calculateHessianMean_arma(const arma::dcolvec& hap_frq, const arma::dmat& fts_mat, const double& rec_rat){
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the marginal and mean fitness
  arma::dcolvec fts_vec = fts_mat * hap_frq;
  double fts_val = arma::as_scalar(hap_frq.t() * fts_mat * hap_frq);

  //
  arma::dcolvec q = (hap_frq % fts_vec) / fts_val;

  // calculate the Hessian matrix over the mean vector
  arma::dmat Jacobian_q = arma::zeros<arma::dmat>(4, 4);
  for(arma::uword i = 0; i < 4; i++) {
    for(arma::uword j = 0; j < 4; j++) {
      Jacobian_q(i, j) = (fts_mat(i, j) * hap_frq(i) + fts_vec(i) * ((i == j) ? 1.0 : 0.0)) / fts_val - 2 * fts_vec(i) * fts_vec(j) * hap_frq(i) / fts_val / fts_val;
    }
  }
  arma::cube Hessian_q = arma::zeros<arma::dcube>(4, 4, 4);
  for(arma::uword i = 0; i < 4; i++) {
    for(arma::uword j = 0; j < 4; j++) {
      for(arma::uword k = 0; k < 4; k++) {
        Hessian_q(i, j, k) = fts_mat(k, j) * (((k == i) ? 1.0 : 0.0) + ((i == j) ? 1.0 : 0.0)) / fts_val - 2 * (fts_mat(i, j) * fts_vec(k) * hap_frq(i) + fts_mat(i, k) * fts_vec(j) * hap_frq(i) + fts_mat(j, k) * fts_vec(i) * hap_frq(i) + fts_vec(j) * fts_vec(k) * ((j == i) ? 1.0 : 0.0) + ((i == k) ? 1.0 : 0.0))) / fts_val / fts_val + 8 * (fts_vec(i) * fts_vec(j) * fts_vec(k) * hap_frq(i)) / fts_val / fts_val / fts_val;
      }
    }
  }
  arma::cube Hessian_mu = arma::zeros<arma::dcube>(4, 4, 4);
  for(arma::uword i = 0; i < 4; i++) {
    for(arma::uword j = 0; j < 4; j++) {
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
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    arma::dmat Jacobian_mu = calculateJacobianMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    mu.col(k) = calculateMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    sigma.slice(k) = (1.0 / 2 / pop_siz) * (arma::diagmat(mu.col(k)) - mu.col(k) * mu.col(k).t()) + Jacobian_mu * sigma.slice(k - 1) * Jacobian_mu.t();
  }

  // return the approximations for the mean vector and variance matrix of the Wright-Fisher model
  return List::create(Named("mean", mu.shed_col(0)),
                      Named("variance", sigma.shed_slice(0)));
}



// Approximate the first two moments of the Wright-Fisher model using the extension of Terhorst et al. (2015)
// [[Rcpp::export]]
List approximateMoment_Terhorst_arma(const arma::dmat& fts_mat, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the mean vector and variance matrix
  arma::dmat x = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dmat epsilon = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dmat mu = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube epsilon_2ndOrder = arma::zeros<arma::dcube>(4, 4, arma::uword(lst_gen - int_gen) + 1);

  // calculate the mean vector and variance matrix
  x.col(0) = int_frq;
  // epsilon.col(0) = arma::zeros<arma::dcolvec>(4);
  mu.col(0) = int_frq;
  // sigma.slice(0) = arma::zeros<arma::dmat>(4, 4);
  // epsilon_2ndOrder.slice(0) = arma::zeros<arma::dmat>(4, 4);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    arma::dmat Jacobian_mu = calculateJacobianMean_arma(x.col(k - 1), fts_mat, rec_rat);
    arma::dcube Hessian_mu = calculateHessianMean_arma(x.col(k - 1), fts_mat, rec_rat);
    arma::dcolvec SecndOrder_term = arma::zeros<arma::dcolvec>(4);
    arma::dmat epsilon_2ndOrder_k_minus_1 = epsilon_2ndOrder.slice(k - 1);
    for(arma::uword i = 1; i < 4; i++) {
      arma::dmat Hessian_mu_i = Hessian_mu.row(i);
      SecndOrder_term(i) = 0.5 * arma::as_scalar(sum(sum(Hessian_mu_i % epsilon_2ndOrder_k_minus_1, 0), 1));
    }

    x.col(k) = calculateMean_arma(x.col(k - 1), fts_mat, rec_rat);
    epsilon.col(k) = Jacobian_mu * epsilon.col(k - 1) + SecndOrder_term;
    mu.col(k) = x.col(k) + epsilon.col(k);
    sigma.slice(k) = 1.0 / 2 / pop_siz * arma::diagmat(x.col(k) + Jacobian_mu * epsilon.col(k - 1)) - 1.0 / 2 / pop_siz * (x.col(k) * x.col(k).t() + Jacobian_mu * epsilon.col(k - 1) * x.col(k).t() + x.col(k) * (Jacobian_mu * epsilon.col(k - 1)).t() + Jacobian_mu * epsilon_2ndOrder_k_minus_1 * Jacobian_mu) + (1 - 1.0 / 2 / pop_siz) * Jacobian_mu * sigma.slice(k - 1) * Jacobian_mu.t();
    epsilon_2ndOrder.slice(k) = sigma.slice(k) + epsilon.col(k) * epsilon.col(k).t();
  }

  // return the approximations for the mean vector and variance matrix of the Wright-Fisher model
  return List::create(Named("mean", mu.shed_col(0)),
                      Named("variance", sigma.shed_slice(0)));
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
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    arma::dmat Jacobian_mu = calculateJacobianMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    mu.col(k) = calculateMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    sigma.slice(k) = (1.0 / 2 / pop_siz) * (arma::diagmat(mu.col(k)) - mu.col(k) * mu.col(k).t()) + (1 - 1.0 / 2 / pop_siz) * Jacobian_mu * sigma.slice(k - 1) * Jacobian_mu.t();
  }

  // return the approximations for the mean vector and variance matrix of the Wright-Fisher model
  return List::create(Named("mean", mu.shed_col(0)),
                      Named("variance", sigma.shed_slice(0)));
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
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    arma::dmat Jacobian_mu = calculateJacobianMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    mu.col(k) = calculateMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    arma::dcube Hessian_mu = calculateHessianMean_arma(mu.col(k - 1), fts_mat, rec_rat);
    arma::dcolvec SecndOrder_term = arma::zeros<arma::dcolvec>(4);
    arma::dmat sigma_k_minus_1 = sigma.slice(k - 1);
    for(arma::uword i = 1; i < 4; i++) {
      arma::dmat Hessian_mu_i = Hessian_mu.row(i);

      SecndOrder_term(i) = 0.5 * arma::as_scalar(sum(sum(Hessian_mu_i % sigma_k_minus_1, 0), 1));
    }
    sigma.slice(k) = (arma::diagmat(mu.col(k) + SecndOrder_term) - (mu.col(k) + SecndOrder_term) * (mu.col(k) + SecndOrder_term).t()) / 2 / pop_siz + (1 - 1.0 / 2 / pop_siz) * Jacobian_mu * sigma.slice(k - 1) * Jacobian_mu.t();
  }

  // return the approximations for the mean vector and variance matrix of the Wright-Fisher model
  return List::create(Named("mean", mu.shed_col(0)),
                      Named("variance", sigma.shed_slice(0)));
}
/*************************/


/****** ParamApprox ******/
// Generate the samples under the normal approximation of the Wright-Fisher model
// [[Rcpp::export]]
arma::dcube simulateWFM_norm_arma(const arma::dmat& mean, const arma::dcube& variance, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // simulate the haplotype frequency trajectories under the normal approximation
  arma::dcube frq_pth = arma::zeros<arma::dcube>(4, arma::uword(lst_gen - int_gen), sim_num);
  for(arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
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
arma::dcube calculateJacobianLocation_arma(const arma::dmat& hap_frq){
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
  for(arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    scalesq.slice(k) = Jacobian_location.slice(k) * variance.slice(k) * Jacobian_location.slice(k).t();
  }

  // return the parameters of the logistic normal approximation
  return List::create(Named("location", location),
                      Named("scalesq", scalesq));
}

// Calculate the parameters for the logistic normal approximation of the Wright-Fisher model
// [[Rcpp::export]]
arma::dmat calculateALT_arma(const arma::dmat& phi_frq){
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
  for(arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    arma::dmat phi_frq = arma::mvnrnd(location.col(k), scalesq.slice(k), sim_num);
    frq_pth.col(k) = calculateALT_arma(phi_frq);
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
  for(arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    arma::dmat frq_smp = frq_pth.col(k);
    mean.col(k) = arma::mean(frq_smp.t(), 0);
    variance.slice(k) = arma::cov(frq_smp.t());
  }

  // return the approximations for the mean vector and variance matrix of the logistic normal approximation
  return List::create(Named("mean", mean),
                      Named("variance", variance));
}

// Calculate the parameters for the hierarchical beta approximation of the Wright-Fisher model
// [[Rcpp::export]]
List calculateParam_HierarchicalBeta_arma(const arma::dmat& mean, const arma::dcube& variance, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the parameters of the hierarchical beta approximation
  arma::dmat alpha = arma::zeros<arma::dmat>(3, arma::uword(lst_gen - int_gen));
  arma::dmat beta = arma::zeros<arma::dmat>(3, arma::uword(lst_gen - int_gen));
  for(arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
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

// Generate the samples under the hierarchical beta approximation of the Wright-Fisher model
// [[Rcpp::export]]
arma::dcube generateSample_HierarchicalBeta_arma(const arma::dmat& alpha, const arma::dmat& beta, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // simulate the haplotype frequency trajectories
  arma::dcube frq_pth(4, arma::uword(lst_gen - int_gen), sim_num);
  for(arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    arma::dmat psi_frq = arma::zeros<arma::dmat>(3, sim_num);
    psi_frq.row(0) = as<arma::drowvec>(Rcpp::rbeta(sim_num, alpha(0, k), beta(0, k)));
    psi_frq.row(1) = as<arma::drowvec>(Rcpp::rbeta(sim_num, alpha(1, k), beta(1, k)));
    psi_frq.row(2) = as<arma::drowvec>(Rcpp::rbeta(sim_num, alpha(2, k), beta(2, k)));

    frq_pth.tube(0, k) = psi_frq.row(0) % psi_frq.row(1);
    frq_pth.tube(1, k) = psi_frq.row(0) % (1 - psi_frq.row(1));
    frq_pth.tube(2, k) = (1 - psi_frq.row(0)) % psi_frq.row(2);
    frq_pth.tube(3, k) = (1 - psi_frq.row(0)) % (1 - psi_frq.row(2));
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
  for(arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    // calculate the mean and variance of the beta distributions
    arma::dcolvec m = alpha.col(k) / (alpha.col(k) + beta.col(k));
    arma::dcolvec V = alpha.col(k) % beta.col(k) / (alpha.col(k) + beta.col(k)) / (alpha.col(k) + beta.col(k)) / (alpha.col(k) + beta.col(k) + 1);

    // calculate the mean vector
    mean(0, k) = m(0) * m(1);
    mean(1, k) = m(0) * (1.0 - m(1));
    mean(2, k) = (1.0 - m(0)) * m(2);
    mean(3, k) = (1.0 - m(0)) * (1.0 - m(2));

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
/*************************/


/********* ECDF **********/
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
        frq_grd.row(0) = hap_frq /grd_num;
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
  for(int j = 0; j < 4; j++){
    frq_grd(_, j) = rgamma(grd_num, 1.0, 1.0);
  }
  for(int i = 0; i < grd_num; i++){
    frq_grd(i, _) = frq_grd(i, _) / sum(frq_grd(i, _));
  }

  return as<arma::dmat>(frq_grd);
}

// calculate the empirical cumulative distribution function
// [[Rcpp::export]]
arma::dmat calculateECDF_2L_arma(const arma::dcube& frq_pth, const arma::uword& sim_num, const arma::dmat& frq_grd) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat ECDF = arma::zeros<arma::dmat>(frq_grd.n_rows, frq_pth.n_cols - 1);
  for(arma::uword k = 0; k < frq_pth.n_cols - 1; k++) {
    arma::dmat frq_smp = frq_pth.col(k + 1);
    arma::dcolvec cdf = arma::zeros<arma::dcolvec>(frq_grd.n_rows);

    //
    for(arma::uword i = 0; i < frq_grd.n_rows; i++) {
      arma::dcolvec cond = arma::zeros<arma::dcolvec>(sim_num);
      cond = (frq_smp.col(0) - frq_grd(i, 0)) / abs(frq_smp.col(0) - frq_grd(i, 0)) + (frq_smp.col(1) - frq_grd(i, 1)) / abs(frq_smp.col(1) - frq_grd(i, 1)) + (frq_smp.col(2) - frq_grd(i, 2)) / abs(frq_smp.col(2) - frq_grd(i, 2));
      for(arma::uword j = 0; j < sim_num; j++) {
        ECDF(i, k) = ECDF(i, k) + ((cond(j) == -3) ? 1.0 : 0.0);
      }
    }
  }

  return ECDF;
}
/*************************/
