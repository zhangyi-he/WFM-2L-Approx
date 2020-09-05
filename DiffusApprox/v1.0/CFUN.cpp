// Numerical simulation of the two-locus Wright-Fisher diffusion with application to approximating transition probability densities
// Zhangyi He, Mark Beaumont and Feng Yu

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
// Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::drowvec simulateWFM_1L_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const double& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the fitness
  arma::dcolvec fts_vec = arma::ones<arma::dcolvec>(3);
  // fts(0) = 1.0;
  fts_vec(1) = 1.0 - sel_cof * dom_par;
  fts_vec(2) = 1.0 - sel_cof;

  // declare the mutant allele frequency trajectory
  arma::drowvec frq_pth = arma::zeros<arma::drowvec>(arma::uword(lst_gen - int_gen) + 1);

  // initialise the mutant allele frequency in generation 0
  frq_pth(0) = int_frq;

  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    // calculate the sampling probability
    arma::dcolvec gen_frq = arma::zeros<arma::dcolvec>(3);
    gen_frq(0) = frq_pth(t - 1) * frq_pth(t - 1);
    gen_frq(1) = 2 * frq_pth(t - 1) * (1 - frq_pth(t - 1));
    gen_frq(2) = (1 - frq_pth(t - 1)) * (1 - frq_pth(t - 1));
    gen_frq = (fts_vec % gen_frq) / arma::accu(fts_vec % gen_frq);
    double prob = gen_frq(0) + gen_frq(1) / 2;

    // proceed the Wright-Fisher sampling
    frq_pth(t) = R::rbinom(2 * pop_siz, prob) / 2 / pop_siz;
  }

  // return the mutant allele frequency trajectories under the Wright-Fisher model
  return frq_pth;
}

// Calculate the fitness matrix for the two-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::dmat calculateFitnessMat_2L_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B) {
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

// Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::dmat simulateWFM_2L_arma(const arma::dmat& fts_mat, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
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


/********** WFD **********/
// Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
// [[Rcpp::export]]
arma::drowvec simulateWFD_1L_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const double& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // rescale the selection coefficient
  double scl_sel_cof = 2 * pop_siz * sel_cof;

  // calculate delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  // generate delta W
  arma::drowvec dW = pow(dt, 0.5) * arma::randn<arma::drowvec>(arma::uword(lst_gen - int_gen) * ptn_num);

  // declare the mutant allele frequency trajectory
  arma::drowvec frq_pth = arma::zeros<arma::drowvec>(arma::uword(lst_gen - int_gen) * ptn_num + 1);

  // initialise the mutant allele frequency in generation 0
  frq_pth(0) = int_frq;

  // simulate the mutant allele frequency trajectory
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) * ptn_num + 1; t++) {
    // calculate the drift coefficient
    double mu = scl_sel_cof * frq_pth(t - 1) * (1 - frq_pth(t - 1)) * ((1 - dom_par) - (1 - 2 * dom_par) * frq_pth(t - 1));

    // calculate the diffusion coefficient
    double sigma = pow(frq_pth(t - 1) * (1 - frq_pth(t - 1)), 0.5);

    // proceed the Euler-Maruyama scheme
    frq_pth(t) = frq_pth(t - 1) + mu * dt + sigma * dW(t - 1);

    // remove the noise from the numerical techniques
    if (frq_pth(t) < 0) {
      frq_pth(t) = 0;
    }
    if (frq_pth(t) > 1) {
      frq_pth(t) = 1;
    }
  }

  // return the mutant allele frequency trajectory under the Wright-Fisher diffusion
  return frq_pth;
}

// Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
// [[Rcpp::export]]
arma::dmat simulateWFD_2L_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
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
  arma::dmat frq_pth = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) * ptn_num + 1);

  // initialise the haplotype frequencies in generation 0
  frq_pth.col(0) = int_frq;

  // simulate the haplotype frequency trajectories
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


/********* ECDF **********/
// Generate the sample mutant allele frequency trajectories according to the one-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::dmat generateSample_WFM_1L_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const double& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat frq_smp = arma::zeros<arma::dmat>(sim_num, arma::uword(lst_gen - int_gen) + 1);
  for (arma::uword i = 0; i < sim_num; i++) {
    frq_smp.row(i) = simulateWFM_1L_arma(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen);
  }

  return frq_smp;
}

// Generate the sample mutant allele frequency trajectories according to the one-locus Wright-Fisher diffusion with selection
// [[Rcpp::export]]
arma::dmat generateSample_WFD_1L_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const double& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat frq_smp = arma::zeros<arma::dmat>(sim_num, arma::uword(lst_gen - int_gen) * ptn_num + 1);
  for (arma::uword i = 0; i < sim_num; i++) {
    arma::dmat frq_pth = simulateWFD_1L_arma(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num);
    arma::ucolvec sel_gen = arma::linspace<arma::ucolvec>(0, arma::uword(lst_gen - int_gen), arma::uword(lst_gen - int_gen) + 1);
    frq_smp.row(i) = frq_pth.cols(sel_gen);
  }

  return frq_smp;
}

// Generate the sample haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::dcube generateSample_WFM_2L_arma(const arma::dmat& fts_mat, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dcube frq_smp = arma::zeros<arma::dcube>(sim_num, 4, arma::uword(lst_gen - int_gen) + 1);
  for (arma::uword i = 0; i < sim_num; i++) {
    frq_smp.row(i) = simulateWFM_2L_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen);
  }

  return frq_smp;
}

// Generate the sample haplotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection
// [[Rcpp::export]]
arma::dcube generateSample_WFD_2L_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dcube frq_smp = arma::zeros<arma::dcube>(sim_num, 4, arma::uword(lst_gen - int_gen) + 1);
  for (arma::uword i = 0; i < sim_num; i++) {
    arma::dmat frq_pth = simulateWFD_2L_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num);
    arma::ucolvec sel_gen = arma::linspace<arma::ucolvec>(0, arma::uword(lst_gen - int_gen), arma::uword(lst_gen - int_gen) + 1);
    frq_smp.row(i) = frq_pth.cols(sel_gen);
  }

  return frq_smp;
}

// Generate the grids for empirical cumulative distribution function
// [[Rcpp::export]]
arma::dmat generateGrid_2L_arma(const arma::uword& grd_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat frq_grd = arma::zeros<arma::dmat>(1, 4);
  for (arma::uword i = 0; i < grd_num + 1; i++) {
    for (arma::uword j = 0; j < grd_num + 1 - i; j++) {
      for (arma::uword k = 0; k < grd_num + 1 - i - j; k++) {
        arma::drowvec frq = {double(i), double(j), double(k), double(grd_num - i - j - k)};
        frq_grd.insert_rows(0, 1);
        frq_grd.row(0) = frq / double(grd_num);
      }
    }
  }
  frq_grd.shed_row(frq_grd.n_rows - 1);

  return frq_grd;
}
/*************************/


/****** GuidedProc *******/
// Simulate the mutant allele frequency trajectory according to the guided process of Fearnhead (2008) using the Euler-Maruyama method



// Simulate the mutant allele frequency trajectory according to the guided process of He et al. (2020) using the Euler-Maruyama method



/*************************/


/******* TransPDF ********/
// Approximate the transition probability density of the two-locus Wright-Fisher model with selection using Monte Carlo integration



// Calculate the importance weight



// Approximate the transition probability density of the two-locus Wright-Fisher model with selection using Monte Carlo integration with importance sampling



/*************************/
