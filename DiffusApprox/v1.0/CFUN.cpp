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
arma::drowvec simulateOLWFMS_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const double& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // declare the mutant allele frequency trajectory
  arma::drowvec ale_frq_pth(arma::uword(lst_gen - int_gen) + 1);
  
  // initialise the mutant allele frequency in generation 0
  ale_frq_pth(0) = int_frq;
  
  // declare the fitness
  arma::dcolvec fts = arma::zeros<arma::dcolvec>(3);
  fts(0) = 1.0;
  fts(1) = 1.0 - sel_cof * dom_par;
  fts(2) = 1.0 - sel_cof;
  
  // declare and initialise the genotype frequencies during a single generation of the life cycle
  arma::dcolvec gen_frq = arma::zeros<arma::dcolvec>(3);
  
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    // random union of gametes
    gen_frq(0) = ale_frq_pth(t - 1) * ale_frq_pth(t - 1);
    gen_frq(1) = 2 * ale_frq_pth(t - 1) * (1 - ale_frq_pth(t - 1));
    gen_frq(2) = (1 - ale_frq_pth(t - 1)) * (1 - ale_frq_pth(t - 1));
    
    // viability selection
    gen_frq = (fts % gen_frq) / arma::accu(fts % gen_frq);
    
    // meiosis (calculate the sampling probability)
    double prob = gen_frq(0) + gen_frq(1) / 2;
    
    // reproduction (the Wright-Fisher sampling)
    ale_frq_pth(t) = R::rbinom(2 * pop_siz, prob) / 2 / pop_siz;
  }
  
  return ale_frq_pth;
}

// Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::dmat simulateTLWFMS_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // declare the haplotype frequency trajectories
  arma::dmat hap_frq_pth(4, arma::uword(lst_gen - int_gen) + 1);
  
  // initialise the haplotype frequencies in generation 0
  hap_frq_pth.col(0) = int_frq;
  
  // declare the fitness
  arma::dmat fts(4, 4);
  fts(0, 0) = 1;
  fts(1, 0) = (1 - dom_par_B * sel_cof_B);
  fts(2, 0) = (1 - dom_par_A * sel_cof_A);
  fts(3, 0) = (1 - dom_par_A * sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts(0, 1) = (1 - dom_par_B * sel_cof_B);
  fts(1, 1) = (1 - sel_cof_B);
  fts(2, 1) = (1 - dom_par_A * sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts(3, 1) = (1 - dom_par_A * sel_cof_A) * (1 - sel_cof_B);
  fts(0, 2) = (1 - dom_par_A * sel_cof_A);
  fts(1, 2) = (1 - dom_par_A * sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts(2, 2) = (1 - sel_cof_A);
  fts(3, 2) = (1 - sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts(0, 3) = (1 - dom_par_A * sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts(1, 3) = (1 - dom_par_A * sel_cof_A) * (1 - sel_cof_B);
  fts(2, 3) = (1 - sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts(3, 3) = (1 - sel_cof_A) * (1 - sel_cof_B);
  
  // declare and initialise the haplotype frequencies during a single generation of the life cycle
  arma::dcolvec hap_frq = int_frq;
  
  // declare and initialise the genotype frequencies during a single generation of the life cycle
  arma::dmat gen_frq = hap_frq * hap_frq.t();
  
  // declare eta
  arma::dcolvec eta(4);
  eta(0) = -1;
  eta(1) = 1;
  eta(2) = 1;
  eta(3) = -1;
  
  for(arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    // random union of gametes
    gen_frq = hap_frq * hap_frq.t();
    
    // viability selection
    gen_frq = (fts % gen_frq) / arma::as_scalar(sum(sum(fts % gen_frq, 0), 1));
    
    // meiosis (calculate the sampling probability)
    arma::dcolvec prob = arma::zeros<arma::dcolvec>(4);
    for(arma::uword i = 0; i < 4; i++) {
      prob(i) = (sum(gen_frq.row(i)) + sum(gen_frq.col(i))) / 2;
    }
    prob = prob + eta * rec_rat * (prob(0) * prob(3) - prob(1) * prob(2));
    
    // reproduction (the Wright-Fisher sampling)
    IntegerVector hap_cnt(4);
    R::rmultinom(2 * pop_siz, prob.begin(), 4, hap_cnt.begin());
    hap_frq = as<arma::dcolvec>(hap_cnt) / 2 / pop_siz;
    
    hap_frq_pth.col(t) = hap_frq;
  }
  
  return hap_frq_pth;
}
/*************************/


/********** WFD **********/
// Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
// [[Rcpp::export]]
arma::drowvec simulateOLWFDS_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const double& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // rescale the selection coefficient
  double scl_sel_cof = 2 * pop_siz * sel_cof;
  
  // declare delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  // declare delta W
  arma::drowvec dW = pow(dt, 0.5) * arma::randn<arma::drowvec>(arma::uword(lst_gen - int_gen) * ptn_num);
  
  // declare the mutant allele frequency trajectory
  arma::drowvec ale_frq_pth(arma::uword(lst_gen - int_gen) * ptn_num + 1);
  
  // initialise the mutant allele frequency in generation 0
  ale_frq_pth(0) = int_frq;
  
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) * ptn_num + 1; t++) {
    // calculate the drift coefficient
    double mu = scl_sel_cof * ale_frq_pth(t - 1) * (1 - ale_frq_pth(t - 1)) * ((1 - dom_par) - (1 - 2 * dom_par) * ale_frq_pth(t - 1));
    
    // calculate the diffusion coefficient
    double sigma = pow(ale_frq_pth(t - 1) * (1 - ale_frq_pth(t - 1)), 0.5);
    
    // proceed the Euler-Maruyama scheme
    ale_frq_pth(t) = ale_frq_pth(t - 1) + mu * dt + sigma * dW(t - 1);
    
    // remove the noise from the numerical techniques
    if (ale_frq_pth(t) < 0) {
      ale_frq_pth(t) = 0;
    }
    if (ale_frq_pth(t) > 1) {
      ale_frq_pth(t) = 1;
    }
  }
  
  return ale_frq_pth;
}

// Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
// [[Rcpp::export]]
arma::dmat simulateTLWFDS_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
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
  arma::dmat hap_frq_pth(4, arma::uword(lst_gen - int_gen) * ptn_num + 1);
  
  // initialise the haplotype frequencies in generation 0
  hap_frq_pth.col(0) = int_frq;
  
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) * ptn_num + 1; t++) {
    // calculate the drift coefficient vector
    arma::dcolvec mu(4);
    mu(0) =  scl_sel_cof_A * hap_frq_pth(0, t - 1) * (hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) * dom_par_A + (hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_A)) + scl_sel_cof_B * hap_frq_pth(0, t - 1) * (hap_frq_pth(1, t - 1) + hap_frq_pth(3, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) * dom_par_B + (hap_frq_pth(1, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_B)) - 0.5 * scl_rec_rat * (hap_frq_pth(0, t - 1) * hap_frq_pth(3, t - 1) - hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1));
    mu(1) =  scl_sel_cof_A * hap_frq_pth(1, t - 1) * (hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) * dom_par_A + (hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_A)) - scl_sel_cof_B * hap_frq_pth(1, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) * dom_par_B + (hap_frq_pth(1, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_B)) + 0.5 * scl_rec_rat * (hap_frq_pth(0, t - 1) * hap_frq_pth(3, t - 1) - hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1));
    mu(2) = -scl_sel_cof_A * hap_frq_pth(2, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) * dom_par_A + (hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_A)) + scl_sel_cof_B * hap_frq_pth(2, t - 1) * (hap_frq_pth(1, t - 1) + hap_frq_pth(3, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) * dom_par_B + (hap_frq_pth(1, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_B)) + 0.5 * scl_rec_rat * (hap_frq_pth(0, t - 1) * hap_frq_pth(3, t - 1) - hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1));
    mu(3) = -scl_sel_cof_A * hap_frq_pth(3, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) * dom_par_A + (hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_A)) - scl_sel_cof_B * hap_frq_pth(3, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) * dom_par_B + (hap_frq_pth(1, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_B)) - 0.5 * scl_rec_rat * (hap_frq_pth(0, t - 1) * hap_frq_pth(3, t - 1) - hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1));
    
    // calculate the diffusion coefficient matrix
    arma::dmat sigma(4, 6);
    sigma(0, 0) = pow(hap_frq_pth(0, t - 1) * hap_frq_pth(1, t - 1), 0.5);
    sigma(0, 1) = pow(hap_frq_pth(0, t - 1) * hap_frq_pth(2, t - 1), 0.5);
    sigma(0, 2) = pow(hap_frq_pth(0, t - 1) * hap_frq_pth(3, t - 1), 0.5);
    sigma(0, 3) = 0;
    sigma(0, 4) = 0;
    sigma(0, 5) = 0;
    sigma(1, 0) = -pow(hap_frq_pth(1, t - 1) * hap_frq_pth(0, t - 1), 0.5);
    sigma(1, 1) = 0;
    sigma(1, 2) = 0;
    sigma(1, 3) = pow(hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1), 0.5);
    sigma(1, 4) = pow(hap_frq_pth(1, t - 1) * hap_frq_pth(3, t - 1), 0.5);
    sigma(1, 5) = 0;
    sigma(2, 0) = 0;
    sigma(2, 1) = -pow(hap_frq_pth(2, t - 1) * hap_frq_pth(0, t - 1), 0.5);
    sigma(2, 2) = 0;
    sigma(2, 3) = -pow(hap_frq_pth(2, t - 1) * hap_frq_pth(1, t - 1), 0.5);
    sigma(2, 4) = 0;
    sigma(2, 5) = pow(hap_frq_pth(2, t - 1) * hap_frq_pth(3, t - 1), 0.5);
    sigma(3, 0) = 0;
    sigma(3, 1) = 0;
    sigma(3, 2) = -pow(hap_frq_pth(3, t - 1) * hap_frq_pth(0, t - 1), 0.5);
    sigma(3, 3) = 0;
    sigma(3, 4) = -pow(hap_frq_pth(3, t - 1) * hap_frq_pth(1, t - 1), 0.5);
    sigma(3, 5) = -pow(hap_frq_pth(3, t - 1) * hap_frq_pth(2, t - 1), 0.5);
    
    // proceed the Euler-Maruyama scheme
    hap_frq_pth.col(t) = hap_frq_pth.col(t - 1) + mu * dt + sigma * dW.col(t - 1);
    
    // remove the noise from the numerical techniques
    for(arma::uword i = 0; i < 4; i++) {
      if(hap_frq_pth(i, t) < 0) {
        hap_frq_pth(i, t) = 0;
      } 
      if(hap_frq_pth(i, t) > 1) {
        hap_frq_pth(i, t) = 1;
      }
    }
    hap_frq_pth.col(t) = hap_frq_pth.col(t) / sum(hap_frq_pth.col(t));
  }
  
  return hap_frq_pth;
}
/*************************/


/*********** GP **********/
// Simulate the mutant allele frequency trajectory according to the guided process of Delyon & Hu (2006) using the Euler-Maruyama method



// Simulate the mutant allele frequency trajectory according to the guided process of Fearnhead (2008) using the Euler-Maruyama method



// Simulate the mutant allele frequency trajectory according to the guided process of He et al. (2020) using the Euler-Maruyama method



/*************************/


/********** TPD **********/
// Approximate the transition probability density of the two-locus Wright-Fisher model with selection using Monte Carlo integration



// Approximate the transition probability density of the two-locus Wright-Fisher model with selection using Monte Carlo integration with importance sampling



/*************************/
