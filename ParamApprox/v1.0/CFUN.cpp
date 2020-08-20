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
// Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::dmat simulateWFM_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
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


/***** MomentApprox ******/
// Approximate the first two moments of the two-locus Wright-Fisher model with selection using Monte Carlo simulation 

// Calculate the fitness matrix for the Wright-Fisher model
// [[Rcpp::export]]
arma::dmat calculateFitnessMat_arma(const double& dom_par_A, const double& dom_par_B, const double& sel_cof_A, const double& sel_cof_B) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // calculate the fitness
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
  
  return fts;
}

arma::dcolvec calculate_p(const double& rec_rat, const arma::dmat fts, const arma::dcolvec& hap_frq){
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // declare eta
  arma::dcolvec eta(4);
  eta(0) = -1;
  eta(1) = 1;
  eta(2) = 1;
  eta(3) = -1;
  arma::dmat gen_frq = hap_frq * hap_frq.t();
  
  // viability selection
  gen_frq = (fts % gen_frq) / arma::as_scalar(sum(sum(fts % gen_frq, 0), 1));
  
  // meiosis (calculate the sampling probability)
  arma::dcolvec prob = arma::zeros<arma::dcolvec>(4);
  for(arma::uword i = 0; i < 4; i++) {
    prob(i) = (sum(gen_frq.row(i)) + sum(gen_frq.col(i))) / 2;
  }
  prob = prob + eta * rec_rat * (prob(0) * prob(3) - prob(1) * prob(2));
  
  return prob;
}

arma::dmat calculate_grad_mu(const arma::dmat fts, const arma::dcolvec& hap_frq){
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::dmat grad_mu(4, 4);
  arma::dmat gen_frq = hap_frq * hap_frq.t();
  
  // viability selection
  arma::dcolvec q = sum(fts % gen_frq, 1) / arma::as_scalar(sum(sum(fts % gen_frq, 0), 1));
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
  grad_mu = (iden_mat + Q) * (arma::diagmat(hap_frq) * (fts / arma::as_scalar(sum(sum(fts % gen_frq, 0), 1)) - 2 * q / hap_frq * (q / hap_frq).t() ) + arma::diagmat(q / hap_frq) );
  return grad_mu;
}

List approximateMoment_MonteCarlo_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  arma::dmat mu(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma(16, 16, arma::uword(lst_gen - int_gen) + 1);
  mu.col(0) = int_frq;
  sigma.slice(0) = arma::zeros<arma::dmat>(4, 4);;
  int num_sim = 1e1; // number of simulations per generation
  arma::dcube path(4, arma::uword(lst_gen - int_gen) + 1, num_sim);
  for(arma::uword i = 1; i < num_sim; i++) {
    path.slice(i) = simulateWFM_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen);
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    arma::dmat pth_k = path.col(k);
    mu.col(k) = arma::mean(pth_k, 0);
    sigma.slice(k) = arma::cov(pth_k.t(), pth_k.t());
    }

  // return the approximations for the mean and variance of the Wright-Fisher model at each generation from int_gen to lst_gen
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}

// Approximate the first two moments of the two-locus Wright-Fisher model with selection using the extension of Lacerda & Seoighe (2014) 
List approximateMoment_Lacerda_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::dmat fts) {
  // ensure RNG gets set/reset
  RNGScope scope;
  arma::dmat mu(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma(16, 16, arma::uword(lst_gen - int_gen) + 1);
  mu.col(0) = int_frq;
  sigma.slice(0) = arma::zeros<arma::dmat>(4, 4);
  int num_sim = 1e1; // number of simulations per generation
  
  for(arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    mu.col(t) = calculate_p(rec_rat, fts, mu.col(t-1));
    arma::dmat grad_mu = calculate_grad_mu(fts, mu.col(t-1));
    sigma.slice(t) = (arma::diagmat(mu.col(t)) - mu.col(t) * (mu.col(t)).t() ) / 2 / pop_siz + grad_mu * sigma.slice(t-1) * grad_mu.t();
    }

  // return the approximations for the mean and variance of the Wright-Fisher model at each generation from int_gen to lst_gen
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}

// Approximate the first two moments of the two-locus Wright-Fisher model with selection using the extension of Terhorst et al. (2015) 
List approximateMoment_Terhorst_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::dmat fts) {
  // ensure RNG gets set/reset
  RNGScope scope;
  arma::dmat x(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dmat mu(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dmat epsilon(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma(16, 16, arma::uword(lst_gen - int_gen) + 1);
  mu.col(0) = int_frq;
  x.col(0) = int_frq;
  sigma.slice(0) = arma::zeros<arma::dmat>(4, 4);;
  int num_sim = 1e1; // number of simulations per generation
  
  for(arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    arma::dmat path(4, num_sim);
    for(arma::uword i = 1; i < num_sim; i++) {
      path.col(i) = simulateWFM_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, x.col(t-1), t-1, t);
    x.col(t) = arma::mean(path, 0);
    arma::dmat grad_mu = calculate_grad_mu(fts, x.col(t-1));
    epsilon.col(t) = grad_mu * epsilon.col(t-1);
    mu.col(t) = x.col(t) + epsilon.col(t);
    sigma.slice(t) = 1 / 2 / pop_siz * mu.col(t) * arma::diagmat(mu.col(t)) - 1 / 2 / pop_siz * mu.col(t) * (mu.col(t)).t() + (1 - 1 / 2 / pop_siz) * grad_mu * sigma.slice(t-1) * grad_mu.t();
    }

  // return the approximations for the mean and variance of the Wright-Fisher model at each generation from int_gen to lst_gen
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}

// Approximate the first two moments of the two-locus Wright-Fisher model with selection using the extension of Paris et al. (2019) 
List approximateMoment_Paris_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::dmat fts) {
  // ensure RNG gets set/reset
  RNGScope scope;
  arma::dmat mu(4, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube sigma(16, 16, arma::uword(lst_gen - int_gen) + 1);
  mu.col(0) = int_frq;
  sigma.slice(0) = arma::zeros<arma::dmat>(4, 4);;
  int num_sim = 1e1; // number of simulations per generation
  
  for(arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    mu.col(t) = calculate_p(rec_rat, fts, mu.col(t-1));
    arma::dmat grad_mu = calculate_grad_mu(fts, mu.col(t-1));
    sigma.slice(t) = (arma::diagmat(mu.col(t)) - mu.col(t) * (mu.col(t)).t() ) / 2 / pop_siz + (1 - 1 / 2 / pop_siz) * grad_mu * sigma.slice(t-1) * grad_mu.t();
    }

  // return the approximations for the mean and variance of the Wright-Fisher model at each generation from int_gen to lst_gen
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}

/*************************/


/****** ParamApprox ******/
// Approximate the two-locus Wright-Fisher model with selection using the normal distribution
List approximatWFM_norm_arma(const dmat& mean, const dcube& variance) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // return the parameters of the normal distribution
  return List::create(Named("mean", mu),
                      Named("variance", sigma));
}

  arma::dcolvec calculate_phi(const arma::dcolvec& hap_frq){
    // ensure RNG gets set/reset
    RNGScope scope;

    arma::dcolvec phi(3);
    phi(0) = arma::log(hap_frq(0) / hap_frq(3));
    phi(1) = arma::log(hap_frq(1) / hap_frq(3));
    phi(2) = arma::log(hap_frq(2) / hap_frq(3));
    return phi;
  }
  
  arma::dmat calculate_grad_phi(const arma::dcolvec& hap_frq){
    // ensure RNG gets set/reset
    RNGScope scope;
    
    arma::dmat grad_phi(3, 4);
    grad_phi(0, 0) = 1 / hap_frq(0);
    grad_phi(1, 1) = 1 / hap_frq(1);
    grad_phi(2, 2) = 1 / hap_frq(2);
    grad_phi(0, 3) = -1 / hap_frq(3);
    grad_phi(1, 3) = -1 / hap_frq(3);
    grad_phi(2, 3) = -1 / hap_frq(3);
    
    return grad_phi;
  }
  
// Approximate the two-locus Wright-Fisher model with selection using the logistic normal distribution
List approximatWFM_logisticnorm_arma(const dmat& mean, const dcube& variance, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::dmat m(3, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube V(9, 9, arma::uword(lst_gen - int_gen) + 1);
  for(arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    m.col(t) = calculate_phi(mean.col(t));
    V.slice(t) = calculate_grad_phi(mean.col(t)) * variance.slice(t) * calculate_grad_phi(mean.col(t)).t();
  }

  // return the parameters of the logit-normal distribution
  return List::create(Named("mean", m),
                      Named("variance", V));
}

// Approximate the first two moments of the logistic normal approximation of the two-locus Wright-Fisher model with selection
arma::dmat approximateMoment_logisticnorm_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;


  // return the approximations for the mean and variance of the Wright-Fisher model at each generation from int_gen to lst_gen
}

// Approximate the two-locus Wright-Fisher model with selection using the hierarchical beta distribution
List approximatWFM_hierarchicalbeta_arma(const dmat& mean, const dcube& variance) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::dmat m(3, arma::uword(lst_gen - int_gen) + 1);
  arma::dcube V(9, 9, arma::uword(lst_gen - int_gen) + 1);
  for(arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    m.col(t) = ??;
    V.slice(t) = ??;
  }


  // return the parameters of the hierarchical beta distribution
  return List::create(Named("mean", m),
                      Named("variance", V));
}

// Approximate the first two moments of the hierarchical beta approximation of the two-locus Wright-Fisher model with selection
arma::dmat approximateMoment_hierarchicalbeta_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& sim_num) {
  // ensure RNG gets set/reset
  RNGScope scope;


  // return the approximations for the mean and variance of the Wright-Fisher model at each generation from int_gen to lst_gen
}

/*************************/
