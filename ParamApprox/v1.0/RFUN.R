#' @title Moment-based approximations for the Wright-Fisher model of population dynamics under natural selection at two linked loci
#' @author Zhangyi He, Wenyang Lyu, Mark Beaumont and Feng Yu

#' version 1.0

#' R functions

#install.packages("MASS")
library("MASS")

#install.packages("coda")
library("coda")

#install.packages("inline")
library("inline")
#install.packages("Rcpp")
library("Rcpp")
#install.packages("RcppArmadillo")
library("RcppArmadillo")

#install.packages("compiler")
library("compiler")
#enableJIT(1)

# call C++ functions
sourceCpp("./ParamApprox/v1.0/CFUN.cpp")

################################################################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories

#' Standard version
simulateWFM <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]

  fts_mat <- calculateFitnessMat_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B)

  frq_pth <- simulateWFM_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen)

  return(frq_pth)
}
#' Compiled version
cmpsimulateWFM <- cmpfun(simulateWFM)

########################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

#' Standard version
simulateDiffusApprox <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]

  frq_pth <- simulateDiffusApprox_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)

  if (dat_aug == FALSE) {
    return(frq_pth[, (0:(lst_gen - int_gen)) * ptn_num + 1])
  } else {
    return(frq_pth)
  }
}
#' Compiled version
cmpsimulateDiffusApprox <- cmpfun(simulateDiffusApprox)

########################################

#' Approximate the first two moments of the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param sim_num the number of the samples in Monte Carlo simulation

#' Standard version
approximateMoment <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]



}
#' Compiled version
cmpapproximateMoment <- cmpfun(approximateMoment)

########################################

#' Approximate the two-locus Wright-Fisher model with selection using the normal
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param sim_num the number of the samples in Monte Carlo simulation

#' Standard version
approximateWFM_Norm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]


  # return the parameters of the normal approximation and the corresponding first two moments

}
#' Compiled version
cmpapproximateWFM_Norm <- cmpfun(approximateWFM_Norm)

########################################

#' Simulate the two-locus Wright-Fisher model with selection using the logistic normal
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param sim_num the number of the samples in Monte Carlo simulation
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param smp_siz the number of the Monte Carlo samples

#' Standard version
simulateWFM_LogisticNorm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, ...) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]



}
#' Compiled version
cmpsimulateWFM_LogisticNorm <- cmpfun(simulateWFM_LogisticNorm)

########################################

#' Approximate the two-locus Wright-Fisher model with selection using the logistic normal
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param smp_siz the number of the Monte Carlo samples

#' Standard version
approximateWFM_LogisticNorm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]


  # return the parameters of the logistic normal approximation and the corresponding first two moments

}
#' Compiled version
cmpapproximateWFM_LogisticNorm <- cmpfun(approximateWFM_LogisticNorm)

########################################

#' Simulate the two-locus Wright-Fisher model with selection using the logistic normal
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param sim_num the number of the samples in Monte Carlo simulation
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param smp_siz the number of the Monte Carlo samples

#' Standard version
simulateWFM_LogisticNorm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, ...) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]



}
#' Compiled version
cmpsimulateWFM_LogisticNorm <- cmpfun(simulateWFM_LogisticNorm)

########################################

#' Approximate the two-locus Wright-Fisher model with selection using the hierarchical beta
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param smp_siz the number of the Monte Carlo samples

#' Standard version
approximateWFM_HierarchicalBeta <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]


  # return the parameters of the hierarchical beta approximation and the corresponding first two moments

}
#' Compiled version
cmpapproximateWFM_HierarchicalBeta <- cmpfun(approximateWFM_HierarchicalBeta)

########################################

#' Simulate the two-locus Wright-Fisher model with selection using the hierarchical beta
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param sim_num the number of the samples in Monte Carlo simulation
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param smp_siz the number of the Monte Carlo samples

#' Standard version
simulateWFM_HierarchicalBeta <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, ...) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]



}
#' Compiled version
cmpsimulateWFM_HierarchicalBeta <- cmpfun(simulateWFM_HierarchicalBeta)

################################################################################
