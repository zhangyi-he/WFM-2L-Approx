#' @title Numerical simulation of the two-locus Wright-Fisher diffusion with application to approximating transition probability densities
#' @author Zhangyi He, Mark Beaumont and Feng Yu

#' version 1.0

#' R functions

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
sourceCpp("./DiffusApprox/v1.0/CFUN.cpp")

################################################################################

#' Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial mutant allele frequency of the population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory

#' Standard version
simulateWFM_1L <- function(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen) {
  frq_pth <- simulateWFM_1L_arma(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen)
  frq_pth <- as.vector(frq_pth)

  return(frq_pth)
}
#' Compiled version
cmpsimulateWFM_1L <- cmpfun(simulateWFM_1L)

########################################

#' Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial mutant allele frequency of the population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

#' Standard version
simulateDiffusApprox_1L <- function(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE) {
  frq_pth <- simulateDiffusApprox_1L_arma(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num)
  frq_pth <- as.vector(frq_pth)

  if (dat_aug == FALSE) {
    return(frq_pth[(0:(lst_gen - int_gen)) * ptn_num + 1])
  } else {
    return(frq_pth)
  }
}
#' Compiled version
cmpsimulateDiffusApprox_1L <- cmpfun(simulateDiffusApprox_1L)

########################################

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
simulateWFM_2L <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]

  fts_mat <- calculateFitnessMat_2L_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B)

  frq_pth <- simulateWFM_2L_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen)

  return(frq_pth)
}
#' Compiled version
cmpsimulateWFM_2L <- cmpfun(simulateWFM_2L)

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
simulateDiffusApprox_2L <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]

  frq_pth <- simulateDiffusApprox_2L_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)

  if (dat_aug == FALSE) {
    return(frq_pth[, (0:(lst_gen - int_gen)) * ptn_num + 1])
  } else {
    return(frq_pth)
  }
}
#' Compiled version
cmpsimulateDiffusApprox_2L <- cmpfun(simulateDiffusApprox_2L)

########################################

#' Compute the empirical cumulative distribution function for the Wright-Fisher model or the Wright-Fisher diffusion
#' Parameter setting
#' @param model = WFM/WFD (return the empirical cumulative distribution function for the Wright-Fisher model or the Wright-Fisher diffusion)
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param sim_num the number of the samples in Monte Carlo simulation
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method



########################################

#' Calculate the distance between the empirical cumulative distribution functions for the Wright-Fisher model or the Wright-Fisher diffusion
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param sim_num the number of the samples in Monte Carlo simulation
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method

sel_cof_A <- sel_cof[1]
sel_cof_B <- sel_cof[2]
dom_par_A <- dom_par[1]
dom_par_B <- dom_par[2]

fts_mat <- calculateFitnessMat_2L_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B)
WFM_smp <- generateSample_WFM_2L_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num)

WFD_smp <- generateSample_WFD_2L_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num)
WFD_smp <- WFD_smp[, , (0:(lst_gen - int_gen)) * ptn_num + 1]

for (k in 1:(lst_gen - int_gen + 1)) {
  statistic[k] <- npdeneqtest(WFM_smp[, , k], WFD_smp[, , k], boot.num = 1e+02)
}


########################################

#' Simulate the haplotype frequency trajectories according to the guided process using the Euler-Maruyama method
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
#' @param guided_proc the guided process (Fearnhead (2008) or He et al. (2020))
#' @param modification = TRUE/FALSE
#' @param rho



########################################

#' Approximate the transition probability density of the two-locus Wright-Fisher model with selection using Monte Carlo integration
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param sim_num the number of Monte Carlo simulations



########################################

#' Calculate the importance weight
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param guided_proc the guided process (Fearnhead (2008) or He et al. (2020))
#' @param modification = TRUE/FALSE
#' @param rho



########################################

#' Approximate the transition probability density of the two-locus Wright-Fisher model with selection using Monte Carlo integration with importance sampling
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param sim_num the number of Monte Carlo simulations
#' @param guided_proc the guided process (Fearnhead (2008) or He et al. (2020))
#' @param modification = TRUE/FALSE
#' @param rho



################################################################################
