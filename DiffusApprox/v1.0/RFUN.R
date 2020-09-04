#' @title Numerical simulation of the two-locus Wright-Fisher diffusion with application to approximating transition probability densities
#' @author Zhangyi He, Mark Beaumont and Feng Yu

#' version 1.0

#' R functions

# install.packages("np")
library("np")

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
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

#' Standard version
simulateWFD_1L <- function(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE) {
  frq_pth <- simulateWFD_1L_arma(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num)
  frq_pth <- as.vector(frq_pth)

  if (dat_aug == FALSE) {
    return(frq_pth[(0:(lst_gen - int_gen)) * ptn_num + 1])
  } else {
    return(frq_pth)
  }
}
#' Compiled version
cmpsimulateWFD_1L <- cmpfun(simulateWFD_1L)

########################################

#' Generate the sample mutant allele frequency trajectories according to the one-locus Wright-Fisher model/diffusion with selection
#' Parameter setting
#' @param model = WFM/WFD (return the sample mutant allele frequency trajectories according to the one-locus Wright-Fisher model/diffusion with selection)
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial mutant allele frequency of the population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory
#' @param sim_num the number of the samples in Monte Carlo simulation
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method

#' Standard version
generateSampleTraj_1L <- function(model, sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, sim_num, ...) {
  if (model == "WFM") {
    smp_pth <- generateSample_WFM_1L_arma(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, sim_num)
  } else {
    smp_pth <- generateSample_WFD_1L_arma(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num)
    smp_pth <- smp_pth[, (0:(lst_gen - int_gen)) * ptn_num + 1]
  }

  return(smp_pth)
}
#' Compiled version
cmpgenerateSampleTraj_1L <- cmpfun(generateSampleTraj_1L)

########################################

#' Calculate the Hellinger distance between the empirical probability distribution functions for the one-locus Wright-Fisher model/diffusion with selection
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial mutant allele frequency of the population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory
#' @param sim_num the number of the samples in Monte Carlo simulation
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method

#' Standard version
calculateHellingerDist_1L <- function(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, sim_num, ptn_num) {
  model <- "WFM"
  smp_WFM <- generateSampleTraj_1L(model, sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, sim_num)

  model <- "WFD"
  smp_WFD <- generateSampleTraj_1L(model, sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num)
  smp_WFD <- smp_WFD[(0:(lst_gen - int_gen)) * ptn_num + 1]

  dist <- rep(NA, length.out = lst_gen - int_gen + 1)
  for (k in 1:(lst_gen - int_gen + 1)) {
    pdf_WFM <- hist(smp_WFM[, k], breaks = (0:(2 * pop_siz)) / (2 * pop_siz), plot = FALSE)$counts / sim_num
    pdf_WFD <- hist(smp_WFD[, k], breaks = (0:(2 * pop_siz)) / (2 * pop_siz), plot = FALSE)$counts / sim_num
    dist[k] <- sqrt(sum((sqrt(pdf_WFM) - sqrt(pdf_WFD))^2) / 2)
  }

  return(dist)
}
#' Compiled version
cmpcalculateHellingerDist_1L <- cmpfun(calculateHellingerDist_1L)

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
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

#' Standard version
simulateWFD_2L <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]

  frq_pth <- simulateWFD_2L_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)

  if (dat_aug == FALSE) {
    return(frq_pth[, (0:(lst_gen - int_gen)) * ptn_num + 1])
  } else {
    return(frq_pth)
  }
}
#' Compiled version
cmpsimulateWFD_2L <- cmpfun(simulateWFD_2L)

########################################

#' Generate the sample haplotype frequency trajectories according to the two-locus Wright-Fisher model/diffusion with selection
#' Parameter setting
#' @param model = WFM/WFD (return the sample haplotype frequency trajectories according to the two-locus Wright-Fisher model/diffusion with selection)
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param sim_num the number of the samples in Monte Carlo simulation
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method

#' Standard version
generateSampleTraj_2L <- function(model, sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, ...) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]

  if (model == "WFM") {
    fts_mat <- calculateFitnessMat_2L_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B)
    smp_pth <- generateSample_WFM_2L_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num)
  } else {
    smp_pth <- generateSample_WFD_2L_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num)
    smp_pth <- smp_pth[, , (0:(lst_gen - int_gen)) * ptn_num + 1]
  }

  return(smp_pth)
}
#' Compiled version
cmpgenerateSampleTraj_2L <- cmpfun(generateSampleTraj_2L)

########################################

#' Calculate the Hellinger distance between the empirical (marginal) probability distribution functions for the two-locus Wright-Fisher model/diffusion with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param sim_num the number of the samples in Monte Carlo simulation
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method

#' Standard version
calculateHellingerDist_2L <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, ptn_num) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]

  model <- "WFM"
  smp_WFM <- generateSampleTraj_2L(model, sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num)

  model <- "WFD"
  smp_WFD <- generateSampleTraj_2L(model, sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num)
  smp_WFD <- smp_WFD[, , (0:(lst_gen - int_gen)) * ptn_num + 1]

  dist <- matrix(NA, nrow = 4, ncol = lst_gen - int_gen + 1)
  for (k in 1:(lst_gen - int_gen + 1)) {
    for (i in 1:4) {
      pdf_WFM <- hist(smp_WFM[, i, k], breaks = (0:(2 * pop_siz)) / (2 * pop_siz), plot = FALSE)$counts / sim_num
      pdf_WFD <- hist(smp_WFD[, i, k], breaks = (0:(2 * pop_siz)) / (2 * pop_siz), plot = FALSE)$counts / sim_num
      dist[i, k] <- sqrt(sum((sqrt(pdf_WFM) - sqrt(pdf_WFD))^2) / 2)
    }
  }

  return(dist)
}
#' Compiled version
cmpcalculateHellingerDist_2L <- cmpfun(calculateHellingerDist_2L)

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
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method
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
