#' @title Moment-based approximations for the Wright-Fisher model of population dynamics under natural selection at two linked loci
#' @author Zhangyi He, Wenyang Lyu, Mark Beaumont and Feng Yu

#' version 1.0

#' R functions

#install.packages("mltools")
library("mltools")
#install.packages("data.table")
library("data.table")

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
sourceCpp("./Code/Code v1.0/CFUN.cpp")

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
  fts_mat <- calculateFitnessMat_arma(sel_cof[1], dom_par[1], sel_cof[2], dom_par[2])
  frq_pth <- simulateWFM_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
  frq_pth <- as.matrix(frq_pth)

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
  frq_pth <- simulateDiffusApprox_arma(sel_cof[1], dom_par[1], sel_cof[2], dom_par[2], rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num)
  frq_pth <- as.matrix(frq_pth)

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
#' @param smp_siz the number of the Monte Carlo samples for the first two moment approximations

#' Standard version
approximateMoment_WFM <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  fts_mat <- calculateFitnessMat_arma(sel_cof[1], dom_par[1], sel_cof[2], dom_par[2])

  # Approximate the first two moments of the Wright-Fisher model using
  if (mnt_apx == "MC") {
    # use Monte Carlo simulations
    mnt <- approximateMoment_MonteCarlo_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen, smp_siz)
  }
  if (mnt_apx == "Lacerda") {
    # use the extension of Lacerda & Seoighe (2014)
    mnt <- approximateMoment_Lacerda_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
  }
  if (mnt_apx == "Terhorst") {
    # use the extension of Terhorst et al. (2015)
    mnt <- approximateMoment_Terhorst_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
  }
  if (mnt_apx == "Paris1") {
    # use the extension of Paris et al. (2019) with the first-order Taylor expansion
    mnt <- approximateMoment_Paris1_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
  }
  if (mnt_apx == "Paris2") {
    # use the extension of Paris et al. (2019) with the second-order Taylor expansion
    mnt <- approximateMoment_Paris2_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
  }

  return(list(mean = as.matrix(mnt$mean),
              variance = as.array(mnt$variance)))
}
#' Compiled version
cmpapproximateMoment_WFM <- cmpfun(approximateMoment_WFM)

########################################

#' Generate the samples under the normal approximation of the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param sim_num the number of the Monte Carlo samples
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param smp_siz the number of the Monte Carlo samples for the first two moment approximations

#' Standard version
generateSample_Norm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    param <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, smp_siz)
  } else {
    param <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  frq_pth <- simulateWFM_norm_arma(param$mean, param$variance, int_gen, lst_gen, sim_num)
  frq_pth <- as.array(frq_pth)

  return(frq_pth)
}
#' Compiled version
cmpgenerateSample_Norm <- cmpfun(generateSample_Norm)

########################################

#' Calculate the parameters for the logistic normal approximation of the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param smp_siz the number of the Monte Carlo samples for the first two moment approximations

#' Standard version
calculateParam_LogisticNorm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, smp_siz)
  } else {
    mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  param = calculateParam_LogisticNorm_arma(mnt$mean, mnt$variance, int_gen, lst_gen)

  return(list(location = as.matrix(mnt$location),
              scalesq = as.array(mnt$scalesq)))
}
#' Compiled version
cmpcalculateParam_LogisticNorm <- cmpfun(calculateParam_LogisticNorm)

########################################

#' Generate the samples under the logistic normal approximation of the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param sim_num the number of the Monte Carlo samples
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param smp_siz the number of the Monte Carlo samples for the first two moment approximations

#' Standard version
generateSample_LogisticNorm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    param <- cmpcalculateParam_LogisticNorm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, smp_siz) 
  } else {
    param <- cmpcalculateParam_LogisticNorm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx) 
  }
  frq_pth = generateSample_LogisticNorm_arma(param$location, param$scalesq, int_gen, lst_gen, sim_num)
  frq_pth <- as.array(frq_pth)

  return(frq_pth)
}
#' Compiled version
cmpgenerateSample_LogisticNorm <- cmpfun(generateSample_LogisticNorm)

########################################

#' Approximate the moments of the logistic normal approximation of the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param sim_num the number of the Monte Carlo samples
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param smp_siz the number of the Monte Carlo samples for the first two moment approximations

#' Standard version
approximateMoment_LogisticNorm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    param <- cmpcalculateParam_LogisticNorm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, smp_siz) 
  } else {
    param <- cmpcalculateParam_LogisticNorm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx) 
  }
  mnt = approximateMoment_LogisticNorm_arma(param$location, param$scalesq, int_gen, lst_gen, sim_num)

  return(list(mean = as.matrix(mnt$mean),
              variance = as.array(mnt$variance)))
}
#' Compiled version
cmpapproximateMoment_LogisticNorm <- cmpfun(approximateMoment_LogisticNorm)

########################################

#' Calculate the parameters for the hierarchical beta approximation of the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param smp_siz the number of the Monte Carlo samples for the first two moment approximations

#' Standard version
calculateParam_HierarchicalBeta <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, smp_siz)
  } else {
    mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  param = calculateParam_HierarchicalBeta_arma(mnt$mean, mnt$variance, int_gen, lst_gen)

  return(list(alpha = as.matrix(param$alpha),
              beta = as.matrix(param$beta)))
}
#' Compiled version
cmpcalculateParam_HierarchicalBeta <- cmpfun(calculateParam_HierarchicalBeta)

########################################

#' Generate the samples under the hierarchical beta approximation of the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param sim_num the number of the Monte Carlo samples
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param smp_siz the number of the Monte Carlo samples for the first two moment approximations

#' Standard version
generateSample_HierarchicalBeta <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    param <- cmpapproximateWFM_HierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, smp_siz)
  } else {
    param <- cmpapproximateWFM_HierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  frq_pth = generateSample_HierarchicalBeta_arma(param$alpha, param$beta, int_gen, lst_gen, sim_num)
  frq_pth <- as.array(frq_pth)

  return(frq_pth)
}
#' Compiled version
cmpgenerateSample_HierarchicalBeta <- cmpfun(generateSample_HierarchicalBeta)

########################################

#' Approximate the moments of the hierarchical beta approximation of the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param smp_siz the number of the Monte Carlo samples for the first two moment approximations

#' Standard version
approximateMoment_HierarchicalBeta <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    param <- cmpapproximateWFM_HierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, smp_siz)
  } else {
    param <- cmpapproximateWFM_HierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  mnt = approximateMoment_HierarchicalBeta_arma(param$alpha, param$beta, int_gen, lst_gen)

  return(list(mean = as.matrix(mnt$mean),
              variance = as.array(mnt$variance)))
}
#' Compiled version
cmpapproximateMoment_HierarchicalBeta <- cmpfun(approximateMoment_HierarchicalBeta)

########################################





#' Calculate the root mean square deviation between the empirical cumulative distribution functions for the Wright-Fisher model/diffusion
#' Parameter setting
#' @param smp_WFM the sample trajectories generated under the Wright-Fisher model
#' @param smp_WFD the sample trajectories generated under the Wright-Fisher diffusion
#' @param grd_num the grid number for the empirical probability distribution function
#' @param rnd_grd = TRUE/FALSE (return the root mean square deviation between the empirical cumulative distribution functions with random grid points or not)

#' Standard version
calculateRMSD <- function(smp_WFM, smp_WFD, grd_num, rnd_grd = TRUE) {
  if (rnd_grd == TRUE) {
    #frq_grd <- generateRndGrid_2L_arma(grd_num)
    #frq_grd <- t(frq_grd)
    frq_grd <- generateRndGridUnif_2L_arma(grd_num)
  } else {
    frq_grd <- generateFixGrid_2L_arma(grd_num)

  }
  frq_grd <- data.table(frq_grd)
  setnames(frq_grd, c("A1B1", "A1B2", "A2B1", "A2B2"))

  dist <- rep(NA, length.out = dim(smp_WFM)[3])
  for (k in 1:dim(smp_WFM)[3]) {
    cdf_WFM <- empirical_cdf(data.table(A1B1 = smp_WFM[, 1, k], A1B2 = smp_WFM[, 2, k], A2B1 = smp_WFM[, 3, k], A2B2 = smp_WFM[, 4, k]), ubounds = frq_grd)
    #cdf_WFM <- empirical_cdf(data.table(A1B1 = smp_WFM[, 1, k], A1B2 = smp_WFM[, 2, k], A2B1 = smp_WFM[, 3, k], A2B2 = smp_WFM[, 4, k]), ubounds = CJ(A1B1=1:3/3, A1B2=1:3/3, A2B1=1:3/3, A2B2=1:3/3))
    print(cdf_WFM)
    cdf_WFD <- empirical_cdf(data.table(A1B1 = smp_WFD[, 1, k], A1B2 = smp_WFD[, 2, k], A2B1 = smp_WFD[, 3, k], A2B2 = smp_WFD[, 4, k]), ubounds = frq_grd)
    dist[k] <- sqrt(sum((cdf_WFM - cdf_WFD)^2) / dim(smp_WFM)[1])
  }

  return(dist)
}
#' Compiled version
cmpcalculateRMSD <- cmpfun(calculateRMSD)

#' Standard version
calculateRMSD_2nd <- function(smp_WFM, smp_WFD, int_gen, lst_gen, sim_num, grd_num, rnd_grd = TRUE) {
  if (rnd_grd == TRUE) {
    frq_grd <- generateRndGrid_2L_arma(grd_num)
    #frq_grd <- t(frq_grd)
    #frq_grd <- generateRndGridUnif_2L_arma(grd_num)
  } else {
    frq_grd <- generateFixGrid_2L_arma(grd_num)

  }

  dist <- rep(NA, length.out = dim(smp_WFM)[3])
  cdf_WFM = calculate_empirical_cdf_2L_arma(smp_WFM, int_gen, lst_gen, sim_num, frq_grd)
  cdf_WFD = calculate_empirical_cdf_2L_arma(smp_WFD, int_gen, lst_gen, sim_num, frq_grd)
  for (k in 1:dim(smp_WFM)[3]) {
    dist[k] <- sqrt(sum((cdf_WFM[, k] - cdf_WFD[, k])^2) / dim(smp_WFM)[1])
  }

  return(dist)
}
#' Compiled version
cmpcalculateRMSD_2nd <- cmpfun(calculateRMSD_2nd)

################################################################################
