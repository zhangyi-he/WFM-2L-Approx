#' @title Moment-based approximations for the Wright-Fisher model of population dynamics under natural selection at two linked loci
#' @author Zhangyi He, Wenyang Lyu, Mark Beaumont and Feng Yu

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
#' @param mnt_num the number of the Monte Carlo samples for the moment approximations

#' Standard version
approximateMoment_WFM <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  fts_mat <- calculateFitnessMat_arma(sel_cof[1], dom_par[1], sel_cof[2], dom_par[2])
  if (mnt_apx == "MC") {
    # use Monte Carlo simulations
    mnt <- approximateMoment_MonteCarlo_arma(fts_mat, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_num)
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

#' Calculate the parameters for the normal approximation of the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param mnt_num the number of the Monte Carlo samples for the moment approximations

#' Standard version
calculateParam_Norm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num)
  } else {
    mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  param <- calculateParam_Norm_arma(mnt$mean, mnt$variance, int_gen, lst_gen)

  return(list(mean = as.matrix(param$mean),
              variance = as.array(param$variance)))
}
#' Compiled version
cmpcalculateParam_Norm <- cmpfun(calculateParam_Norm)

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
#' @param mnt_num the number of the Monte Carlo samples for the moment approximations

#' Standard version
generateSample_Norm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    param <- cmpcalculateParam_Norm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num)
  } else {
    param <- cmpcalculateParam_Norm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  frq_pth <- generateSample_Norm_arma(param$mean, param$variance, int_gen, lst_gen, sim_num)
  frq_pth <- as.array(frq_pth)

  return(frq_pth)
}
#' Compiled version
cmpgenerateSample_Norm <- cmpfun(generateSample_Norm)

########################################

#' Approximate the moments of the normal approximation of the two-locus Wright-Fisher model with selection
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
#' @param mnt_num the number of the Monte Carlo samples for the moment approximations

#' Standard version
approximateMoment_Norm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    param <- cmpcalculateParam_Norm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num)
  } else {
    param <- cmpcalculateParam_Norm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  mnt <- approximateMoment_Norm_arma(param$mean, param$variance, int_gen, lst_gen, sim_num)

  return(list(mean = as.matrix(mnt$mean),
              variance = as.array(mnt$variance)))
}
#' Compiled version
cmpapproximateMoment_Norm <- cmpfun(approximateMoment_Norm)

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
#' @param mnt_num the number of the Monte Carlo samples for the moment approximations

#' Standard version
calculateParam_LogisticNorm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num)
  } else {
    mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  param <- calculateParam_LogisticNorm_arma(mnt$mean, mnt$variance, int_gen, lst_gen)

  return(list(location = as.matrix(param$location),
              scalesq = as.array(param$scalesq)))
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
#' @param mnt_num the number of the Monte Carlo samples for the moment approximations

#' Standard version
generateSample_LogisticNorm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    param <- cmpcalculateParam_LogisticNorm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num)
  } else {
    param <- cmpcalculateParam_LogisticNorm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  frq_pth <- generateSample_LogisticNorm_arma(param$location, param$scalesq, int_gen, lst_gen, sim_num)
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
#' @param mnt_num the number of the Monte Carlo samples for the moment approximations

#' Standard version
approximateMoment_LogisticNorm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    param <- cmpcalculateParam_LogisticNorm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num)
  } else {
    param <- cmpcalculateParam_LogisticNorm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  mnt <- approximateMoment_LogisticNorm_arma(param$location, param$scalesq, int_gen, lst_gen, sim_num)

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
#' @param mnt_num the number of the Monte Carlo samples for the moment approximations

#' Standard version
calculateParam_HierarchicalBeta <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num)
  } else {
    mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  param <- calculateParam_HierarchicalBeta_arma(mnt$mean, mnt$variance, int_gen, lst_gen)

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
#' @param mnt_num the number of the Monte Carlo samples for the moment approximations

#' Standard version
generateSample_HierarchicalBeta <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    param <- cmpcalculateParam_HierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num)
  } else {
    param <- cmpcalculateParam_HierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  frq_pth <- generateSample_HierarchicalBeta_arma(param$alpha, param$beta, int_gen, lst_gen, sim_num)
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
#' @param mnt_num the number of the Monte Carlo samples for the moment approximations

#' Standard version
approximateMoment_HierarchicalBeta <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    param <- cmpcalculateParam_HierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num)
  } else {
    param <- cmpcalculateParam_HierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  mnt <- approximateMoment_HierarchicalBeta_arma(param$alpha, param$beta, int_gen, lst_gen)

  return(list(mean = as.matrix(mnt$mean),
              variance = as.array(mnt$variance)))
}
#' Compiled version
cmpapproximateMoment_HierarchicalBeta <- cmpfun(approximateMoment_HierarchicalBeta)

########################################

#' Calculate the parameters for the pyramidal hierarchical beta approximation of the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param mnt_num the number of the Monte Carlo samples for the moment approximations

#' Standard version
calculateParam_PyramidalHierarchicalBeta <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num)
  } else {
    mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  param <- calculateParam_PyramidalHierarchicalBeta_arma(mnt$mean, mnt$variance, int_gen, lst_gen)

  return(list(alpha = as.matrix(param$alpha),
              beta = as.matrix(param$beta)))
}
#' Compiled version
cmpcalculateParam_PyramidalHierarchicalBeta <- cmpfun(calculateParam_PyramidalHierarchicalBeta)

########################################

#' Generate the samples under the pyramidal hierarchical beta approximation of the two-locus Wright-Fisher model with selection
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
#' @param mnt_num the number of the Monte Carlo samples for the moment approximations

#' Standard version
generateSample_PyramidalHierarchicalBeta <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    param <- cmpcalculateParam_PyramidalHierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num)
  } else {
    param <- cmpcalculateParam_PyramidalHierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  frq_pth <- generateSample_PyramidalHierarchicalBeta_arma(param$alpha, param$beta, int_gen, lst_gen, sim_num)
  frq_pth <- as.array(frq_pth)

  return(frq_pth)
}
#' Compiled version
cmpgenerateSample_PyramidalHierarchicalBeta <- cmpfun(generateSample_PyramidalHierarchicalBeta)

########################################

#' Approximate the moments of the pyramidal hierarchical beta approximation of the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param mnt_apx the moment approximation (Monte Carlo, Lacerda & Seoighe (2014), Terhorst et al. (2015) or Paris et al. (2019))
#' @param mnt_num the number of the Monte Carlo samples for the moment approximations

#' Standard version
approximateMoment_PyramidalHierarchicalBeta <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, ...) {
  if (mnt_apx == "MC") {
    param <- cmpcalculateParam_PyramidalHierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num)
  } else {
    param <- cmpcalculateParam_PyramidalHierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx)
  }
  mnt <- approximateMoment_PyramidalHierarchicalBeta_arma(param$alpha, param$beta, int_gen, lst_gen)

  return(list(mean = as.matrix(mnt$mean),
              variance = as.array(mnt$variance)))
}
#' Compiled version
cmpapproximateMoment_PyramidalHierarchicalBeta <- cmpfun(approximateMoment_PyramidalHierarchicalBeta)

########################################

#' Calculate the root mean square deviation between the two empirical cumulative distribution functions
#' Parameter setting
#' @param smp_mod the sample trajectories generated under the Wright-Fisher model
#' @param smp_apx the sample trajectories generated under the approximation of the Wright-Fisher model
#' @param grd_num the grid number for the empirical probability distribution function
#' @param rnd_grd = TRUE/FALSE (return the root mean square deviation between the empirical cumulative distribution functions with random grid points or not)

#' Standard version
calculateRMSD <- function(smp_mod, smp_apx, sim_num, grd_num, rnd_grd = TRUE) {
  if (rnd_grd == TRUE) {
    frq_grd <- generateRndGrid_2L_arma(grd_num)
  } else {
    frq_grd <- generateFixGrid_2L_arma(grd_num)
  }

  dist <- rep(NA, length.out = dim(smp_mod)[2])
  cdf_mod <- calculateECDF_2L_arma(smp_mod, sim_num, frq_grd)
  cdf_apx <- calculateECDF_2L_arma(smp_apx, sim_num, frq_grd)
  for (k in 1:dim(smp_mod)[2]) {
    dist[k] <- sqrt(sum((cdf_mod[, k] - cdf_apx[, k])^2) / grd_num)
  }

  return(dist)
}
#' Compiled version
cmpcalculateRMSD <- cmpfun(calculateRMSD)

########################################

#' Calculate the root mean square deviation between the two empirical cumulative distribution functions
#' Parameter setting
#' @param smp_mod the sample trajectories generated under the Wright-Fisher model

#' Standard version
generateInitFreq <- function(smp_siz) {
  int_frq <- generateInitFreq_arma(smp_siz)
  int_frq <- as.matrix(int_frq)

  return(int_frq)
}
#' Compiled version
cmpgenerateInitFreq <- cmpfun(generateInitFreq)

################################################################################
