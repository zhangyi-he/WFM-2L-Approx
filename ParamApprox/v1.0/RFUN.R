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

  fts_mat <- calculateFitnessMat_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B)

  # Approximate the first two moments of the Wright-Fisher model using
  if (mnt_apx == "MC") {
    # Monte Carlo simulation
    Moments <- approximateMoment_MonteCarlo_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num)
  }
  if (mnt_apx == "Lacerda") {
    # the extension of Lacerda & Seoighe (2014)
    Moments <- approximateMoment_Lacerda_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Terhorst") {
    # the extension of Terhorst et al. (2015)
    Moments <- approximateMoment_Terhorst_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Paris") {
    # the extension of Paris et al. (2019)
    Moments <- approximateMoment_Paris_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Terhorst2nd") {
    # the extension of Terhorst et al. (2015)
    Moments <- approximateMoment_Terhorst2ndOrder_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Paris2nd") {
    # the extension of Paris et al. (2019)
    Moments <- approximateMoment_Paris2ndOrder_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }

  mean = Moments$mean
  variance = Moments$variance

  return(list(mean = mean,
              variance = variance))
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

  fts_mat <- calculateFitnessMat_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B)

  # Approximate the first two moments of the Wright-Fisher model using
  if (mnt_apx == "MC") {
    # Monte Carlo simulation
    Moments <- approximateMoment_MonteCarlo_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num)
  }
  if (mnt_apx == "Lacerda") {
    # the extension of Lacerda & Seoighe (2014)
    Moments <- approximateMoment_Lacerda_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Terhorst") {
    # the extension of Terhorst et al. (2015)
    Moments <- approximateMoment_Terhorst_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Paris") {
    # the extension of Paris et al. (2019)
    Moments <- approximateMoment_Paris_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }

  mean = Moments$mean
  variance = Moments$variance

  Moments_norm = approximatWFM_norm_arma(mean, variance)
  mean_norm = Moments_norm$mean
  var_norm = Moments_norm$variance
  # return the parameters of the normal approximation and the corresponding first two moments
  return(list(mean = mean_norm,
              variance = var_norm))
}
#' Compiled version
cmpapproximateWFM_Norm <- cmpfun(approximateWFM_Norm)

########################################

#' Simulate the two-locus Wright-Fisher model with selection using the normal
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
simulateWFM_Norm <- function(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, ...) {
  sel_cof_A <- sel_cof[1]
  sel_cof_B <- sel_cof[2]
  dom_par_A <- dom_par[1]
  dom_par_B <- dom_par[2]
  fts_mat <- calculateFitnessMat_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B)

  # Approximate the first two moments of the Wright-Fisher model using
  if (mnt_apx == "MC") {
    # Monte Carlo simulation
    Moments <- approximateMoment_MonteCarlo_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num)
  }
  if (mnt_apx == "Lacerda") {
    # the extension of Lacerda & Seoighe (2014)
    Moments <- approximateMoment_Lacerda_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Terhorst") {
    # the extension of Terhorst et al. (2015)
    Moments <- approximateMoment_Terhorst_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Paris") {
    # the extension of Paris et al. (2019)
    Moments <- approximateMoment_Paris_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }

  mean = Moments$mean
  variance = Moments$variance

  Moments_norm = approximatWFM_norm_arma(mean, variance)
  mean_norm = Moments_norm$mean
  var_norm = Moments_norm$variance

  frq_pth = simulateWFM_norm_arma(mean_norm, var_norm, int_gen, lst_gen, smp_siz)
  return(frq_pth)
}
#' Compiled version
cmpsimulateWFM_Norm <- cmpfun(simulateWFM_Norm)

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

  fts_mat <- calculateFitnessMat_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B)

  # Approximate the first two moments of the Wright-Fisher model using
  if (mnt_apx == "MC") {
    # Monte Carlo simulation
    Moments <- approximateMoment_MonteCarlo_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num)
  }
  if (mnt_apx == "Lacerda") {
    # the extension of Lacerda & Seoighe (2014)
    Moments <- approximateMoment_Lacerda_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Terhorst") {
    # the extension of Terhorst et al. (2015)
    Moments <- approximateMoment_Terhorst_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Paris") {
    # the extension of Paris et al. (2019)
    Moments <- approximateMoment_Paris_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }

  mean = Moments$mean
  variance = Moments$variance

  Moments_loginorm = approximatWFM_LogisticNorm_arma(mean, variance, int_gen, lst_gen)
  location = Moments_loginorm$location
  squared_scale = Moments_loginorm$squared_scale
  # return the parameters of the logistic normal approximation and the corresponding first two moments
  return(list(location = location,
              squared_scale = squared_scale))
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

  fts_mat <- calculateFitnessMat_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B)

  # Approximate the first two moments of the Wright-Fisher model using
  if (mnt_apx == "MC") {
    # Monte Carlo simulation
    Moments <- approximateMoment_MonteCarlo_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num)
  }
  if (mnt_apx == "Lacerda") {
    # the extension of Lacerda & Seoighe (2014)
    Moments <- approximateMoment_Lacerda_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Terhorst") {
    # the extension of Terhorst et al. (2015)
    Moments <- approximateMoment_Terhorst_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Paris") {
    # the extension of Paris et al. (2019)
    Moments <- approximateMoment_Paris_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }

  mean = Moments$mean
  variance = Moments$variance

  Moments_loginorm = approximatWFM_LogisticNorm_arma(mean, variance, int_gen, lst_gen)
  location = Moments_loginorm$location
  squared_scale = Moments_loginorm$squared_scale
  frq_pth = simulateWFM_LogisticNorm_arma(location, squared_scale, int_gen, lst_gen, smp_siz)
  return(frq_pth)

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
  fts_mat <- calculateFitnessMat_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B)

  # Approximate the first two moments of the Wright-Fisher model using
  if (mnt_apx == "MC") {
    # Monte Carlo simulation
    Moments <- approximateMoment_MonteCarlo_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num)
  }
  if (mnt_apx == "Lacerda") {
    # the extension of Lacerda & Seoighe (2014)
    Moments <- approximateMoment_Lacerda_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Terhorst") {
    # the extension of Terhorst et al. (2015)
    Moments <- approximateMoment_Terhorst_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Paris") {
    # the extension of Paris et al. (2019)
    Moments <- approximateMoment_Paris_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }

  mean = Moments$mean
  variance = Moments$variance

  Moments_hierabeta = approximatWFM_HierarchicalBeta_arma(mean, variance, int_gen, lst_gen)
  alpha = Moments_hierabeta$alpha
  beta = Moments_hierabeta$beta

  # return the parameters of the hierarchical beta approximation and the corresponding first two moments
  return(list(alpha = alpha,
              beta = beta))
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
  fts_mat <- calculateFitnessMat_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B)

  # Approximate the first two moments of the Wright-Fisher model using
  if (mnt_apx == "MC") {
    # Monte Carlo simulation
    Moments <- approximateMoment_MonteCarlo_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num)
  }
  if (mnt_apx == "Lacerda") {
    # the extension of Lacerda & Seoighe (2014)
    Moments <- approximateMoment_Lacerda_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Terhorst") {
    # the extension of Terhorst et al. (2015)
    Moments <- approximateMoment_Terhorst_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }
  if (mnt_apx == "Paris") {
    # the extension of Paris et al. (2019)
    Moments <- approximateMoment_Paris_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_frq, int_gen, lst_gen, fts_mat)
  }

  mean = Moments$mean
  variance = Moments$variance

  Moments_hierabeta = approximatWFM_HierarchicalBeta_arma(mean, variance, int_gen, lst_gen)
  alpha = Moments_hierabeta$alpha
  beta = Moments_hierabeta$beta
  frq_pth = simulateWFM_HierarchicalBeta_arma(alpha, beta, int_gen, lst_gen, smp_siz)
  return(frq_pth)

}
#' Compiled version
cmpsimulateWFM_HierarchicalBeta <- cmpfun(simulateWFM_HierarchicalBeta)

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
