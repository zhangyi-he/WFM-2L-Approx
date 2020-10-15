#' @title Moment-based approximations for the Wright-Fisher model of population dynamics under natural selection at two linked loci
#' @author Zhangyi He, Wenyang Lyu and Feng Yu

#' version 1.0

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2017/HE2020-WFM-2L-ParamApprox-TheorPopulBiol")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("viridis")
library("viridis")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

#' call R functions
source("./Code/Code v1.0/RFUN.R")

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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500





################################################################################

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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00





################################################################################

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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05





################################################################################

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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
sim_num <- 1e+05
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05





################################################################################

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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05





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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
sim_num <- 1e+05
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05





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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
sim_num <- 1e+05
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05





################################################################################

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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05





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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
sim_num <- 1e+05
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05





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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05





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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500

frq_pth <- cmpsimulateTLWFMS(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)

k <- int_gen:lst_gen
plot(k, frq_pth[1, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A1B1 haplotype generated with the Wright-Fisher model")
plot(k, frq_pth[2, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A1B2 haplotype generated with the Wright-Fisher model")
plot(k, frq_pth[3, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A2B1 haplotype generated with the Wright-Fisher model")
plot(k, frq_pth[4, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A2B2 haplotype generated with the Wright-Fisher model")

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
#' @param data_augmentation = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00

frq_pth <- cmpsimulateTLWFDS(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = TRUE)

t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / pop_siz
plot(t, frq_pth[1, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A1B1 haplotype generated with the Wright-Fisher diffusion")
plot(t, frq_pth[2, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A1B2 haplotype generated with the Wright-Fisher diffusion")
plot(t, frq_pth[3, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A2B1 haplotype generated with the Wright-Fisher diffusion")
plot(t, frq_pth[4, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A2B2 haplotype generated with the Wright-Fisher diffusion")

########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+02
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 50
ptn_num <- 5e+00
sim_num <- 1e+04

sim_frq_WFM <- matrix(NA, nrow = 4, ncol = sim_num)
sim_frq_WFD <- matrix(NA, nrow = 4, ncol = sim_num)
for (i in 1:sim_num) {
  print(i)
  sim_frq_WFM[, i] <- cmpsimulateWFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)[, (lst_gen - int_gen) + 1]
  sim_frq_WFD[, i] <- cmpsimulateDiffusApprox(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = FALSE)[, (lst_gen - int_gen) + 1]
}

save(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num, sim_frq_WFM, sim_frq_WFD,
     file = "./Output/Output v2.1/Test v2.1/TEST_2L_WFM_vs_WFD.rda")

load("./Output/Output v2.1/Test v2.1/TEST_2L_WFM_vs_WFD.rda")

pdf(file = "./Desktop/TEST_2L_WFM_vs_WFD.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), max(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), max(sim_frq_WFM[1, ], sim_frq_WFD[1, ])),
     xlab = "Haplotype frequency", main = "Haplotype A1B1")
hist(sim_frq_WFD[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), max(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), max(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), max(sim_frq_WFM[2, ], sim_frq_WFD[2, ])),
     xlab = "Haplotype frequency", main = "Haplotype A1B2")
hist(sim_frq_WFD[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), max(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), max(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), max(sim_frq_WFM[3, ], sim_frq_WFD[3, ])),
     xlab = "Haplotype frequency", main = "Haplotype A2B1")
hist(sim_frq_WFD[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), max(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), max(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), max(sim_frq_WFM[4, ], sim_frq_WFD[4, ])),
     xlab = "Haplotype frequency", main = "Haplotype A2B2")
hist(sim_frq_WFD[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), max(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
title(paste("Histograms of the haplotype frequencies in generation", lst_gen, "under the Wright-Fisher model and the Wright-Fisher diffusion"), outer = TRUE)
dev.off()

################################################################################

#' Compare the moments of the Wright-Fisher model using different methods
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
#' @param mnt_num the number of the Monte Carlo samples

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 100
sim_num <- 1e+05
mnt_num <- 5e+00

Moments_MC <- approximateMoment(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, "MC", sim_num)
mean_MC = Moments_MC$mean
var_MC = Moments_MC$variance
Moments_Lacerda <- approximateMoment(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, "Lacerda", sim_num)
mean_Lacerda = Moments_Lacerda$mean
var_Lacerda = Moments_Lacerda$variance
Moments_Terhorst <- approximateMoment(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, "Terhorst", sim_num)
mean_Terhorst = Moments_Terhorst$mean
var_Terhorst = Moments_Terhorst$variance
Moments_Paris <- approximateMoment(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, "Paris", sim_num)
mean_Paris = Moments_Paris$mean
var_Paris = Moments_Paris$variance
Moments_Terhorst2nd <- approximateMoment(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, "Terhorst2nd", sim_num)
mean_Terhorst2nd = Moments_Terhorst2nd$mean
var_Terhorst2nd = Moments_Terhorst2nd$variance
Moments_Paris2nd <- approximateMoment(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, "Paris2nd", sim_num)
mean_Paris2nd = Moments_Paris2nd$mean
var_Paris2nd = Moments_Paris2nd$variance
#################### mean
k <- int_gen:lst_gen
plot(k, mean_MC[1, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency", col = "red", ylim = c(0, 1),
     main = "")
lines(k, mean_Lacerda[1, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency", col = "green", ylim = c(0, 1),
     main = "")
lines(k, mean_Terhorst[1, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(0, 1),
     main = "")
lines(k, mean_Paris[1, ], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(0, 1),
      main = "")
lines(k, mean_Terhorst2nd[1, ], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(0, 1),
      main = "")
lines(k, mean_Paris2nd[1, ], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(0, 1),
      main = "")

lines(k, mean_MC[2, ], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "red", ylim = c(0, 1),
      main = "")
lines(k, mean_Lacerda[2, ], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "green", ylim = c(0, 1),
      main = "")
lines(k, mean_Terhorst[2, ], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(0, 1),
      main = "")
lines(k, mean_Paris[2, ], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(0, 1),
      main = "")
lines(k, mean_Terhorst2nd[2, ], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(0, 1),
      main = "")
lines(k, mean_Paris2nd[2, ], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(0, 1),
      main = "")

lines(k, mean_MC[3, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency", col = "red", ylim = c(0, 1),
     main = "")
lines(k, mean_Lacerda[3, ], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "green", ylim = c(0, 1),
      main = "")
lines(k, mean_Terhorst[3, ], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(0, 1),
      main = "")
lines(k, mean_Paris[3, ], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(0, 1),
      main = "")
lines(k, mean_Terhorst2nd[3, ], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(0, 1),
      main = "")
lines(k, mean_Paris2nd[3, ], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(0, 1),
      main = "")

lines(k, mean_MC[4, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency",col = "red", ylim = c(0, 1),
     main = "")
lines(k, mean_Lacerda[4, ], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "green", ylim = c(0, 1),
      main = "")
lines(k, mean_Terhorst[4, ], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "blue", ylim = c(0, 1),
      main = "")
lines(k, mean_Paris[4, ], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "yellow", ylim = c(0, 1),
      main = "")
lines(k, mean_Terhorst2nd[4, ], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "blue", ylim = c(0, 1),
      main = "")
lines(k, mean_Paris2nd[4, ], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "yellow", ylim = c(0, 1),
      main = "")
######################################## variance

plot(k, var_MC[1, 1,], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency", col = "red", ylim = c(0, .002),
     main = "")
lines(k, var_Lacerda[1, 1,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "green", ylim = c(0, .002),
      main = "")
lines(k, var_Terhorst[1, 1,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(0, .002),
      main = "")
lines(k, var_Paris[1, 1,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(0, .002),
      main = "")
lines(k, var_Terhorst2nd[1, 1,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(0, .002),
      main = "")
lines(k, var_Paris2nd[1, 1,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(0, .002),
      main = "")

lines(k, var_MC[2, 2,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "red", ylim = c(0, .002),
      main = "")
lines(k, var_Lacerda[2, 2,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "green", ylim = c(0, .002),
      main = "")
lines(k, var_Terhorst[2, 2,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(0, .002),
      main = "")
lines(k, var_Paris[2, 2,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(0, .002),
      main = "")
lines(k, var_Terhorst2nd[2, 2,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(0, .002),
      main = "")
lines(k, var_Paris2nd[2, 2,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(0, .002),
      main = "")

lines(k, var_MC[3, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "red", ylim = c(0, .002),
      main = "")
lines(k, var_Lacerda[3, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "green", ylim = c(0, .002),
      main = "")
lines(k, var_Terhorst[3, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(0, .002),
      main = "")
lines(k, var_Paris[3, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(0, .002),
      main = "")
lines(k, var_Terhorst2nd[3, 3,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(0, .002),
      main = "")
lines(k, var_Paris2nd[3, 3,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(0, .002),
      main = "")

lines(k, var_MC[4, 4,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "red", ylim = c(0, .002),
      main = "")
lines(k, var_Lacerda[4, 4,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "green", ylim = c(0, .002),
      main = "")
lines(k, var_Terhorst[4, 4,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "blue", ylim = c(0, .002),
      main = "")
lines(k, var_Paris[4, 4,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "yellow", ylim = c(0, .002),
      main = "")
lines(k, var_Terhorst2nd[4, 4,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "blue", ylim = c(0, .002),
      main = "")
lines(k, var_Paris2nd[4, 4,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "yellow", ylim = c(0, .002),
      main = "")

################################### covariance
plot(k, var_MC[1, 2,], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency", col = "red", ylim = c(-0.001, .0001),
     main = "")
lines(k, var_Lacerda[1, 2,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "green", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Terhorst[1, 2,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Paris[1, 2,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Terhorst2nd[1, 2,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Paris2nd[1, 2,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(-0.001, .0001),
      main = "")

lines(k, var_MC[1, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "red", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Lacerda[1, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "green", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Terhorst[1, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Paris[1, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Terhorst2nd[1, 3,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Paris2nd[1, 3,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(-0.001, .0001),
      main = "")

lines(k, var_MC[1, 4,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "red", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Lacerda[1, 4,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "green", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Terhorst[1, 4,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Paris[1, 4,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Terhorst2nd[1, 4,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "blue", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Paris2nd[1, 4,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency", col = "yellow", ylim = c(-0.001, .0001),
      main = "")

lines(k, var_MC[2, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "red", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Lacerda[2, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "green", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Terhorst[2, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "blue", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Paris[2, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "yellow", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Terhorst2nd[2, 3,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "blue", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Paris2nd[2, 3,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "yellow", ylim = c(-0.001, .0001),
      main = "")

lines(k, var_MC[2, 4,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "red", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Lacerda[2, 4,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "green", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Terhorst[2, 4,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "blue", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Paris[2, 4,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "yellow", ylim = c(-0.001, .0001),
      main = "")

lines(k, var_MC[4, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "red", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Lacerda[4, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "green", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Terhorst[4, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "blue", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Paris[4, 3,], type = "l", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "yellow", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Terhorst2nd[4, 3,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "blue", ylim = c(-0.001, .0001),
      main = "")
lines(k, var_Paris2nd[4, 3,], type = "p", lwd = 1.5,
      xlab = "Time", ylab = "Haplotype frequency",col = "yellow", ylim = c(-0.001, .0001),
      main = "")

################################################################################

#' Simulate the Wright-Fisher model with different methods using Normal approximation
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
#' @param mnt_num the number of the Monte Carlo samples

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 150
sim_num <- 1e+04
mnt_num <- 1e+04

frq_pth = simulateWFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
frq_pth_MC = simulateWFM_Norm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, "MC", mnt_num)
frq_pth_Lacerda = simulateWFM_Norm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, "Lacerda", mnt_num)
frq_pth_Terhorst = simulateWFM_Norm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, "Terhorst", mnt_num)
frq_pth_Paris = simulateWFM_Norm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, "Paris", mnt_num)

####################
sim_frq_WFM <- matrix(NA, 4, mnt_num)
sim_frq_MC <- matrix(NA, 4, mnt_num)
sim_frq_Lacerda <- matrix(NA, 4, mnt_num)
sim_frq_Terhorst <- matrix(NA, 4, mnt_num)
sim_frq_Paris <- matrix(NA, 4, mnt_num)

for (i in 1:sim_num) {
  print(i)
  sim_frq_WFM[, i] <- simulateWFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)[, (lst_gen - int_gen) + 1]
}
sim_frq_MC = frq_pth_MC[, (lst_gen - int_gen) + 1, ]
sim_frq_Lacerda <- frq_pth_Lacerda[, (lst_gen - int_gen) + 1, ]
sim_frq_Terhorst <- frq_pth_Terhorst[, (lst_gen - int_gen) + 1, ]
sim_frq_Paris <- frq_pth_Paris[, (lst_gen - int_gen) + 1, ]

save(sel_cof, sel_cof, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num, sim_frq_WFM, sim_frq_MC, sim_frq_Lacerda, sim_frq_Terhorst, sim_frq_Paris,
     file = "./Desktop/Code v1.0/TEST_1L_WFM_vs_Normal.rda")

load("./Desktop/Code v1.0/TEST_1L_WFM_vs_Normal.rda")

pdf(file = "./Desktop/Code v1.0/TEST_2L_WFM_vs_Normal_MC.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_MC[1, ]), max(sim_frq_WFM[1, ], sim_frq_MC[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_MC[1, ]), max(sim_frq_WFM[1, ], sim_frq_MC[1, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_MC[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_MC[1, ]), max(sim_frq_WFM[1, ], sim_frq_MC[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_MC[2, ]), max(sim_frq_WFM[2, ], sim_frq_MC[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_MC[2, ]), max(sim_frq_WFM[2, ], sim_frq_MC[2, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_MC[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_MC[2, ]), max(sim_frq_WFM[2, ], sim_frq_MC[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_MC[3, ]), max(sim_frq_WFM[3, ], sim_frq_MC[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_MC[3, ]), max(sim_frq_WFM[3, ], sim_frq_MC[3, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_MC[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_MC[3, ]), max(sim_frq_WFM[3, ], sim_frq_MC[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_MC[4, ]), max(sim_frq_WFM[4, ], sim_frq_MC[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_MC[4, ]), max(sim_frq_WFM[4, ], sim_frq_MC[4, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_MC[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_MC[4, ]), max(sim_frq_WFM[4, ], sim_frq_MC[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

dev.off()

pdf(file = "./Desktop/Code v1.0/TEST_2L_WFM_vs_Normal_Lacerda.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), max(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), max(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_Lacerda[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), max(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), max(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), max(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_Lacerda[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), max(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), max(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), max(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_Lacerda[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), max(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), max(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), max(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_Lacerda[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), max(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

dev.off()

pdf(file = "./Desktop/Code v1.0/TEST_2L_WFM_vs_Normal_Terhorst.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), max(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), max(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_Terhorst[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), max(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), max(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), max(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_Terhorst[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), max(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), max(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), max(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_Terhorst[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), max(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), max(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), max(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_Terhorst[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), max(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

dev.off()

pdf(file = "./Desktop/Code v1.0/TEST_2L_WFM_vs_Normal_Paris.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), max(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), max(sim_frq_WFM[1, ], sim_frq_Paris[1, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_Paris[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), max(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), max(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), max(sim_frq_WFM[2, ], sim_frq_Paris[2, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_Paris[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), max(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), max(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), max(sim_frq_WFM[3, ], sim_frq_Paris[3, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_Paris[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), max(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), max(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), max(sim_frq_WFM[4, ], sim_frq_Paris[4, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Normal approximation"))
hist(sim_frq_Paris[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), max(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

dev.off()


################################################################################

#' Simulate the Wright-Fisher model with different methods using LogitNormal approximation
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
#' @param mnt_num the number of the Monte Carlo samples

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 150
sim_num <- 1e+04
mnt_num <- 1e+04

frq_pth = simulateWFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
frq_pth_MC = simulateWFM_LogisticNorm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, "MC", mnt_num)
frq_pth_Lacerda = simulateWFM_LogisticNorm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, "Lacerda", mnt_num)
frq_pth_Terhorst = simulateWFM_LogisticNorm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, "Terhorst", mnt_num)
frq_pth_Paris = simulateWFM_LogisticNorm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, "Paris", mnt_num)

####################
sim_frq_WFM <- matrix(NA, 4, mnt_num)
sim_frq_MC <- matrix(NA, 4, mnt_num)
sim_frq_Lacerda <- matrix(NA, 4, mnt_num)
sim_frq_Terhorst <- matrix(NA, 4, mnt_num)
sim_frq_Paris <- matrix(NA, 4, mnt_num)

for (i in 1:sim_num) {
  print(i)
  sim_frq_WFM[, i] <- simulateWFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)[, (lst_gen - int_gen) + 1]
}
sim_frq_MC = frq_pth_MC[, (lst_gen - int_gen) + 1, ]
sim_frq_Lacerda <- frq_pth_Lacerda[, (lst_gen - int_gen) + 1, ]
sim_frq_Terhorst <- frq_pth_Terhorst[, (lst_gen - int_gen) + 1, ]
sim_frq_Paris <- frq_pth_Paris[, (lst_gen - int_gen) + 1, ]

save(sel_cof, sel_cof, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num, sim_frq_WFM, sim_frq_MC, sim_frq_Lacerda, sim_frq_Terhorst, sim_frq_Paris,
     file = "./Desktop/Code v1.0/TEST_1L_WFM_vs_LogitNormal.rda")

load("./Desktop/Code v1.0/TEST_1L_WFM_vs_LogitNormal.rda")

pdf(file = "./Desktop/Code v1.0/TEST_2L_WFM_vs_LogitNormal_MC.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_MC[1, ]), max(sim_frq_WFM[1, ], sim_frq_MC[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_MC[1, ]), max(sim_frq_WFM[1, ], sim_frq_MC[1, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_MC[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_MC[1, ]), max(sim_frq_WFM[1, ], sim_frq_MC[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_MC[2, ]), max(sim_frq_WFM[2, ], sim_frq_MC[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_MC[2, ]), max(sim_frq_WFM[2, ], sim_frq_MC[2, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_MC[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_MC[2, ]), max(sim_frq_WFM[2, ], sim_frq_MC[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_MC[3, ]), max(sim_frq_WFM[3, ], sim_frq_MC[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_MC[3, ]), max(sim_frq_WFM[3, ], sim_frq_MC[3, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_MC[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_MC[3, ]), max(sim_frq_WFM[3, ], sim_frq_MC[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_MC[4, ]), max(sim_frq_WFM[4, ], sim_frq_MC[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_MC[4, ]), max(sim_frq_WFM[4, ], sim_frq_MC[4, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_MC[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_MC[4, ]), max(sim_frq_WFM[4, ], sim_frq_MC[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

dev.off()

pdf(file = "./Desktop/Code v1.0/TEST_2L_WFM_vs_LogitNormal_Lacerda.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), max(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), max(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_Lacerda[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), max(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), max(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), max(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_Lacerda[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), max(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), max(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), max(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_Lacerda[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), max(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), max(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), max(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_Lacerda[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), max(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

dev.off()

pdf(file = "./Desktop/Code v1.0/TEST_2L_WFM_vs_LogitNormal_Terhorst.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), max(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), max(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_Terhorst[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), max(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), max(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), max(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_Terhorst[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), max(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), max(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), max(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_Terhorst[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), max(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), max(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), max(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_Terhorst[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), max(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

dev.off()

pdf(file = "./Desktop/Code v1.0/TEST_2L_WFM_vs_LogitNormal_Paris.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), max(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), max(sim_frq_WFM[1, ], sim_frq_Paris[1, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_Paris[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), max(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), max(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), max(sim_frq_WFM[2, ], sim_frq_Paris[2, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_Paris[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), max(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), max(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), max(sim_frq_WFM[3, ], sim_frq_Paris[3, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_Paris[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), max(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), max(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), max(sim_frq_WFM[4, ], sim_frq_Paris[4, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the LogitNormal approximation"))
hist(sim_frq_Paris[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), max(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

dev.off()

#' Simulate the Wright-Fisher model with different methods using HierarchicalBeta approximation
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
#' @param mnt_num the number of the Monte Carlo samples

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 200
sim_num <- 1e+04
mnt_num <- 1e+04

frq_pth_MC = simulateWFM_HierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, "MC", mnt_num)
frq_pth_Lacerda = simulateWFM_HierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, "Lacerda", mnt_num)
frq_pth_Terhorst = simulateWFM_HierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, "Terhorst", mnt_num)
frq_pth_Paris = simulateWFM_HierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, "Paris", mnt_num)

####################
sim_frq_WFM <- matrix(NA, 4, mnt_num)
sim_frq_MC <- matrix(NA, 4, mnt_num)
sim_frq_Lacerda <- matrix(NA, 4, mnt_num)
sim_frq_Terhorst <- matrix(NA, 4, mnt_num)
sim_frq_Paris <- matrix(NA, 4, mnt_num)

for (i in 1:sim_num) {
  print(i)
  sim_frq_WFM[, i] <- simulateWFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)[, (lst_gen - int_gen) + 1]
}
sim_frq_MC = frq_pth_MC[, (lst_gen - int_gen) + 1, ]
sim_frq_Lacerda <- frq_pth_Lacerda[, (lst_gen - int_gen) + 1, ]
sim_frq_Terhorst <- frq_pth_Terhorst[, (lst_gen - int_gen) + 1, ]
sim_frq_Paris <- frq_pth_Paris[, (lst_gen - int_gen) + 1, ]

save(sel_cof, sel_cof, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num, sim_frq_WFM, sim_frq_MC, sim_frq_Lacerda, sim_frq_Terhorst, sim_frq_Paris,
     file = "./Desktop/Code v1.0/TEST_1L_WFM_vs_HieraBeta.rda")

load("./Desktop/Code v1.0/TEST_1L_WFM_vs_HieraBeta.rda")

pdf(file = "./Desktop/Code v1.0/TEST_2L_WFM_vs_HieraBeta_MC.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_MC[1, ]), max(sim_frq_WFM[1, ], sim_frq_MC[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_MC[1, ]), max(sim_frq_WFM[1, ], sim_frq_MC[1, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_MC[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_MC[1, ]), max(sim_frq_WFM[1, ], sim_frq_MC[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_MC[2, ]), max(sim_frq_WFM[2, ], sim_frq_MC[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_MC[2, ]), max(sim_frq_WFM[2, ], sim_frq_MC[2, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_MC[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_MC[2, ]), max(sim_frq_WFM[2, ], sim_frq_MC[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_MC[3, ]), max(sim_frq_WFM[3, ], sim_frq_MC[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_MC[3, ]), max(sim_frq_WFM[3, ], sim_frq_MC[3, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_MC[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_MC[3, ]), max(sim_frq_WFM[3, ], sim_frq_MC[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_MC[4, ]), max(sim_frq_WFM[4, ], sim_frq_MC[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_MC[4, ]), max(sim_frq_WFM[4, ], sim_frq_MC[4, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_MC[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_MC[4, ]), max(sim_frq_WFM[4, ], sim_frq_MC[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

dev.off()

pdf(file = "./Desktop/Code v1.0/TEST_2L_WFM_vs_HieraBeta_Lacerda.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), max(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), max(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_Lacerda[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), max(sim_frq_WFM[1, ], sim_frq_Lacerda[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), max(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), max(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_Lacerda[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), max(sim_frq_WFM[2, ], sim_frq_Lacerda[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), max(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), max(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_Lacerda[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), max(sim_frq_WFM[3, ], sim_frq_Lacerda[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), max(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), max(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_Lacerda[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), max(sim_frq_WFM[4, ], sim_frq_Lacerda[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

dev.off()

pdf(file = "./Desktop/Code v1.0/TEST_2L_WFM_vs_HieraBeta_Terhorst.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), max(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), max(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_Terhorst[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), max(sim_frq_WFM[1, ], sim_frq_Terhorst[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), max(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), max(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_Terhorst[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), max(sim_frq_WFM[2, ], sim_frq_Terhorst[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), max(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), max(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_Terhorst[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), max(sim_frq_WFM[3, ], sim_frq_Terhorst[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), max(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), max(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_Terhorst[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), max(sim_frq_WFM[4, ], sim_frq_Terhorst[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

dev.off()

pdf(file = "./Desktop/Code v1.0/TEST_2L_WFM_vs_HieraBeta_Paris.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), max(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), max(sim_frq_WFM[1, ], sim_frq_Paris[1, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_Paris[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), max(sim_frq_WFM[1, ], sim_frq_Paris[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), max(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), max(sim_frq_WFM[2, ], sim_frq_Paris[2, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_Paris[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), max(sim_frq_WFM[2, ], sim_frq_Paris[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), max(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), max(sim_frq_WFM[3, ], sim_frq_Paris[3, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_Paris[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), max(sim_frq_WFM[3, ], sim_frq_Paris[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), max(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), max(sim_frq_WFM[4, ], sim_frq_Paris[4, ])),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the HierarchicalBeta approximation"))
hist(sim_frq_Paris[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), max(sim_frq_WFM[4, ], sim_frq_Paris[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

dev.off()

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+02
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 5
ptn_num <- 5e+02
sim_num <- 1e+04

sim_frq_WFM <- array(NA, c(sim_num, 4, lst_gen - int_gen + 1))
sim_frq_WFD <- array(NA, c(sim_num, 4, lst_gen - int_gen + 1))
for (i in 1:sim_num) {
  print(i)
  sim_frq_WFM[i, ,] <- cmpsimulateWFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
  sim_frq_WFD[i, ,] <- cmpsimulateDiffusApprox(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = FALSE)
}

dist = cmpcalculateRMSD_2nd(sim_frq_WFM, sim_frq_WFD, int_gen, lst_gen, sim_num, 100, rnd_grd = TRUE)

