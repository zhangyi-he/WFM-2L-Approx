#' @title Numerical simulation of the two-locus Wright-Fisher diffusion with application to approximating transition probability densities
#' @author Zhangyi He, Mark Beaumont and Feng Yu

#' version 1.0

# set the directory
setwd("~/Documents/GitHub/WFM-2L-Approx")

source("./DiffusApprox/v1.0/RFUN.R")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("viridis")
library("viridis")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

setwd("~/Dropbox/Jeffery He/iResearch/Publications/2012/HE2020-WFM-2L-DiffusApprox-TheorPopulBiol")

################################################################################

#' Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial mutant allele frequency of the population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory

sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_frq <- 2e-01
int_gen <- 0
lst_gen <- 500

frq_pth <- cmpsimulateWFM_1L(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen)

k <- int_gen:lst_gen
plot(k, frq_pth, type = 'l', lwd = 1.5,
     xlab = "Generation", ylab = "Allele frequency",
     main = "A frequency trajectory of the mutant allele generated with the Wright-Fisher model")

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

sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_frq <- 2e-01
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00

frq_pth <- simulateWFD_1L(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE)

t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / pop_siz
plot(t, frq_pth, type = 'l', lwd = 1.5,
     xlab = "Time", ylab = "Allele frequency",
     main = "A frequency trajectory of the mutant allele generated with the Wright-Fisher diffusion")

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

sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_frq <- 2e-01
int_gen <- 0
lst_gen <- 500
sim_num <- 5e+05
ptn_num <- 5e+00

model <- "WFM"
smp_WFM <- cmpgenerateSampleTraj_1L(model, sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, sim_num)[, (lst_gen - int_gen) + 1]
model <- "WFD"
smp_WFD <- cmpgenerateSampleTraj_1L(model, sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, sim_num, ptn_num)[, (lst_gen - int_gen) + 1]

hist(smp_WFM, breaks = seq(min(smp_WFM, smp_WFD), max(smp_WFM, smp_WFD), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_WFM, smp_WFD), max(smp_WFM, smp_WFD)),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen))
hist(smp_WFD, breaks = seq(min(smp_WFM, smp_WFD), max(smp_WFM, smp_WFD), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

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

frq_pth <- cmpsimulateWFM_2L(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)

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
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00

frq_pth <- cmpsimulateWFD_2L(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE)

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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
sim_num <- 5e+05
ptn_num <- 5e+00

model <- "WFM"
smp_WFM <- cmpgenerateSampleTraj_2L(model, sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num)[, , (lst_gen - int_gen) + 1]
model <- "WFD"
smp_WFD <- cmpgenerateSampleTraj_2L(model, sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, ptn_num)[, , (lst_gen - int_gen) + 1]

hist(smp_WFM[, 1], breaks = seq(min(smp_WFM[, 1], smp_WFD[, 1]), max(smp_WFM[, 1], smp_WFD[, 1]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_WFM[, 1], smp_WFD[, 1]), max(smp_WFM[, 1], smp_WFD[, 1])),
     xlab = "Haplotype frequency", main = paste("Histogram of the A1B1 haplotype", lst_gen))
hist(smp_WFD[, 1], breaks = seq(min(smp_WFM[, 1], smp_WFD[, 1]), max(smp_WFM[, 1], smp_WFD[, 1]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(smp_WFM[, 2], breaks = seq(min(smp_WFM[, 2], smp_WFD[, 2]), max(smp_WFM[, 2], smp_WFD[, 2]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_WFM[, 2], smp_WFD[, 2]), max(smp_WFM[, 2], smp_WFD[, 2])),
     xlab = "Haplotype frequency", main = paste("Histogram of the A1B2 haplotype", lst_gen))
hist(smp_WFD[, 2], breaks = seq(min(smp_WFM[, 2], smp_WFD[, 2]), max(smp_WFM[, 2], smp_WFD[, 2]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(smp_WFM[, 3], breaks = seq(min(smp_WFM[, 3], smp_WFD[, 3]), max(smp_WFM[, 3], smp_WFD[, 3]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_WFM[, 3], smp_WFD[, 3]), max(smp_WFM[, 3], smp_WFD[, 3])),
     xlab = "Haplotype frequency", main = paste("Histogram of the A2B1 haplotype", lst_gen))
hist(smp_WFD[, 3], breaks = seq(min(smp_WFM[, 3], smp_WFD[, 3]), max(smp_WFM[, 3], smp_WFD[, 3]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(smp_WFM[, 4], breaks = seq(min(smp_WFM[, 4], smp_WFD[, 4]), max(smp_WFM[, 4], smp_WFD[, 4]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_WFM[, 4], smp_WFD[, 4]), max(smp_WFM[, 4], smp_WFD[, 4])),
     xlab = "Haplotype frequency", main = paste("Histogram of the A2B2 haplotype", lst_gen))
hist(smp_WFD[, 4], breaks = seq(min(smp_WFM[, 4], smp_WFD[, 4]), max(smp_WFM[, 4], smp_WFD[, 4]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

################################################################################

#' Calculate the Hellinger distance between the empirical probability distribution functions for the Wright-Fisher model/diffusion
#' Parameter setting
#' @param smp_WFM the sample trajectories generated under the Wright-Fisher model
#' @param smp_WFD the sample trajectories generated under the Wright-Fisher diffusion
#' @param grd_num the grid number for the empirical probability distribution function

sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_frq <- 2e-01
int_gen <- 0
lst_gen <- 500
sim_num <- 5e+05
ptn_num <- 5e+00

model <- "WFM"
system.time(smp_WFM <- cmpgenerateSampleTraj_1L(model, sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, sim_num))
model <- "WFD"
system.time(smp_WFD <- cmpgenerateSampleTraj_1L(model, sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, sim_num, ptn_num))
grd_num <- 2 * pop_siz

system.time(dist <- cmpcalculateHellingerDist(smp_WFM, smp_WFD, grd_num))

plot(int_gen:lst_gen, dist, type = "l", xlab = "Generation", ylab = "Hellinger distance",
     main = "Hellinger distance between the empirical probability distribution functions for the mutant allele")


sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
sim_num <- 5e+05
ptn_num <- 5e+00

model <- "WFM"
system.time(smp_WFM <- cmpgenerateSampleTraj_2L(model, sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num))
model <- "WFD"
system.time(smp_WFD <- cmpgenerateSampleTraj_2L(model, sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, ptn_num))
grd_num <- 2 * pop_siz

system.time(dist <- cmpcalculateHellingerDist(smp_WFM, smp_WFD, grd_num))

plot(int_gen:lst_gen, dist[1, ], type = "l", xlab = "Generation", ylab = "Hellinger distance",
     main = "Hellinger distance between the empirical probability distribution functions for the A1B1 haplotype")
plot(int_gen:lst_gen, dist[2, ], type = "l", xlab = "Generation", ylab = "Hellinger distance",
     main = "Hellinger distance between the empirical probability distribution functions for the A1B2 haplotype")
plot(int_gen:lst_gen, dist[3, ], type = "l", xlab = "Generation", ylab = "Hellinger distance",
     main = "Hellinger distance between the empirical probability distribution functions for the A2B1 haplotype")
plot(int_gen:lst_gen, dist[4, ], type = "l", xlab = "Generation", ylab = "Hellinger distance",
     main = "Hellinger distance between the empirical probability distribution functions for the A2B2 haplotype")

################################################################################
