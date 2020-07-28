#' @title Numerical simulation of the two-locus Wright-Fisher diffusion with application to approximating transition probability densities
#' @author Zhangyi He, Mark Beaumont and Feng Yu

#' version 1.0

# set the directory
setwd("~/Documents/GitHub/WFM-2L-Approx")

source("./DiffusApprox/v1.0/RFUN_2L.R")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("viridis")
library("viridis")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

setwd("~/Dropbox/Jeffery He/iResearch/Publications/2017/HE2020-WFM-2L-DiffusApprox-TheorPopulBiol")

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

frq_pth <- cmpsimulateOLWFMS(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen)

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
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param data_augmentation = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_frq <- 2e-01
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00

frq_pth <- cmpsimulateOLWFDS(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = TRUE)

t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / pop_siz
plot(t, frq_pth, type = 'l', lwd = 1.5,
     xlab = "Generation", ylab = "Allele frequency",
     main = "A frequency trajectory of the mutant allele generated with the Wright-Fisher diffusion")

########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_frq <- 2e-01
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00
sim_num <- 1e+06

sim_frq_WFM <- numeric(sim_num)
sim_frq_WFD <- numeric(sim_num)
for (i in 1:sim_num) {
  print(i)
  sim_frq_WFM[i] <- cmpsimulateOLWFMS(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen)[(lst_gen - int_gen) + 1]
  sim_frq_WFD[i] <- cmpsimulateOLWFDS(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = FALSE)[(lst_gen - int_gen) + 1]
}

save(sel_cof, sel_cof, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num, sim_frq_WFM, sim_frq_WFD,
     file = "./Output/Output v2.1/Test v2.1/TEST_1L_WFM_vs_WFD.rda")

load("./Output/Output v2.1/Test v2.1/TEST_1L_WFM_vs_WFD.rda")

pdf(file = "./Output/Output v2.1/Test v2.1/TEST_1L_WFM_vs_WFD.pdf", width = 20, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
hist(sim_frq_WFM, breaks = seq(min(sim_frq_WFM, sim_frq_WFD), max(sim_frq_WFM, sim_frq_WFD), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(sim_frq_WFM, sim_frq_WFD), max(sim_frq_WFM, sim_frq_WFD)),
     xlab = "Allele frequency",
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Wright-Fisher diffusion"))
hist(sim_frq_WFD, breaks = seq(min(sim_frq_WFM, sim_frq_WFD), max(sim_frq_WFM, sim_frq_WFD), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
dev.off()

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
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00
sim_num <- 1e+06

sim_frq_WFM <- matrix(NA, nrow = 4, ncol = sim_num)
sim_frq_WFD <- matrix(NA, nrow = 4, ncol = sim_num)
for (i in 1:sim_num) {
  print(i)
  sim_frq_WFM[, i] <- cmpsimulateTLWFMS(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)[, (lst_gen - int_gen) + 1]
  sim_frq_WFD[, i] <- cmpsimulateTLWFDS(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = FALSE)[, (lst_gen - int_gen) + 1]
}

save(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num, sim_frq_WFM, sim_frq_WFD,
     file = "./Output/Output v2.1/Test v2.1/TEST_2L_WFM_vs_WFD.rda")

load("./Output/Output v2.1/Test v2.1/TEST_2L_WFM_vs_WFD.rda")

pdf(file = "./Output/Output v2.1/Test v2.1/TEST_2L_WFM_vs_WFD.pdf", width = 20, height = 10)
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
