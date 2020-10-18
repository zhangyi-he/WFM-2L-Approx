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

# sel_cof <- c(1e-02, 5e-03)
# dom_par <- c(5e-01, 5e-01)
# rec_rat <- 1e-05
# pop_siz <- 5e+03
# int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
# int_gen <- 0
# lst_gen <- 500
#
# frq_pth <- cmpsimulateWFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)
#
# k <- int_gen:lst_gen
# plot(k, frq_pth[1, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Haplotype frequency",
#      main = "WFM: the A1B1 haplotype frequency trajectory")
# plot(k, frq_pth[2, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Haplotype frequency",
#      main = "WFM: the A1B2 haplotype frequency trajectory")
# plot(k, frq_pth[3, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Haplotype frequency",
#      main = "WFM: the A2B1 haplotype frequency trajectory")
# plot(k, frq_pth[4, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Haplotype frequency",
#      main = "WFM: the A2B2 haplotype frequency trajectory")

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
rec_rat <- 1e-05
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00

# frq_pth <- cmpsimulateDiffusApprox(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE)
#
# t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / pop_siz
# plot(t, frq_pth[1, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Haplotype frequency",
#      main = "WFD: the A1B1 haplotype frequency trajectory")
# plot(t, frq_pth[2, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Haplotype frequency",
#      main = "WFD: the A1B2 haplotype frequency trajectory")
# plot(t, frq_pth[3, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Haplotype frequency",
#      main = "WFD: the A2B1 haplotype frequency trajectory")
# plot(t, frq_pth[4, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Haplotype frequency",
#      main = "WFD: the A2B2 haplotype frequency trajectory")

sim_num <- 1e+05
smp_mod <- array(NA, dim = c(4, lst_gen - int_gen, sim_num))
smp_apx <- array(NA, dim = c(4, lst_gen - int_gen, sim_num))
for (i in 1:sim_num) {
  print(i)
  smp_mod[, , i] <- cmpsimulateWFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)[, 2:(lst_gen - int_gen + 1)]
  smp_apx[, , i] <- cmpsimulateDiffusApprox(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = FALSE)[, 2:(lst_gen - int_gen + 1)]
}

save(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num, smp_mod, smp_apx,
     file = "./Output/Output v1.0/Test v1.0/TEST_Diffus.rda")

load("./Output/Output v1.0/Test v1.0/TEST_Diffus.rda")

gen <- 500
smp_mod <- smp_mod[, gen, ]
smp_apx <- smp_apx[, gen, ]

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_Diffus.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist_mod <- hist(smp_mod[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), plot = FALSE)
hist(smp_mod[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A1B1 haplotype at generation", gen))
hist(smp_apx[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist_mod <- hist(smp_mod[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), plot = FALSE)
hist(smp_mod[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A1B2 haplotype at generation", gen))
hist(smp_apx[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist_mod <- hist(smp_mod[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), plot = FALSE)
hist(smp_mod[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A2B1 haplotype at generation", gen))
hist(smp_apx[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist_mod <- hist(smp_mod[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), plot = FALSE)
hist(smp_mod[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A2B2 haplotype at generation", gen))
hist(smp_apx[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
dev.off()

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
rec_rat <- 1e-05
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05

mean <- array(NA, dim = c(5, 4, lst_gen - int_gen))
variance <- array(NA, dim = c(5, 4, lst_gen - int_gen))
covariance <- array(NA, dim = c(5, 6, lst_gen - int_gen))
system.time(mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx[1], mnt_num))
mean[1, , ] <- mnt$mean
variance[1, 1, ] <- mnt$variance[1, 1, ]
variance[1, 2, ] <- mnt$variance[2, 2, ]
variance[1, 3, ] <- mnt$variance[3, 3, ]
variance[1, 4, ] <- mnt$variance[4, 4, ]
covariance[1, 1, ] <- mnt$variance[1, 2, ]
covariance[1, 2, ] <- mnt$variance[1, 3, ]
covariance[1, 3, ] <- mnt$variance[2, 3, ]
covariance[1, 4, ] <- mnt$variance[1, 4, ]
covariance[1, 5, ] <- mnt$variance[2, 4, ]
covariance[1, 6, ] <- mnt$variance[3, 4, ]
system.time(mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx[2]))
mean[2, , ] <- mnt$mean
variance[2, 1, ] <- mnt$variance[1, 1, ]
variance[2, 2, ] <- mnt$variance[2, 2, ]
variance[2, 3, ] <- mnt$variance[3, 3, ]
variance[2, 4, ] <- mnt$variance[4, 4, ]
covariance[2, 1, ] <- mnt$variance[1, 2, ]
covariance[2, 2, ] <- mnt$variance[1, 3, ]
covariance[2, 3, ] <- mnt$variance[2, 3, ]
covariance[2, 4, ] <- mnt$variance[1, 4, ]
covariance[2, 5, ] <- mnt$variance[2, 4, ]
covariance[2, 6, ] <- mnt$variance[3, 4, ]
system.time(mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx[3]))
mean[3, , ] <- mnt$mean
variance[3, 1, ] <- mnt$variance[1, 1, ]
variance[3, 2, ] <- mnt$variance[2, 2, ]
variance[3, 3, ] <- mnt$variance[3, 3, ]
variance[3, 4, ] <- mnt$variance[4, 4, ]
covariance[3, 1, ] <- mnt$variance[1, 2, ]
covariance[3, 2, ] <- mnt$variance[1, 3, ]
covariance[3, 3, ] <- mnt$variance[2, 3, ]
covariance[3, 4, ] <- mnt$variance[1, 4, ]
covariance[3, 5, ] <- mnt$variance[2, 4, ]
covariance[3, 6, ] <- mnt$variance[3, 4, ]
system.time(mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx[4]))
mean[4, , ] <- mnt$mean
variance[4, 1, ] <- mnt$variance[1, 1, ]
variance[4, 2, ] <- mnt$variance[2, 2, ]
variance[4, 3, ] <- mnt$variance[3, 3, ]
variance[4, 4, ] <- mnt$variance[4, 4, ]
covariance[4, 1, ] <- mnt$variance[1, 2, ]
covariance[4, 2, ] <- mnt$variance[1, 3, ]
covariance[4, 3, ] <- mnt$variance[2, 3, ]
covariance[4, 4, ] <- mnt$variance[1, 4, ]
covariance[4, 5, ] <- mnt$variance[2, 4, ]
covariance[4, 6, ] <- mnt$variance[3, 4, ]
system.time(mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx[5]))
mean[5, , ] <- mnt$mean
variance[5, 1, ] <- mnt$variance[1, 1, ]
variance[5, 2, ] <- mnt$variance[2, 2, ]
variance[5, 3, ] <- mnt$variance[3, 3, ]
variance[5, 4, ] <- mnt$variance[4, 4, ]
covariance[5, 1, ] <- mnt$variance[1, 2, ]
covariance[5, 2, ] <- mnt$variance[1, 3, ]
covariance[5, 3, ] <- mnt$variance[2, 3, ]
covariance[5, 4, ] <- mnt$variance[1, 4, ]
covariance[5, 5, ] <- mnt$variance[2, 4, ]
covariance[5, 6, ] <- mnt$variance[3, 4, ]

save(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num, mean, variance, covariance,
     file = "./Output/Output v1.0/Test v1.0/TEST_MomentApprox.rda")

load("./Output/Output v1.0/Test v1.0/TEST_MomentApprox.rda")

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_MomentApprox.pdf", width = 8, height = 18)
par(mfrow = c(3, 1), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
k <- (int_gen + 1):lst_gen
plot(0, type = "n", xlim = c(int_gen + 1, lst_gen), ylim = c(min(mean), max(mean)),
     xlab = "Generation", ylab = "Mean",
     main = "Moment approximation of the Wright-Fisher model: mean")
for (i in 1:5) {
  for (j in 1:4) {
    lines(k, mean[i, j, ], type = "l", lty = i, lwd = 1.5, col = brewer.pal(n = 4, name = 'Set1')[j])
  }
}
legend("topleft", legend = mnt_apx, col = "black", lty = 1:5, lwd = 1.5, bty = "n", cex = 1.5)
legend("bottomleft", legend = c("A1B1", "A1B2", "A2B1", "A2B2"), col = brewer.pal(n = 4, name = 'Set1'), lty = 1, lwd = 1.5, bty = "n", cex = 1.5)
plot(0, type = "n", xlim = c(int_gen + 1, lst_gen), ylim = c(min(variance), max(variance)),
     xlab = "Generation", ylab = "Variance",
     main = "Moment approximation of the Wright-Fisher model: variance")
for (i in 1:5) {
  for (j in 1:4) {
    lines(k, variance[i, j, ], type = "l", lty = i, lwd = 1.5, col = brewer.pal(n = 4, name = 'Set1')[j])
  }
}
legend("topleft", legend = mnt_apx, col = "black", lty = 1:5, lwd = 1.5, bty = "n", cex = 1.5)
legend("bottomleft", legend = c("A1B1", "A1B2", "A2B1", "A2B2"), col = brewer.pal(n = 4, name = 'Set1'), lty = 1, lwd = 1.5, bty = "n", cex = 1.5)
plot(0, type = "n", xlim = c(int_gen + 1, lst_gen), ylim = c(min(covariance), max(covariance)),
     xlab = "Generation", ylab = "Covariance",
     main = "Moment approximation of the Wright-Fisher model: covariance")
for (i in 1:5) {
  for (j in 1:6) {
    lines(k, covariance[i, j, ], type = "l", lty = i, lwd = 1.5, col = brewer.pal(n = 6, name = 'Set1')[j])
  }
}
legend("topleft", legend = mnt_apx, col = "black", lty = 1:5, lwd = 1.5, bty = "n", cex = 1.5)
legend("bottomleft", legend = c("A1B1/A1B2", "A1B1/A2B1", "A1B2/A2B1", "A1B1/A2B2", "A1B2/A2B2", "A2B1/A2B2"), col = brewer.pal(n = 6, name = 'Set1'), lty = 1, lwd = 1.5, bty = "n", cex = 1.5)
dev.off()

################################################################################

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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-05
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05





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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-05
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
sim_num <- 1e+05
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05

system.time(smp_apx <- cmpgenerateSample_Norm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx[1], mnt_num))

smp_mod <- array(NA, dim = c(4, lst_gen - int_gen, sim_num))
for (i in 1:sim_num) {
  print(i)
  smp_mod[, , i] <- cmpsimulateWFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)[, 2:(lst_gen - int_gen + 1)]
}

save(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, mnt_num, smp_mod, smp_apx,
     file = "./Output/Output v1.0/Test v1.0/TEST_Norm_Approx.rda")

load("./Output/Output v1.0/Test v1.0/TEST_Norm_Approx.rda")

gen <- 500
smp_mod <- smp_mod[, gen, ]
smp_apx <- smp_apx[, gen, ]

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_Norm_Approx.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist_mod <- hist(smp_mod[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), plot = FALSE)
hist(smp_mod[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A1B1 haplotype at generation", gen))
hist(smp_apx[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist_mod <- hist(smp_mod[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), plot = FALSE)
hist(smp_mod[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A1B2 haplotype at generation", gen))
hist(smp_apx[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist_mod <- hist(smp_mod[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), plot = FALSE)
hist(smp_mod[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A2B1 haplotype at generation", gen))
hist(smp_apx[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist_mod <- hist(smp_mod[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), plot = FALSE)
hist(smp_mod[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A2B2 haplotype at generation", gen))
hist(smp_apx[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
dev.off()

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

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-05
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
sim_num <- 1e+05
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05

mean <- array(NA, dim = c(2, 4, lst_gen - int_gen))
variance <- array(NA, dim = c(2, 4, lst_gen - int_gen))
covariance <- array(NA, dim = c(2, 6, lst_gen - int_gen))
system.time(mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx[1], mnt_num))
mean[1, , ] <- mnt$mean
variance[1, 1, ] <- mnt$variance[1, 1, ]
variance[1, 2, ] <- mnt$variance[2, 2, ]
variance[1, 3, ] <- mnt$variance[3, 3, ]
variance[1, 4, ] <- mnt$variance[4, 4, ]
covariance[1, 1, ] <- mnt$variance[1, 2, ]
covariance[1, 2, ] <- mnt$variance[1, 3, ]
covariance[1, 3, ] <- mnt$variance[2, 3, ]
covariance[1, 4, ] <- mnt$variance[1, 4, ]
covariance[1, 5, ] <- mnt$variance[2, 4, ]
covariance[1, 6, ] <- mnt$variance[3, 4, ]
system.time(mnt <- cmpapproximateMoment_Norm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx[1], mnt_num))
mean[2, , ] <- mnt$mean
variance[2, 1, ] <- mnt$variance[1, 1, ]
variance[2, 2, ] <- mnt$variance[2, 2, ]
variance[2, 3, ] <- mnt$variance[3, 3, ]
variance[2, 4, ] <- mnt$variance[4, 4, ]
covariance[2, 1, ] <- mnt$variance[1, 2, ]
covariance[2, 2, ] <- mnt$variance[1, 3, ]
covariance[2, 3, ] <- mnt$variance[2, 3, ]
covariance[2, 4, ] <- mnt$variance[1, 4, ]
covariance[2, 5, ] <- mnt$variance[2, 4, ]
covariance[2, 6, ] <- mnt$variance[3, 4, ]

save(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num, mean, variance, covariance,
     file = "./Output/Output v1.0/Test v1.0/TEST_Norm_Moment.rda")

load("./Output/Output v1.0/Test v1.0/TEST_Norm_Moment.rda")

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_Norm_Moment.pdf", width = 8, height = 18)
par(mfrow = c(3, 1), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
k <- (int_gen + 1):lst_gen
plot(0, type = "n", xlim = c(int_gen + 1, lst_gen), ylim = c(min(mean), max(mean)),
     xlab = "Generation", ylab = "Mean",
     main = "Moment approximation of the normal approximation: mean")
for (i in 1:2) {
  for (j in 1:4) {
    lines(k, mean[i, j, ], type = "l", lty = i, lwd = 1.5, col = brewer.pal(n = 4, name = 'Set1')[j])
  }
}
legend("topleft", legend = c("Wright-Fisher", "normal"), col = "black", lty = 1:2, lwd = 1.5, bty = "n", cex = 1.5)
legend("bottomleft", legend = c("A1B1", "A1B2", "A2B1", "A2B2"), col = brewer.pal(n = 4, name = 'Set1'), lty = 1, lwd = 1.5, bty = "n", cex = 1.5)
plot(0, type = "n", xlim = c(int_gen + 1, lst_gen), ylim = c(min(variance), max(variance)),
     xlab = "Generation", ylab = "Variance",
     main = "Moment approximation of the normal approximation: variance")
for (i in 1:2) {
  for (j in 1:4) {
    lines(k, variance[i, j, ], type = "l", lty = i, lwd = 1.5, col = brewer.pal(n = 4, name = 'Set1')[j])
  }
}
legend("topleft", legend = c("Wright-Fisher", "normal"), col = "black", lty = 1:2, lwd = 1.5, bty = "n", cex = 1.5)
legend("bottomleft", legend = c("A1B1", "A1B2", "A2B1", "A2B2"), col = brewer.pal(n = 4, name = 'Set1'), lty = 1, lwd = 1.5, bty = "n", cex = 1.5)
plot(0, type = "n", xlim = c(int_gen + 1, lst_gen), ylim = c(min(covariance), max(covariance)),
     xlab = "Generation", ylab = "Covariance",
     main = "Moment approximation of the normal approximation: covariance")
for (i in 1:2) {
  for (j in 1:6) {
    lines(k, covariance[i, j, ], type = "l", lty = i, lwd = 1.5, col = brewer.pal(n = 6, name = 'Set1')[j])
  }
}
legend("topleft", legend = c("Wright-Fisher", "normal"), col = "black", lty = 1:2, lwd = 1.5, bty = "n", cex = 1.5)
legend("bottomleft", legend = c("A1B1/A1B2", "A1B1/A2B1", "A1B2/A2B1", "A1B1/A2B2", "A1B2/A2B2", "A2B1/A2B2"), col = brewer.pal(n = 6, name = 'Set1'), lty = 1, lwd = 1.5, bty = "n", cex = 1.5)
dev.off()

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
rec_rat <- 1e-05
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
rec_rat <- 1e-05
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
sim_num <- 1e+05
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05

system.time(smp_apx <- generateSample_LogisticNorm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx[1], mnt_num))

smp_mod <- array(NA, dim = c(4, lst_gen - int_gen, sim_num))
for (i in 1:sim_num) {
  print(i)
  smp_mod[, , i] <- cmpsimulateWFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)[, 2:(lst_gen - int_gen + 1)]
}

save(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, mnt_num, smp_mod, smp_apx,
     file = "./Output/Output v1.0/Test v1.0/TEST_LogisticNorm_Approx.rda")

load("./Output/Output v1.0/Test v1.0/TEST_LogisticNorm_Approx.rda")

gen <- 500
smp_mod <- smp_mod[, gen, ]
smp_apx <- smp_apx[, gen, ]

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_LogisticNorm_Approx.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist_mod <- hist(smp_mod[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), plot = FALSE)
hist(smp_mod[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A1B1 haplotype at generation", gen))
hist(smp_apx[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist_mod <- hist(smp_mod[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), plot = FALSE)
hist(smp_mod[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A1B2 haplotype at generation", gen))
hist(smp_apx[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist_mod <- hist(smp_mod[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), plot = FALSE)
hist(smp_mod[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A2B1 haplotype at generation", gen))
hist(smp_apx[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist_mod <- hist(smp_mod[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), plot = FALSE)
hist(smp_mod[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A2B2 haplotype at generation", gen))
hist(smp_apx[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
dev.off()

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
rec_rat <- 1e-05
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
sim_num <- 1e+05
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05

mean <- array(NA, dim = c(2, 4, lst_gen - int_gen))
variance <- array(NA, dim = c(2, 4, lst_gen - int_gen))
covariance <- array(NA, dim = c(2, 6, lst_gen - int_gen))
system.time(mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx[1], mnt_num))
mean[1, , ] <- mnt$mean
variance[1, 1, ] <- mnt$variance[1, 1, ]
variance[1, 2, ] <- mnt$variance[2, 2, ]
variance[1, 3, ] <- mnt$variance[3, 3, ]
variance[1, 4, ] <- mnt$variance[4, 4, ]
covariance[1, 1, ] <- mnt$variance[1, 2, ]
covariance[1, 2, ] <- mnt$variance[1, 3, ]
covariance[1, 3, ] <- mnt$variance[2, 3, ]
covariance[1, 4, ] <- mnt$variance[1, 4, ]
covariance[1, 5, ] <- mnt$variance[2, 4, ]
covariance[1, 6, ] <- mnt$variance[3, 4, ]
system.time(mnt <- cmpapproximateMoment_LogisticNorm(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx[1], mnt_num))
mean[2, , ] <- mnt$mean
variance[2, 1, ] <- mnt$variance[1, 1, ]
variance[2, 2, ] <- mnt$variance[2, 2, ]
variance[2, 3, ] <- mnt$variance[3, 3, ]
variance[2, 4, ] <- mnt$variance[4, 4, ]
covariance[2, 1, ] <- mnt$variance[1, 2, ]
covariance[2, 2, ] <- mnt$variance[1, 3, ]
covariance[2, 3, ] <- mnt$variance[2, 3, ]
covariance[2, 4, ] <- mnt$variance[1, 4, ]
covariance[2, 5, ] <- mnt$variance[2, 4, ]
covariance[2, 6, ] <- mnt$variance[3, 4, ]

save(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num, mean, variance, covariance,
     file = "./Output/Output v1.0/Test v1.0/TEST_LogisticNorm_Moment.rda")

load("./Output/Output v1.0/Test v1.0/TEST_LogisticNorm_Moment.rda")

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_LogisticNorm_Moment.pdf", width = 8, height = 18)
par(mfrow = c(3, 1), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
k <- (int_gen + 1):lst_gen
plot(0, type = "n", xlim = c(int_gen + 1, lst_gen), ylim = c(min(mean), max(mean)),
     xlab = "Generation", ylab = "Mean",
     main = "Moment approximation of the logistic normal approximation: mean")
for (i in 1:2) {
  for (j in 1:4) {
    lines(k, mean[i, j, ], type = "l", lty = i, lwd = 1.5, col = brewer.pal(n = 4, name = 'Set1')[j])
  }
}
legend("topleft", legend = c("Wright-Fisher", "logistic normal"), col = "black", lty = 1:2, lwd = 1.5, bty = "n", cex = 1.5)
legend("bottomleft", legend = c("A1B1", "A1B2", "A2B1", "A2B2"), col = brewer.pal(n = 4, name = 'Set1'), lty = 1, lwd = 1.5, bty = "n", cex = 1.5)
plot(0, type = "n", xlim = c(int_gen + 1, lst_gen), ylim = c(min(variance), max(variance)),
     xlab = "Generation", ylab = "Variance",
     main = "Moment approximation of the logistic normal approximation: variance")
for (i in 1:2) {
  for (j in 1:4) {
    lines(k, variance[i, j, ], type = "l", lty = i, lwd = 1.5, col = brewer.pal(n = 4, name = 'Set1')[j])
  }
}
legend("topleft", legend = c("Wright-Fisher", "logistic normal"), col = "black", lty = 1:2, lwd = 1.5, bty = "n", cex = 1.5)
legend("bottomleft", legend = c("A1B1", "A1B2", "A2B1", "A2B2"), col = brewer.pal(n = 4, name = 'Set1'), lty = 1, lwd = 1.5, bty = "n", cex = 1.5)
plot(0, type = "n", xlim = c(int_gen + 1, lst_gen), ylim = c(min(covariance), max(covariance)),
     xlab = "Generation", ylab = "Covariance",
     main = "Moment approximation of the logistic normal approximation: covariance")
for (i in 1:2) {
  for (j in 1:6) {
    lines(k, covariance[i, j, ], type = "l", lty = i, lwd = 1.5, col = brewer.pal(n = 6, name = 'Set1')[j])
  }
}
legend("topleft", legend = c("Wright-Fisher", "logistic normal"), col = "black", lty = 1:2, lwd = 1.5, bty = "n", cex = 1.5)
legend("bottomleft", legend = c("A1B1/A1B2", "A1B1/A2B1", "A1B2/A2B1", "A1B1/A2B2", "A1B2/A2B2", "A2B1/A2B2"), col = brewer.pal(n = 6, name = 'Set1'), lty = 1, lwd = 1.5, bty = "n", cex = 1.5)
dev.off()

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
rec_rat <- 1e-05
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
rec_rat <- 1e-05
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
sim_num <- 1e+05
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05

system.time(smp_apx <- cmpgenerateSample_HierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx[1], mnt_num))

smp_mod <- array(NA, dim = c(4, lst_gen - int_gen, sim_num))
for (i in 1:sim_num) {
  print(i)
  smp_mod[, , i] <- cmpsimulateWFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)[, 2:(lst_gen - int_gen + 1)]
}

save(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, mnt_apx, mnt_num, smp_mod, smp_apx,
     file = "./Output/Output v1.0/Test v1.0/TEST_HierarchicalBeta_Approx.rda")

load("./Output/Output v1.0/Test v1.0/TEST_HierarchicalBeta_Approx.rda")

gen <- 100
smp_mod <- smp_mod[, gen, ]
smp_apx <- smp_apx[, gen, ]

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_HierarchicalBeta_Approx.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist_mod <- hist(smp_mod[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), plot = FALSE)
hist(smp_mod[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A1B1 haplotype at generation", gen))
hist(smp_apx[1, ], breaks = seq(min(smp_mod[1, ], smp_apx[1, ]), max(smp_mod[1, ], smp_apx[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist_mod <- hist(smp_mod[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), plot = FALSE)
hist(smp_mod[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A1B2 haplotype at generation", gen))
hist(smp_apx[2, ], breaks = seq(min(smp_mod[2, ], smp_apx[2, ]), max(smp_mod[2, ], smp_apx[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist_mod <- hist(smp_mod[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), plot = FALSE)
hist(smp_mod[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A2B1 haplotype at generation", gen))
hist(smp_apx[3, ], breaks = seq(min(smp_mod[3, ], smp_apx[3, ]), max(smp_mod[3, ], smp_apx[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist_mod <- hist(smp_mod[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), plot = FALSE)
hist_apx <- hist(smp_apx[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), plot = FALSE)
hist(smp_mod[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ])), ylim = c(0, max(hist_mod$density, hist_apx$density)),
     xlab = "Haplotype frequency", main = paste("Histogram of the A2B2 haplotype at generation", gen))
hist(smp_apx[4, ], breaks = seq(min(smp_mod[4, ], smp_apx[4, ]), max(smp_mod[4, ], smp_apx[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
dev.off()

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
rec_rat <- 1e-05
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
mnt_apx <- c("MC", "Lacerda", "Terhorst", "Paris1", "Paris2")
mnt_num <- 1e+05

mean <- array(NA, dim = c(2, 4, lst_gen - int_gen))
variance <- array(NA, dim = c(2, 4, lst_gen - int_gen))
covariance <- array(NA, dim = c(2, 6, lst_gen - int_gen))
system.time(mnt <- cmpapproximateMoment_WFM(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx[1], mnt_num))
mean[1, , ] <- mnt$mean
variance[1, 1, ] <- mnt$variance[1, 1, ]
variance[1, 2, ] <- mnt$variance[2, 2, ]
variance[1, 3, ] <- mnt$variance[3, 3, ]
variance[1, 4, ] <- mnt$variance[4, 4, ]
covariance[1, 1, ] <- mnt$variance[1, 2, ]
covariance[1, 2, ] <- mnt$variance[1, 3, ]
covariance[1, 3, ] <- mnt$variance[2, 3, ]
covariance[1, 4, ] <- mnt$variance[1, 4, ]
covariance[1, 5, ] <- mnt$variance[2, 4, ]
covariance[1, 6, ] <- mnt$variance[3, 4, ]
system.time(mnt <- cmpapproximateMoment_HierarchicalBeta(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx[1], mnt_num))
mean[2, , ] <- mnt$mean
variance[2, 1, ] <- mnt$variance[1, 1, ]
variance[2, 2, ] <- mnt$variance[2, 2, ]
variance[2, 3, ] <- mnt$variance[3, 3, ]
variance[2, 4, ] <- mnt$variance[4, 4, ]
covariance[2, 1, ] <- mnt$variance[1, 2, ]
covariance[2, 2, ] <- mnt$variance[1, 3, ]
covariance[2, 3, ] <- mnt$variance[2, 3, ]
covariance[2, 4, ] <- mnt$variance[1, 4, ]
covariance[2, 5, ] <- mnt$variance[2, 4, ]
covariance[2, 6, ] <- mnt$variance[3, 4, ]

save(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, mnt_apx, mnt_num, mean, variance, covariance,
     file = "./Output/Output v1.0/Test v1.0/TEST_HierarchicalBeta_Moment.rda")

load("./Output/Output v1.0/Test v1.0/TEST_HierarchicalBeta_Moment.rda")

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_HierarchicalBeta_Moment.pdf", width = 8, height = 18)
par(mfrow = c(3, 1), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
k <- (int_gen + 1):lst_gen
plot(0, type = "n", xlim = c(int_gen + 1, lst_gen), ylim = c(min(mean), max(mean)),
     xlab = "Generation", ylab = "Mean",
     main = "Moment approximation of the hierarchical beta approximation: mean")
for (i in 1:2) {
  for (j in 1:4) {
    lines(k, mean[i, j, ], type = "l", lty = i, lwd = 1.5, col = brewer.pal(n = 4, name = 'Set1')[j])
  }
}
legend("topleft", legend = c("Wright-Fisher", "hierarchical beta"), col = "black", lty = 1:2, lwd = 1.5, bty = "n", cex = 1.5)
legend("bottomleft", legend = c("A1B1", "A1B2", "A2B1", "A2B2"), col = brewer.pal(n = 4, name = 'Set1'), lty = 1, lwd = 1.5, bty = "n", cex = 1.5)
plot(0, type = "n", xlim = c(int_gen + 1, lst_gen), ylim = c(min(variance), max(variance)),
     xlab = "Generation", ylab = "Variance",
     main = "Moment approximation of the hierarchical beta approximation: variance")
for (i in 1:2) {
  for (j in 1:4) {
    lines(k, variance[i, j, ], type = "l", lty = i, lwd = 1.5, col = brewer.pal(n = 4, name = 'Set1')[j])
  }
}
legend("topleft", legend = c("Wright-Fisher", "hierarchical beta"), col = "black", lty = 1:2, lwd = 1.5, bty = "n", cex = 1.5)
legend("bottomleft", legend = c("A1B1", "A1B2", "A2B1", "A2B2"), col = brewer.pal(n = 4, name = 'Set1'), lty = 1, lwd = 1.5, bty = "n", cex = 1.5)
plot(0, type = "n", xlim = c(int_gen + 1, lst_gen), ylim = c(min(covariance), max(covariance)),
     xlab = "Generation", ylab = "Covariance",
     main = "Moment approximation of the hierarchical beta approximation: covariance")
for (i in 1:2) {
  for (j in 1:6) {
    lines(k, covariance[i, j, ], type = "l", lty = i, lwd = 1.5, col = brewer.pal(n = 6, name = 'Set1')[j])
  }
}
legend("topleft", legend = c("Wright-Fisher", "hierarchical beta"), col = "black", lty = 1:2, lwd = 1.5, bty = "n", cex = 1.5)
legend("bottomleft", legend = c("A1B1/A1B2", "A1B1/A2B1", "A1B2/A2B1", "A1B1/A2B2", "A1B2/A2B2", "A2B1/A2B2"), col = brewer.pal(n = 6, name = 'Set1'), lty = 1, lwd = 1.5, bty = "n", cex = 1.5)
dev.off()

################################################################################

#' Calculate the root mean square deviation between the two empirical cumulative distribution functions
#' Parameter setting
#' @param smp_mod the sample trajectories generated under the Wright-Fisher model
#' @param smp_apx the sample trajectories generated under the approximation of the Wright-Fisher model
#' @param grd_num the grid number for the empirical probability distribution function
#' @param rnd_grd = TRUE/FALSE (return the root mean square deviation between the empirical cumulative distribution functions with random grid points or not)





################################################################################
