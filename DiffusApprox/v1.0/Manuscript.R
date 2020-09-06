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

# install.packages("pals")
library("pals")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

#install.packages("MASS")
library("MASS")

#install.packages("coda")
library("coda")

#install.packages("Metrics")
library("Metrics")

setwd("~/Dropbox/Jeffery He/iResearch/Publications/2012/HE2020-WFM-2L-DiffusApprox-TheorPopulBiol")

########################################################################################################################

#' @section Section 1. Introduction

########################################################################################################################

#' @section Section 2. Wright-Fisher diffusion

#' @section Section 2.1. Wright-Fisher diffusion

#' @section Section 2.2. Diffusion approximation

########################################################################################################################

#' @section Section 3. Stochastic Taylor scheme

#' @section Section 3.1. Equivalent SDE formulation

#' @section Section 3.2. Performance evaluation

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(2e-01, 2e-01, 2e-01, 4e-01)
int_gen <- 0
lst_gen <- 100
sim_num <- 1e+05

dist <- array(NA, dim = c(3, 4, lst_gen - int_gen + 1))
model <- "WFM"
smp_WFM <- cmpgenerateSampleTraj_2L(model, sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num)
model <- "WFD"
ptn_num <- 1e+00
smp_WFD <- cmpgenerateSampleTraj_2L(model, sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, ptn_num)
grd_num <- 2 * pop_siz
dist[1, , ] <- cmpcalculateHellingerDist(smp_WFM, smp_WFD, grd_num)
model <- "WFD"
ptn_num <- 1e+01
smp_WFD <- cmpgenerateSampleTraj_2L(model, sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, ptn_num)
grd_num <- 2 * pop_siz
dist[2, , ] <- cmpcalculateHellingerDist(smp_WFM, smp_WFD, grd_num)
model <- "WFD"
ptn_num <- 1e+02
smp_WFD <- cmpgenerateSampleTraj_2L(model, sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, ptn_num)
grd_num <- 2 * pop_siz
dist[3, , ] <- cmpcalculateHellingerDist(smp_WFM, smp_WFD, grd_num)

plot(int_gen:lst_gen, dist[1, 4, ], type = "l")
lines(int_gen:lst_gen, dist[2, 4, ])
lines(int_gen:lst_gen, dist[3, 4, ])

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(5e-03, 5e-03, 9e-02, 9e-01)
int_gen <- 0
lst_gen <- 500
sim_num <- 5e+05
ptn_num <- c(1e+00, 1e+01, 1e+02)

dist <- array(NA, dim = c(3, 4, lst_gen - int_gen + 1))
for (i in 1:3) {
  model <- "WFM"
  smp_WFM <- cmpgenerateSampleTraj_2L(model, sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num)
  model <- "WFD"
  smp_WFD <- cmpgenerateSampleTraj_2L(model, sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, sim_num, ptn_num[i])
  grd_num <- 2 * pop_siz

  dist[i, , ] <- cmpcalculateHellingerDist(smp_WFM, smp_WFD, grd_num)
}


########################################################################################################################

#' @section Section 4. Transition probability density approxiamtion

#' @section Section 3.1. Euler-Maruyama approximation

#' @section Section 3.2. Monte Carlo integration

#' @section Section 3.3. Monte Carlo integration with importance sampling

#' @section Section 3.4. Performance evaluation

########################################################################################################################

#' @section Section 3. Discussion

########################################################################################################################

