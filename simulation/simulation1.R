## Simulation study comparing GAM2d (implemented in the "gam" function)

################################################################################
## Load packages
################################################################################
rm(list = ls())
library(refund)
library(dplyr)
library(mgcv)
library(ggplot2)
library(ggpubr)

source("simu_generate.R") 
source("2DFMM.R") 

################################################################################
## Set parameters
################################################################################
n <- 50 ## number of subjects
S <- 10 ##number of observations per subject
T <- 100 ## dimension of the functional domain
type <- 1 ##beta_true type 
SNR_B <- 1 ## relative importance of random effects
SNR_sigma <- 1 ## signal-to-noise ratio
family <- "gaussian" ## family of longitudinal functional data

grid <- seq(0,1,length.out=T)
## self-defined phi
phi <- array(0, dim=c(2, T))
phi[1,] <- (1.5 - sin(2*grid*pi) - cos(2*grid*pi) )
phi[1,] <- phi[1,] / sqrt(sum(phi[1,]^2))
phi[2,] <- sin(4*grid*pi)
phi[2,] <- phi[2,] / sqrt(sum(phi[2,]^2))

## self-defined phi1
phi1 <- array(0, dim=c(2, T))
phi1[1,] <- cos(2*pi*grid)
phi1[1,] <- phi1[1,] / sqrt(sum(phi1[1,]^2))
phi1[2,] <- sin(2*pi*grid)
phi1[2,] <- phi1[2,] / sqrt(sum(phi1[2,]^2))

## self-defined phi2
phi2 <- array(0, dim=c(2, T))
phi2[1,] <- sqrt(3)*(2*grid-1)
phi2[1,] <- phi2[1,] / sqrt(sum(phi2[1,]^2))
phi2[2,] <- sqrt(5)*(6*grid^2 - 6*grid + 1)
phi2[2,] <- phi2[2,] / sqrt(sum(phi2[2,]^2))


################################################################################
## Do simulations on a local laptop 
################################################################################
nsim <- 100
sim_local <- list() ## store simulation results

for(iter in 1:nsim){
  
  ################################################################################
  ## Generate simulated data
  ################################################################################
  data <- sim_generate(family = "gaussian", n = n, S = S, T = T, 
                       type = type, phi = phi, phi1 = phi1, phi2 = NULL, rho = 0.5)   
  beta_true <- array(0, dim = c(T,S,2))
  beta_true[,,1] <- t(intercept(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))) #intercept
  if(type==1){
    beta_true[,,2] <- t(f1(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))) #beta1
  }else{
    beta_true[,,2] <- t(f2(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))) #beta2
  }

  ################################################################################
  ## Implement different estimation methods
  ################################################################################
  ## fit the fmm2d model 
  ## (S1): S-3; (S2): max(round(S/4),4)
  ptm <- proc.time()
  fit_fmm2d <- fmm2d(formula=Y~X, data=data, S=S, smoother="te", knots=c(max(round(S/4),4), min(round(T/4), 35)),  parallel = TRUE)
  time_fmm2d <- (proc.time() - ptm)[3]
  
  ## fit gam2d model
  T.arg <- rep(seq(0,1, length.out = T), n*S)
  S.arg <- rep(rep(seq(0,1, length.out = S), each=T), n)
  gam.data <- data.frame(T.arg=T.arg, S.arg=S.arg, X=rep(data$X, each=T), Y=as.vector(t(data$Y)))
  
  ptm <- proc.time()
  fit_gam <-  mgcv::gam(Y ~ te(T.arg, S.arg, k = c(min(round(T/4), 35), max(round(S/4),4))) + te(T.arg, S.arg, by=X, k=c(min(round(T/4), 35), max(round(S/4),4))), method="REML", data = gam.data)
  time_gam2d <- (proc.time() - ptm)[3]
  
  library(gratia)
  new.data <- data.frame(T.arg=rep(seq(0,1,length.out=T), S), S.arg=rep(seq(0,1,length.out=S),each=T), X=1)
  fit_gam2d <- list(intercept = smooth_estimates(fit_gam, smooth="te(T.arg,S.arg)", data=new.data), 
                    X= smooth_estimates(fit_gam, smooth="te(T.arg,S.arg):X", data=new.data))
  #b <- plot(fit.gam,scheme=2)
  #draw(fit.gam)

  ################################################################################
  ## Organize simulation results
  ################################################################################
  
  ## fmm2d
  MISE_fmm2d <- rep(0,2); names(MISE_fmm2d) <- c("intercept", "beta")
  for(p in 1:2) MISE_fmm2d[p] <- mean((fit_fmm2d$betaHat[[p]] - beta_true[,,p])^2)
  print(MISE_fmm2d)
  
  cover_pw_fmm2d <- rep(0,2); names(cover_pw_fmm2d) <- c("intercept", "beta")
  MIAW_fmm2d <- rep(0,2); names(MIAW_fmm2d) <- c("intercept", "beta")
  for(p in 1:2){
    betaHat.var <- fit_fmm2d$betaHat.cov[,,p]
    betaHat.up <- as.vector(fit_fmm2d$betaHat[[p]]) + 2*sqrt(diag(betaHat.var))
    betaHat.low <- as.vector(fit_fmm2d$betaHat[[p]]) - 2*sqrt(diag(betaHat.var))
    cover_upper <- which(betaHat.up > as.vector(beta_true[,,p]))
    cover_lower <- which(betaHat.low < as.vector(beta_true[,,p]))
    cover_pw_fmm2d[p] <- length(intersect(cover_upper, cover_lower))/(S*T)
    MIAW_fmm2d[p] <- sum(betaHat.up - betaHat.low)/(S*T)
  }
  
  ## gam2d
  MISE_gam2d <- rep(0,2); names(MISE_gam2d) <- c("intercept", "beta")
  for(p in 1:2) MISE_gam2d[p] <- colMeans((fit_gam2d[[p]][,"est"] - as.vector(beta_true[,,p]))^2)
  print(MISE_gam2d)
  
  cover_pw_gam2d <- rep(0,2); names(cover_pw_gam2d) <- c("intercept", "beta")
  MIAW_gam2d <- rep(0,2); names(MIAW_gam2d) <- c("intercept", "beta")
  for(p in 1:2){
    betaHat.up <- fit_gam2d[[p]][,"est"] + 2*fit_gam2d[[p]][,"se"]
    betaHat.low <- fit_gam2d[[p]][,"est"] - 2*fit_gam2d[[p]][,"se"]
    cover_upper <- which(betaHat.up > as.vector(beta_true[,,p]))
    cover_lower <- which(betaHat.low < as.vector(beta_true[,,p]))
    cover_pw_gam2d[p] <- length(intersect(cover_upper, cover_lower))/(S*T)
    MIAW_gam2d[p] <-sum(betaHat.up - betaHat.low)/(S*T)
  }
  
  sim_result <- list(MISE_fmm2d = MISE_fmm2d, MISE_gam2d = MISE_gam2d,
                     cover_pw_fmm2d = cover_pw_fmm2d, cover_pw_gam2d = cover_pw_gam2d,
                     MIAW_fmm2d = MIAW_fmm2d, MIAW_gam2d = MIAW_gam2d,
                     time_fmm2d = time_fmm2d, time_gam2d = time_gam2d)

  sim_local[[iter]] <- sim_result
  print(iter)
}


MISE_fmm2d <- lapply(sim_local, '[[', 1) %>% bind_rows()
MISE_gam2d <- lapply(sim_local, '[[', 2) %>% bind_rows()

cover_pw_fmm2d <- lapply(sim_local, '[[', 3) %>% bind_rows()
cover_pw_gam2d <- lapply(sim_local, '[[', 4) %>% bind_rows()

MIAW_fmm2d <- lapply(sim_local, '[[', 5) %>% bind_rows()
MIAW_gam2d <- lapply(sim_local, '[[', 6) %>% bind_rows()

time_fmm2d <- lapply(sim_local, '[[', 7) %>% bind_rows()
time_gam2d <- lapply(sim_local, '[[', 8) %>% bind_rows()

colMeans(MISE_fmm2d); colMeans(MISE_gam2d)
colMeans(cover_pw_fmm2d); colMeans(cover_pw_gam2d)
colMeans(time_fmm2d); colMeans(time_gam2d)























