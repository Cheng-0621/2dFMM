## Simulation study comparing FUI (implemented in the "lfosr3s" function), 
## FAMM (implemented in the "pffr" function) and FILF (implemented in "fgee" function)

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
source("lfosr3s.R") ## required packages: lme4, refund, dplyr, mgcv, progress, mvtnorm, parallel
source("fgee.R") 

################################################################################
## Do simulations on a local laptop 
################################################################################
n <- 50 ## number of subjects
S <- 10 ##number of observations per subject
T <- 100 ## dimension of the functional domain
type <- 2 ##beta_true type 
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
## Bivariate persepctive and univariate perspective
################################################################################
for (bi in c("T", "F")){
  for (rho in c(0.5,2,6)){
    
    nsim <- 100
    sim_local <- list() ## store simulation results
    
    for(iter in 1:nsim){
      
      ################################################################################
      ## Generate simulated data
      ################################################################################
      data <- sim_generate(family = "gaussian", bivariate = bi, n = n, S = S, T = T, 
                           type = type, phi = phi, phi1 = phi1, phi2 = NULL, rho = rho)   
      
      beta_true <- array(0, dim = c(T,S,2))
      beta_true[,,1] <- t(intercept(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))) #intercept
      
      if(type==1){
        beta_true[,,2] <- t(f1(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))) #beta1
      }else{
        beta_true[,,2] <- t(f2(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))) #beta2
      }
      
      beta_true_t <- array(0, dim=c(2, T))
      for(p in 1:2) beta_true_t[p,] <- apply(beta_true[,,p], 1, mean)
      
      ################################################################################
      ## Implement different estimation methods
      ################################################################################
      ## fit the fmm2d model 
      ptm <- proc.time()
      fit_fmm2d <- fmm2d(formula=Y ~ X , data=data, S=S, smoother="te", knots=c(max(round(S/4),4), min(round(T/4, 35))),  parallel = FALSE, scb = FALSE)
      time_fmm2d <- (proc.time() - ptm)[3]
      
      ## fit the lfosr3s model
      ptm <- proc.time()
      fit_lfosr3s <- lfosr3s(formula = Y ~ X + Longit + (1 + Longit | ID), data = data, family = "gaussian", 
                             var = TRUE, analytic = TRUE, scb = FALSE, parallel = FALSE, silent = FALSE)
      time_lfosr3s <- (proc.time() - ptm)[3]
      
      ## fit the pffr model
      # we use "bam" for faster computation, and change k = 15 in bs.yindex to capture curvatures in both scenarios
      #+ s(Longit, ID, bs = "re")
      ptm <- proc.time()
      fit_pffr <- pffr(formula = Y ~ X + Longit + s(ID, bs = "re") + s(Longit, ID, bs = "re"), data = data, family = "gaussian", algorithm = "bam",
                       bs.yindex = list(bs = "ps", k = 15 , m = c(2, 1)))
      time_pffr <- (proc.time() - ptm)[3]
      
      ## fit the filf model
      T.arg <- seq(0,1, length.out = T)
      S.arg <- rep(seq(0,1, length.out = S), n)
      
      filf.x <- matrix(data$X); colnames(filf.x) <- "X"
      filf.data <- list(Y=matrix(data$Y, nrow = S*n), X=filf.x, ID=data$ID)
      
      ptm <- proc.time()
      fit_filf <-  fgee(formula = Y ~ X, Y= filf.data$Y, Cov=filf.data$X,
                        s=T.arg, subjID=filf.data$ID, 
                        Tij=S.arg, numLongiPoints=S,
                        control=list(FPCA1.pve=0.95, FPCA2.pve=0.95, FPCA1.knots= min(round(T/4, 35)), FPCA2.knots=max(round(S/4),4)),
                        corstr = "exchangeable")
      time_filf <- (proc.time() - ptm)[3]
      
      ################################################################################
      ## Organize simulation results
      ################################################################################
      ## fmm2d
      #betaHat.s <- array(0, dim=c(2, S))
      betaHat.t <- array(0, dim=c(2, T))     
      #for(p in 1:2) betaHat.s[p,] <- apply(fit_fmm2d$betaHat[[p]], 2, mean)
      for(p in 1:2) betaHat.t[p,] <- apply(fit_fmm2d$betaHat[[p]], 1, mean)
      
      #for(p in 1:2){
      #  for(t in 1:T){
      #    betaHat.t[p,t] <- mean(lm(betaHat[[p]][t,]~betaHat.s[p,])$fitted.values)
      #  }
      #}
      
      betaHat.t.cov <- array(0, dim=c(T, T, 2))
      for(p in 1:2){
        for(t in 1:T){
          for(v in t:T){
            cov.tmp <- 0
            for(s in 1:S){
              for(u in 1:S){
                cov.tmp <-  fit_fmm2d$betaHat.cov[(s-1)*T+t, (u-1)*T+v, p] + cov.tmp
              }
            }
            betaHat.t.cov[v,t,p] <- betaHat.t.cov[t,v,p] <- cov.tmp/(S^2)
          }
        }
      }
      
      MISE_fmm2d <- rep(0,2); names(MISE_fmm2d) <- c("intercept", "beta")
      for(p in 1:2) MISE_fmm2d[p] <- mean((betaHat.t[p,] - beta_true_t[p,])^2)
      print(MISE_fmm2d)
      
      cover_pw_fmm2d <- rep(0,2); names(cover_pw_fmm2d) <- c("intercept", "beta")
      MIAW_fmm2d <- rep(0,2); names(MIAW_fmm2d) <- c("intercept", "beta")
      for(p in 1:2){
        betaHat.var <- betaHat.t.cov[,,p]
        betaHat.up <- betaHat.t[p,] + 1.96*sqrt(diag(betaHat.var))
        betaHat.low <- betaHat.t[p,] - 1.96*sqrt(diag(betaHat.var))
        cover_upper <- which(betaHat.up > beta_true_t[p,])
        cover_lower <- which(betaHat.low < beta_true_t[p,])
        cover_pw_fmm2d[p] <- length(intersect(cover_upper, cover_lower))/T
        MIAW_fmm2d[p] <- sum(betaHat.up - betaHat.low)/T
      }
      print(paste(cover_pw_fmm2d))
      
      ## FUI
      MISE_lfosr3s <- rep(0,2); names(MISE_lfosr3s) <- c("intercept", "beta")
      for(p in 1:2) MISE_lfosr3s[p] <- mean((fit_lfosr3s$betaHat[2,] - beta_true_t[2,])^2)
      print(MISE_lfosr3s)
      
      cover_lfosr3s <- rep(0,2); names(cover_lfosr3s) <- c("intercept", "beta")
      MIAW_lfosr3s <- rep(0,2); names(MIAW_lfosr3s) <- c("intercept", "beta")
      for(p in 1:2){
        betaHat.var <- fit_lfosr3s$betaHat.var[,,p]
        betaHat.up <- fit_lfosr3s$betaHat[p,] + 1.96*sqrt(diag(fit_lfosr3s$betaHat.var[,,p]))
        betaHat.low <- fit_lfosr3s$betaHat[p,] - 1.96*sqrt(diag(fit_lfosr3s$betaHat.var[,,p]))
        cover_upper <- which(betaHat.up > beta_true_t[p,])
        cover_lower <- which(betaHat.low < beta_true_t[p,])
        cover_lfosr3s[p] <- length(intersect(cover_upper, cover_lower))/T
        MIAW_lfosr3s[p] <- sum(betaHat.up - betaHat.low)/T
      }
      print(paste(cover_lfosr3s))
      
      ## FAMM
      coef_pffr <- coef(fit_pffr, n1 = T)
      betaHat_pffr <- betaHat_pffr.var <- array(0, dim=c(2,T))
      betaHat_pffr[1,] <- as.vector(coef_pffr$smterms$`Intercept(yindex)`$value) + coef_pffr$pterms[1]
      betaHat_pffr[2,] <- as.vector(coef_pffr$smterms$`X(yindex)`$value)
      betaHat_pffr.var[1,] <- as.vector(coef_pffr$smterms$`Intercept(yindex)`$se^2)
      betaHat_pffr.var[2,] <- as.vector(coef_pffr$smterms$`X(yindex)`$se^2)
      
      MISE_pffr <- rep(0,2); names(MISE_pffr) <- c("intercept", "beta")
      for(p in 1:2) MISE_pffr[p] <- mean((betaHat_pffr[p,] - beta_true_t[p,])^2)
      print(MISE_pffr)
      
      cover_pffr <- rep(0,2); names(cover_pffr) <- c("intercept", "beta")
      MIAW_pffr <- rep(0,2); names(MIAW_pffr) <- c("intercept", "beta")
      for(p in 1:2){
        betaHat.up <- betaHat_pffr[p,] + 1.96*sqrt(betaHat_pffr.var[p,])
        betaHat.low <- betaHat_pffr[p,] - 1.96*sqrt(betaHat_pffr.var[p,])
        cover_upper <- which(betaHat.up > beta_true_t[p,])
        cover_lower <- which(betaHat.low < beta_true_t[p,])
        cover_pffr[p] <- length(intersect(cover_upper, cover_lower))/T
        MIAW_pffr[p] <- sum(betaHat.up - betaHat.low)/T
      }
      
      print(paste(cover_pffr))
      
      ## FILF
      MISE_filf <- rep(0,2); names(MISE_filf) <- c("intercept", "beta")
      for(p in 1:2) MISE_filf[p] <- mean((fit_filf$beta[,p] - beta_true_t[p,])^2)
      print(MISE_filf)
      
      cover_filf <- rep(0,2); names(cover_filf) <- c("intercept", "beta")
      MIAW_filf <- rep(0,2); names(MIAW_filf) <- c("intercept", "beta")
      for(p in 1:2){
        betaHat.up <- fit_filf$beta[,p] + 1.96*fit_filf$beta.se[,p]
        betaHat.low <- fit_filf$beta[,p] - 1.96*fit_filf$beta.se[,p]
        cover_upper <- which(betaHat.up > beta_true_t[p,])
        cover_lower <- which(betaHat.low < beta_true_t[p,])
        cover_filf[p] <- length(intersect(cover_upper, cover_lower))/T
        MIAW_filf[p] <- sum(betaHat.up - betaHat.low)/T
      }
      print(paste(cover_filf))
      
      sim_result <- list(MISE_fmm2d = MISE_fmm2d, MISE_lfosr3s = MISE_lfosr3s, MISE_pffr = MISE_pffr, MISE_filf = MISE_filf,
                         cover_pw_fmm2d = cover_pw_fmm2d, cover_lfosr3s = cover_lfosr3s, cover_pffr = cover_pffr, cover_filf = cover_filf,
                         MIAW_fmm2d =  MIAW_fmm2d, MIAW_lfosr3s = MIAW_lfosr3s, MIAW_pffr = MIAW_pffr, MIAW_filf =  MIAW_filf,
                         time_fmm2d = time_fmm2d, time_lfosr3s = time_lfosr3s, time_pffr = time_pffr, time_filf = time_filf)
      
      sim_local[[iter]] <- sim_result
      print(iter)
    }
    
    
    MISE_fmm2d <- lapply(sim_local, '[[', 1) %>% bind_rows()
    MISE_lfosr3s <- lapply(sim_local, '[[', 2) %>% bind_rows()
    MISE_pffr <- lapply(sim_local, '[[', 3) %>% bind_rows()
    MISE_filf <- lapply(sim_local, '[[', 4) %>% bind_rows()
    
    cover_pw_fmm2d <- lapply(sim_local, '[[', 5) %>% bind_rows()
    cover_lfosr3s <- lapply(sim_local, '[[', 6) %>% bind_rows()
    cover_pffr <- lapply(sim_local, '[[', 7) %>% bind_rows()
    cover_filf <- lapply(sim_local, '[[', 8) %>% bind_rows()
    
    MIAW_fmm2d <- lapply(sim_local, '[[', 9) %>% bind_rows()
    MIAW_lfosr3s <- lapply(sim_local, '[[', 10) %>% bind_rows()
    MIAW_pffr <- lapply(sim_local, '[[', 11) %>% bind_rows()
    MIAW_filf <- lapply(sim_local, '[[', 12) %>% bind_rows()
    
    time_fmm2d <- lapply(sim_local, '[[', 13) %>% bind_rows()
    time_lfosr3s <- lapply(sim_local, '[[', 14) %>% bind_rows()
    time_pffr <- lapply(sim_local, '[[', 15) %>% bind_rows()
    time_filf <- lapply(sim_local, '[[', 16) %>% bind_rows()
    
    
    colMeans(MISE_fmm2d); colMeans(MISE_lfosr3s); colMeans(MISE_pffr); colMeans(MISE_filf)
    colMeans(cover_pw_fmm2d); colMeans(cover_lfosr3s); colMeans(cover_pffr); colMeans(cover_filf)
    
    colMeans(MIAW_fmm2d); colMeans(MIAW_lfosr3s); colMeans(MIAW_pffr); colMeans(MIAW_filf)
    colMeans(time_fmm2d); colMeans(time_lfosr3s);  colMeans(time_pffr); colMeans(time_filf)
    
  }
}


