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

setwd("~/Documents/CityU/research/Project on Shanghai Actigraph Data")
source("codes/simu_generate.R") 
source("codes/2DFMM.R")
source("reference/codes/FUI/lfosr3s.R") ## required packages: lme4, refund, dplyr, mgcv, progress, mvtnorm, parallel
source("reference/codes/FILF/fgee.R") 

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
      fit_fmm2d <- fmm2d(formula=Y ~ X , data=data, S=S, smoother="te", knots=c(max(round(S/4),4), min(round(T/4, 35))),
                         fpca.opt = list(dataType = 'Dense', methodSelectK = 'FVE'),  parallel = FALSE)
      time_fmm2d <- (proc.time() - ptm)[3]
      
      ## fit the lfosr3s model
      ptm <- proc.time()
      fit_lfosr3s <- lfosr3s(formula = Y ~ X + Longit + (1 + Longit | ID), data = data, family = "gaussian", 
                             var = TRUE, analytic = TRUE, parallel = FALSE, silent = FALSE)
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
    
    
    save(MISE_fmm2d, MISE_lfosr3s, MISE_pffr, MISE_filf,
         cover_pw_fmm2d, cover_lfosr3s, cover_pffr, cover_filf,
         MIAW_fmm2d, MIAW_lfosr3s, MIAW_pffr, MIAW_filf,
         time_fmm2d, time_lfosr3s, time_pffr, time_filf,
         file = ifelse(bi == "T", 
                       paste0("results/simu2S2biT_5010100rho", rho, ".RData"),
                       paste0("results/simu2S2uniT_5010100rho", rho, ".RData")))
  }
}



#load results to draw MISE comparison
load_results <- function(arg, val, dir, iftime = FALSE){
  
  load(paste0("codes/simu_results/", dir[1], ".RData"))
  MISEf1 <- cbind(MISE_fmm2d, method="FMM2d", arg = val[1])
  MISEl1 <- cbind(MISE_lfosr3s, method="FUI", arg = val[1])
  MISEp1 <- cbind(MISE_pffr, method="FAMM", arg = val[1])
  MISEfi1 <- cbind(MISE_filf, method="FILF", arg = val[1])
  
  coverf1 <- cbind(cover_pw_fmm2d, method="FMM2d", arg = val[1])
  coverl1 <- cbind(cover_lfosr3s, method="FUI", arg = val[1])
  coverp1 <- cbind(cover_pffr, method="FAMM", arg = val[1])
  coverfi1 <- cbind(cover_filf, method="FILF", arg = val[1])
  
  MIAWf1 <- cbind(MIAW_fmm2d, method="FMM2d", arg = val[1])
  MIAWl1 <- cbind(MIAW_lfosr3s, method="FUI", arg = val[1])
  MIAWp1 <- cbind(MIAW_pffr, method="FAMM", arg = val[1])
  MIAWfi1 <- cbind(MIAW_filf, method="FILF", arg = val[1])
  
  timef1 <- cbind(time_fmm2d, method="FMM2d", arg = val[1])
  timel1 <- cbind(time_lfosr3s, method="FUI", arg = val[1])
  timep1 <- cbind(time_pffr, method="FAMM", arg = val[1])
  timefi1 <- cbind(time_filf, method="FILF", arg = val[1])
  
  load(paste0("codes/simu_results/", dir[2], ".RData"))
  MISEf2 <- cbind(MISE_fmm2d, method="FMM2d", arg = val[2])
  MISEl2 <- cbind(MISE_lfosr3s, method="FUI", arg = val[2])
  MISEp2 <- cbind(MISE_pffr, method="FAMM", arg = val[2])
  MISEfi2 <- cbind(MISE_filf, method="FILF", arg = val[2])
  
  coverf2 <- cbind(cover_pw_fmm2d, method="FMM2d", arg = val[2])
  coverl2 <- cbind(cover_lfosr3s, method="FUI", arg = val[2])
  coverp2 <- cbind(cover_pffr, method="FAMM", arg = val[2])
  coverfi2 <- cbind(cover_filf, method="FILF", arg = val[2])
  
  MIAWf2 <- cbind(MIAW_fmm2d, method="FMM2d", arg = val[2])
  MIAWl2 <- cbind(MIAW_lfosr3s, method="FUI", arg = val[2])
  MIAWp2 <- cbind(MIAW_pffr, method="FAMM", arg = val[2])
  MIAWfi2 <- cbind(MIAW_filf, method="FILF", arg = val[2])
  
  timef2 <- cbind(time_fmm2d, method="FMM2d", arg = val[2])
  timel2 <- cbind(time_lfosr3s, method="FUI", arg = val[2])
  timep2 <- cbind(time_pffr, method="FAMM", arg = val[2])
  timefi2 <- cbind(time_filf, method="FILF", arg = val[2])
  
  load(paste0("codes/simu_results/", dir[3], ".RData"))
  MISEf3 <- cbind(MISE_fmm2d, method="FMM2d", arg = val[3])
  MISEl3 <- cbind(MISE_lfosr3s, method="FUI", arg = val[3])
  MISEp3 <- cbind(MISE_pffr, method="FAMM", arg = val[3])
  MISEfi3 <- cbind(MISE_filf, method="FILF", arg = val[3])
  
  coverf3 <- cbind(cover_pw_fmm2d, method="FMM2d", arg = val[3])
  coverl3 <- cbind(cover_lfosr3s, method="FUI", arg = val[3])
  coverp3 <- cbind(cover_pffr, method="FAMM", arg = val[3])
  coverfi3 <- cbind(cover_filf, method="FILF", arg = val[3])
  
  MIAWf3 <- cbind(MIAW_fmm2d, method="FMM2d", arg = val[3])
  MIAWl3 <- cbind(MIAW_lfosr3s, method="FUI", arg = val[3])
  MIAWp3 <- cbind(MIAW_pffr, method="FAMM", arg = val[3])
  MIAWfi3 <- cbind(MIAW_filf, method="FILF", arg = val[3])
  
  timef3 <- cbind(time_fmm2d, method="FMM2d", arg = val[3])
  timel3 <- cbind(time_lfosr3s, method="FUI", arg = val[3])
  timep3 <- cbind(time_pffr, method="FAMM", arg = val[3])
  timefi3 <- cbind(time_filf, method="FILF", arg = val[3])
  
  MISE <- rbind(MISEf1, MISEf2, MISEf3, MISEl1, MISEl2, MISEl3,
                MISEp1, MISEp2, MISEp3, MISEfi1, MISEfi2, MISEfi3)
  MISE$method <- factor(MISE$method, levels = c("FMM2d", "FAMM", "FUI", "FILF"))
  colnames(MISE)[4] <- arg
  
  cover <- rbind(coverf1, coverf2, coverf3, coverl1, coverl2, coverl3, 
                 coverp1, coverp2, coverp3, coverfi1, coverfi2, coverfi3)
  cover$method <- factor(cover$method, levels = c("FMM2d", "FAMM", "FUI","FILF"))
  colnames(cover)[4] <- arg
  
  MIAW <- rbind(MIAWf1, MIAWf2, MIAWf3, MIAWl1, MIAWl2, MIAWl3,
                MIAWp1, MIAWp2, MIAWp3, MIAWfi1, MIAWfi2, MIAWfi3)
  MIAW$method <- factor(MIAW$method, levels = c("FMM2d", "FAMM", "FUI", "FILF"))
  colnames(MIAW)[4] <- arg
  
  if(iftime){
    time <- rbind(timef1, timef2, timef3, timel1, timel2, timel3,
                  timep1, timep2, timep3, timefi1, timefi2, timefi3)
    time$method <- factor(MISE$method, levels = c("FMM2d", "FAMM", "FUI", "FILF"))
    colnames(time)[3] <- arg
    ret <- list(MISE=MISE, cover=cover, MIAW=MIAW, time=time)
  }else{
    ret <- list(MISE=MISE, cover=cover, MIAW=MIAW)
  }
  return(ret)
}

ret <- load_results("rho", val=c(0.5, 2, 6), iftime = TRUE, 
                    dir=c("simu2S2biT_5010100rho0.5", 
                          "simu2S2biT_5010100rho2",
                          "simu2S2biT_5010100rho6") )

p1 <- ggplot(ret$MISE, aes(x=factor(rho), y=beta, fill = method)) + 
  geom_boxplot() + 
  theme_bw() +
  scale_fill_discrete(name = "Method") +
  labs(y = expression(paste("MISE(", hat(beta)[1], ")")), x = expression(rho) , 
       title = expression(paste("MISE(", hat(beta)[1], ") (S2)"))) + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size=13),  axis.text.x = element_text(size=13), 
        legend.title = element_text(size=14), legend.text = element_text(size=13),
        plot.title =element_text(size=16, hjust = 0.5))

require(dplyr)
ret$cover %>% dplyr::group_by(method, rho) %>% dplyr::summarise(mean = mean(beta), sd = sd(beta))
ret$time %>% dplyr::group_by(method, rho) %>% dplyr::summarise(mean = mean(elapsed), sd = sd(elapsed))

ret <- load_results("rho", val=c(0.5, 2, 6), iftime = TRUE, 
                    dir=c("simu2S2uniT_5010100rho0.5", 
                          "simu2S2uniT_5010100rho2",
                          "simu2S2uniT_5010100rho6") )

p2 <- ggplot(ret$MISE, aes(x=factor(rho), y=beta, fill = method)) + 
  geom_boxplot() + 
  theme_bw() +
  scale_fill_discrete(name = "Method") +
  labs(y = expression(paste("MISE(", hat(beta)[1], ")")), x = expression(rho) , 
       title = expression(paste("MISE(", hat(beta)[1], ") (S2)"))) + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size=13),  axis.text.x = element_text(size=13), 
        legend.title = element_text(size=14), legend.text = element_text(size=13),
        plot.title =element_text(size=16, hjust = 0.5))


ggarrange(p1, p2, ncol=2, nrow=1, align = "hv", 
          common.legend = T, legend = "right")


load("codes/simu_results/simu2S2uniT_5010100rho6.RData")
colMeans(cover_pw_fmm2d); colMeans(cover_lfosr3s); colMeans(cover_pffr); colMeans(cover_filf)
colMeans(MIAW_fmm2d); colMeans(MIAW_lfosr3s); colMeans(MIAW_pffr); colMeans(MIAW_filf)
apply(MIAW_fmm2d, 2, sd); apply(MIAW_lfosr3s, 2, sd); apply(MIAW_pffr, 2, sd); apply(MIAW_filf, 2, sd) 


################################################################################
## Time evaluation 
################################################################################
require(R.utils)
n <- c(200)
T <- c(100,200,400)
nsim <- 3
time_results <- array(0, dim=c(3,nsim, length(n), length(T)))
for(i in 1:length(n)){
  for(j in 1:length(T)){
    
    S <- 10
    type <- 2 ##beta_true type 
    SNR_B <- 1 ## relative importance of random effects
    SNR_sigma <- 1 ## signal-to-noise ratio
    family <- "gaussian" ## family of longitudinal functional data
    
    grid <- seq(0,1,length.out=T[j])
    ## self-defined phi
    phi <- array(0, dim=c(2, T[j]))
    phi[1,] <- (1.5 - sin(2*grid*pi) - cos(2*grid*pi) )
    phi[1,] <- phi[1,] / sqrt(sum(phi[1,]^2))
    phi[2,] <- sin(4*grid*pi)
    phi[2,] <- phi[2,] / sqrt(sum(phi[2,]^2))
    
    ## self-defined phi1
    phi1 <- array(0, dim=c(2, T[j]))
    phi1[1,] <- cos(2*pi*grid)
    phi1[1,] <- phi1[1,] / sqrt(sum(phi1[1,]^2))
    phi1[2,] <- sin(2*pi*grid)
    phi1[2,] <- phi1[2,] / sqrt(sum(phi1[2,]^2))
    
    ## self-defined phi2
    phi2 <- array(0, dim=c(2, T[j]))
    phi2[1,] <- sqrt(3)*(2*grid-1)
    phi2[1,] <- phi2[1,] / sqrt(sum(phi2[1,]^2))
    phi2[2,] <- sqrt(5)*(6*grid^2 - 6*grid + 1)
    phi2[2,] <- phi2[2,] / sqrt(sum(phi2[2,]^2))
    
    time_fmm2d <- time_lfosr3s <-  time_filf <- rep(0, nsim)
    for(iter in 1:nsim){
      data <- sim_generate(family = "gaussian", bivariate = F, n = n[i], S = S, T = T[j], 
                           type = type, phi = phi, phi1 = phi1, phi2 = NULL, rho = 0.5)  
    #  
    #  ptm <- proc.time()
    #  tryCatch(
    #    expr = {
    #      withTimeout({fit_fmm2d <- fmm2d(formula=Y ~ X , data=data, S=S, smoother="te", knots=c(max(round(S/4),4), min(round(T[j]/4, 35))),
    #                                      fpca.opt = list(dataType = 'Dense', methodSelectK = 'FVE'),  parallel = FALSE)}, timeout = 36000)#stop execution after one second
    #      
    #   }, 
    #    TimeoutException = function(ex) cat("Timeout. Skipping.\n")
    #  )
    #  time_fmm2d[iter] <- (proc.time() - ptm)[3]
    #  
    # if(T[j] != 400 | n[i] < 200){
    #    ptm <- proc.time()
    #    tryCatch(
     #     expr = {
    #        withTimeout({fit_lfosr3s <- lfosr3s(formula = Y ~ X + Longit  + (1 + Longit | ID), data = data, family = "gaussian", 
    #                                            var = TRUE, analytic = TRUE, parallel = FALSE, silent = FALSE)}, timeout = 36000)#stop execution after one second
    #        
    #      }, 
    #      TimeoutException = function(ex) cat("Timeout. Skipping.\n")
    #    )
    #    time_lfosr3s[iter] <- (proc.time() - ptm)[3]
    #  }else{
    #    time_lfosr3s[iter] <- NA
    #  }
      
      
      if(n[i] >=400){
        time_filf[iter] <- NA
      }else{
        T.arg <- seq(0,1, length.out = T[j])
        S.arg <- rep(seq(0,1, length.out = S), n[i])
        
        filf.x <- matrix(data$X); colnames(filf.x) <- "X"
        filf.data <- list(Y=matrix(data$Y, nrow = S*n[i]), X=filf.x, ID=data$ID)
        ptm <- proc.time()
        tryCatch(
          expr = {
            withTimeout({fit_filf <-  fgee(formula = Y ~ X, Y= filf.data$Y, Cov=filf.data$X,
                                           s=T.arg, subjID=filf.data$ID, 
                                           Tij=S.arg, numLongiPoints=S,
                                           control=list(FPCA1.pve=0.95, FPCA2.pve=0.95, FPCA1.knots= min(round(T[j]/4, 35)), FPCA2.knots=max(round(S/4),4)),
                                           corstr = "exchangeable")}, timeout = 36000)#stop execution after one second
            
          }, 
          TimeoutException = function(ex) cat("Timeout. Skipping.\n")
        )
        time_filf[iter] <- (proc.time() - ptm)[3]
      }
      
      
      
    }
    time_results[,,i,j] <- rbind(time_fmm2d, time_lfosr3s, time_filf)
  }
}

apply(time_results[,,1,3], 1, mean)/60
apply(time_results[,,1,3], 1, sd)/60

save(time_results, file="results/time.RData")


################################################################################
## Draw pointwise covariance matrix
################################################################################
n <- 50
S <- 10
T <- 100
type <- 1 ##beta_true type 
SNR_B <- 0.5 ## relative importance of random effects
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


data <- sim_generate(family = "gaussian", bivariate = T, n = n, S = S, T = T, 
                     type = type, phi = phi, phi1 = phi1, phi2 = NULL, rho = 0.001)   

#random intercept model 
fit_lfosr3s1 <- lfosr3s(formula = Y ~ X + (1 | ID), data = data, family = "gaussian", 
                       var = TRUE, analytic = TRUE, parallel = FALSE, silent = FALSE)
cov.lfosr3s <- array(0, dim=c(1000,1000))
for(s in 1:10){
  cov.lfosr3s[((s-1)*100+1):(s*100),((s-1)*100+1):(s*100)] <- fit_lfosr3s1$betaHat.var[1:100,1:100,2]
}
sub.avg <- data.frame(avg=as.vector(cov.lfosr3s), x=rep(1:1000, 1000), y=rep(1:1000, each=1000))
fig1 <- ggplot(sub.avg, aes(x, y, fill= avg)) + 
  geom_tile() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15), axis.title.y = element_blank(), 
        legend.title = element_text(size = 15), legend.key.width = unit(1, 'cm'), 
        plot.title = element_text(size = 19)) + 
  scale_fill_gradient2(low="Steel Blue 3", mid="grey98", midpoint = 0,
                       high="dark red", name = "Cov") + 
  scale_x_continuous(breaks = seq(0,1000, length.out = 6) ,labels = seq(0,1000, length.out = 6)) +
  scale_y_reverse(breaks = seq(1000,0, length.out = 6) ,labels = seq(1000,0, length.out = 6)) 

#random slope model/fmm2d
fit_fmm2d <- fmm2d(formula=Y ~ X , data=data, S=S, smoother="te", knots=c(max(round(S/4),4), min(round(T/4, 35))),
                   fpca.opt = list(dataType = 'Dense', methodSelectK = 'FVE'),  parallel = FALSE)
sub.avg <- data.frame(avg=as.vector(fit_fmm2d$betaHat.cov[1:1000,1:1000,2]), x=rep(1:1000, 1000), y=rep(1:1000, each=1000))
fig2 <- ggplot(sub.avg, aes(x, y, fill= avg)) + 
  geom_tile() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15), axis.title.y = element_blank(), 
        legend.title = element_text(size = 15), legend.key.width = unit(1, 'cm'), 
        plot.title = element_text(size = 19)) + 
  scale_fill_gradient2(low="Steel Blue 3", mid="grey98", midpoint = 0,
                       high="dark red", name = "Cov") + 
  scale_x_continuous(breaks = seq(0,1000, length.out = 6) ,labels = seq(0,1000, length.out = 6)) +
  scale_y_reverse(breaks = seq(1000,0, length.out = 6) ,labels = seq(1000,0, length.out = 6)) 



pdf("figures/2D_simu_betacov.pdf", width = 14, height = 7)
ggarrange(fig1, fig2, nrow=1, align = "hv", common.legend = T)
dev.off()











plot(beta_true_t[1,], type="l", ylim = c(-1,2))
lines(betaHat.t[1,], col="red")
lines(betaHat.t[1,] + 1.96*sqrt(diag(betaHat.t.cov[,,1])), col="red", lty=2 )
lines(betaHat.t[1,] - 1.96*sqrt(diag(betaHat.t.cov[,,1])), col="red", lty=2 )
lines(fit_lfosr3s$betaHat[1,], col="green")
lines(betaHat.up <- fit_lfosr3s$betaHat[1,] + 1.96*sqrt(diag(fit_lfosr3s$betaHat.var[,,1])), lty=2, col="green")
lines(betaHat.low <- fit_lfosr3s$betaHat[1,] - 1.96*sqrt(diag(fit_lfosr3s$betaHat.var[,,1])), lty=2, col="green")
lines(betaHat_pffr[1,], col= "blue")
lines(betaHat_pffr[1,] + 1.96*sqrt(betaHat_pffr.var[2,]), col="blue", lty=2)
lines(betaHat_pffr[1,] - 1.96*sqrt(betaHat_pffr.var[2,]), col="blue", lty=2)


plot(beta_true_t[2,], type ="l")
lines(betaHat.t[2,], col="red")
lines(betaHat.t[2,] + 1.96*sqrt(diag(betaHat.t.cov[,,2])), col="red", lty=2 )
lines(betaHat.t[2,] - 1.96*sqrt(diag(betaHat.t.cov[,,2])), col="red", lty=2 )
lines(fit_lfosr3s$betaHat[2,], col="green")
lines(fit_lfosr3s$betaHat[2,] + 1.96*sqrt(diag(fit_lfosr3s$betaHat.var[,,2])), lty=2, col="green")
lines(fit_lfosr3s$betaHat[2,] - 1.96*sqrt(diag(fit_lfosr3s$betaHat.var[,,2])), lty=2, col="green")
lines(betaHat_pffr[2,], col= "blue")
lines(betaHat_pffr[2,] + 1.96*sqrt(betaHat_pffr.var[2,]), col="blue", lty=2)
lines(betaHat_pffr[2,] - 1.96*sqrt(betaHat_pffr.var[2,]), col="blue", lty=2)
lines(fit_filf$beta[,2], col= "purple")
lines(fit_filf$beta[,2] + 1.96*fit_filf$beta.se[,2], col="purple", lty=2)
lines(fit_filf$beta[,2] - 1.96*fit_filf$beta.se[,2], col="purple", lty=2)




