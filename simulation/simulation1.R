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

setwd("~/Documents/CityU/research/2DFMM/")
source("codes/simu_generate.R") 
source("codes/2DFMM.R") 

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
  fit_fmm2d <- fmm2d(formula=Y~X, data=data, S=S, smoother="te", knots=c(max(round(S/4),4), min(round(T/4), 35)),
                     fpca.opt = list(dataType = 'Dense', methodSelectK = 'FVE'),  parallel = TRUE)
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

save(MISE_fmm2d, MISE_gam2d,
     cover_pw_fmm2d, cover_pw_gam2d, 
     MIAW_fmm2d, MIAW_gam2d,
     time_fmm2d, time_gam2d,
     file = "codes/simu_results/simu1S2_5010150.RData")
load("codes/simu_results/simu1S2_5010150.RData")



est <- data.frame(y = as.vector(beta_true[,,2]), t=rep(1:T, S), s = rep(1:S, each=T)) 
ggplot(est, aes(t, s, fill= y)) + 
  geom_tile() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 16), 
        legend.title = element_text(size = 14), legend.key.width = unit(1, 'cm'), 
        plot.title = element_text(size = 19)) + 
  labs(y = "S", x = "T", title = "") +
  scale_y_continuous(breaks = 1:S ,labels = 1:S) +
  scale_x_continuous(breaks = seq(0, T, length.out = 5)) +
  scale_fill_gradient2(low="blue", mid="grey98",
                       high="red", space ="Lab", name=expression(paste(hat(beta), "(s,t)")))

betaHat <- fit_fmm2d$betaHat
est <- data.frame(y = as.vector(betaHat[[2]]), t=rep(1:T, S), s = rep(1:S, each=T)) 
ggplot(est, aes(t, s, fill= y)) + 
  geom_tile() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 16), 
        legend.title = element_text(size = 14), legend.key.width = unit(1, 'cm'), 
        plot.title = element_text(size = 19)) + 
  labs(y = "S", x = "T", title = "") +
  scale_y_continuous(breaks = 1:S ,labels = 1:S) +
  scale_x_continuous(breaks = seq(0, T, length.out = 5)) +
  scale_fill_gradient2(low="blue", mid="grey98",
                       high="red", space ="Lab", name=expression(paste(hat(beta), "(s,t)")))

betaHat.var <- fit_fmm2d$betaHat.cov
betaHat.up <- as.vector(betaHat[[2]]) + 2*sqrt(diag(betaHat.var[,,2]))
betaHat.low <- as.vector(betaHat[[2]]) - 2*sqrt(diag(betaHat.var[,,2]))

est <- data.frame(betaHat.low=betaHat.low, betaHat.up=betaHat.up, t=rep(1:T, S), s = rep(1:S, each=T)) 
est$sign <- ifelse(est$betaHat.up < 0,  est$betaHat.up, 
                            ifelse(est$betaHat.low > 0, est$betaHat.low, 0))

ggplot(est, aes(t, s, fill= sign)) + 
  geom_tile() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 16), 
        legend.title = element_text(size = 14), legend.key.width = unit(1, 'cm'), 
        plot.title = element_text(size = 19)) + 
  labs(y = "S", x = "T", title = "") +
  scale_y_continuous(breaks = 1:S ,labels = 1:S) +
  scale_x_continuous(breaks = seq(0, T, length.out = 5)) +
  scale_fill_gradient2(low="blue", mid="grey98", 
                       high="red", space ="Lab", name=expression(paste(hat(beta), "(s,t)")))


est <- data.frame(y = fit_gam2d$X[,"est"], t=rep(1:T, S), s = rep(1:S, each=T)) 
ggplot(est, aes(x=t, y=s, fill= est)) + 
  geom_tile() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 16), 
        legend.title = element_text(size = 14), legend.key.width = unit(1, 'cm'), 
        plot.title = element_text(size = 19)) + 
  labs(y = "S", x = "T", title = "") +
  scale_y_continuous(breaks = 1:S ,labels = 1:S) +
  scale_x_continuous(breaks = seq(0, T, length.out = 5)) +
  scale_fill_gradient2(low="blue", mid="grey98",
                       high="red", space ="Lab", name=expression(paste(hat(beta), "(s,t)")))



est[,"betaHat.up"] <- fit_gam2d$X[,"est"] + 2*fit_gam2d$X[,"se"]
est[,"betaHat.low"] <- fit_gam2d$X[,"est"] - 2*fit_gam2d$X[,"se"]
est[,"sign"] <- ifelse(est$betaHat.up < 0,  est$betaHat.up, 
                   ifelse(est$betaHat.low > 0, est$betaHat.low, 0))

ggplot(est, aes(x=t, y=s, fill= sign)) + 
  geom_tile() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 16), 
        legend.title = element_text(size = 14), legend.key.width = unit(1, 'cm'), 
        plot.title = element_text(size = 19)) + 
  labs(y = "S", x = "T", title = "") +
  scale_y_continuous(breaks = 1:S ,labels = 1:S) +
  scale_x_continuous(breaks = seq(0, T, length.out = 5)) +
  scale_fill_gradient2(low="blue", mid="grey98",
                       high="red", space ="Lab", name=expression(paste(hat(beta), "(s,t)")))




#plot the 3D surface
library(plotly)

beta_true <- array(0, dim = c(T,S,2))
beta_true[,,1] <- t(intercept(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))) #intercept
beta_true[,,2] <- t(f2(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))) #beta1

t.dir <- 1:T
s.dir <- 1:S

#true beta(s,t)
fig <- plot_ly(x = s.dir, y = t.dir, z = beta_true[,,1], showscale = FALSE) %>% 
  layout(font=list(size=13), scene=list(xaxis = list(title = "S"), 
                                        yaxis = list(title = "T"),
                                        zaxis = list(title = "Coef"), 
                                        camera = list(eye = list(x = -1.25, y = -1.25, z = 1.5)))) %>% 
  add_surface(colorscale = list(list(0,"#4575b4"), list(0.5,"#fafafa"), list(1,"#d73027")))
fig

#betahat(s,t)
fig <- plot_ly(x = s.dir, y = t.dir, z = fit_fmm2d$betaHat$X, showscale = FALSE) %>% 
  config(
    toImageButtonOptions = list(
      format = "jpeg",
      filename = "",
      width = 640,
      height = 480, 
      scale = 5
    )
  ) %>%
  layout(font=list(size=13), scene=list(xaxis = list(title = "S"), 
                                        yaxis = list(title = "T"),
                                        zaxis = list(title = "Coef"), 
                                        camera = list(eye = list(x = -1.25, y = -1.35, z = 1.5)))) %>% 
  add_surface(colorscale = "Surface") 
fig


#CI betahat(s,t)
#zmean <- array(as.matrix(fit_gam2d$X[,"est"]), dim=c(T,S)) #
#zvar <- array(as.matrix(fit_gam2d$X[,"se"]), dim=c(T,S))^2 #
zmean <- fit_fmm2d$betaHat$X
zvar <- diag(fit_fmm2d$betaHat.cov[,,2])

fig <- plot_ly(x = s.dir, y = t.dir, 
               z = zmean + array(2*sqrt(zvar),dim=c(T,S)), showscale = FALSE) %>% 
  layout(font=list(size=13), scene=list(xaxis = list(title = "S"), 
                                        yaxis = list(title = "T"),
                                        zaxis = list(title = "Coef"), 
                                        camera = list(eye = list(x = -1.25, y = -1.25, z = 1.5)))) %>% 
  add_surface() 
fig <- fig %>% 
  add_surface(x = s.dir, y = t.dir, z = zmean - array(2*sqrt(as.vector(zvar)),dim=c(T,S)), showscale = FALSE) 
fig <- fig %>% add_surface(x = s.dir, y = t.dir, z = beta_true[,,2], showscale = FALSE, colorscale = "Surface") 

fig




#load results
load_results <- function(arg, val, dir, iftime = FALSE){
  
  load(paste0("codes/simu_results/", dir[1], ".RData"))
  MISEf1 <- cbind(MISE_fmm2d, method="FMM2d", arg=val[1])
  MISEg1 <- cbind(MISE_gam2d, method="GAM2d", arg=val[1])
  coverf1 <- cbind(cover_pw_fmm2d, method="FMM2d", arg=val[1])
  coverg1 <- cbind(cover_pw_gam2d, method="GAM2d", arg=val[1])
  timef1 <- cbind(time_fmm2d, method="FMM2d", arg=val[1])
  timeg1 <- cbind(time_gam2d, method="GAM2d", arg=val[1])
  MIAWf1 <- cbind(MIAW_fmm2d, method="FMM2d", arg = val[1])
  MIAWg1 <- cbind(MIAW_gam2d, method="GAM2d", arg = val[1])
  
  load(paste0("codes/simu_results/", dir[2], ".RData"))
  MISEf2 <- cbind(MISE_fmm2d, method="FMM2d", arg=val[2])
  MISEg2 <- cbind(MISE_gam2d, method="GAM2d", arg=val[2])
  coverf2 <- cbind(cover_pw_fmm2d, method="FMM2d", arg=val[2])
  coverg2 <- cbind(cover_pw_gam2d, method="GAM2d", arg=val[2])
  timef2 <- cbind(time_fmm2d, method="FMM2d", arg=val[2])
  timeg2 <- cbind(time_gam2d, method="GAM2d", arg=val[2])
  MIAWf2 <- cbind(MIAW_fmm2d, method="FMM2d", arg = val[2])
  MIAWg2 <- cbind(MIAW_gam2d, method="GAM2d", arg = val[2])
  
  load(paste0("codes/simu_results/", dir[3], ".RData"))
  MISEf3 <- cbind(MISE_fmm2d, method="FMM2d", arg=val[3])
  MISEg3 <- cbind(MISE_gam2d, method="GAM2d", arg=val[3])
  coverf3 <- cbind(cover_pw_fmm2d, method="FMM2d", arg=val[3])
  coverg3 <- cbind(cover_pw_gam2d, method="GAM2d", arg=val[3])
  timef3 <- cbind(time_fmm2d, method="FMM2d", arg=val[3])
  timeg3 <- cbind(time_gam2d, method="GAM2d", arg=val[3])
  MIAWf3 <- cbind(MIAW_fmm2d, method="FMM2d", arg = val[3])
  MIAWg3 <- cbind(MIAW_gam2d, method="GAM2d", arg = val[3])
  
  
  MISE <- rbind(MISEf1, MISEf2, MISEf3, MISEg1, MISEg2, MISEg3)
  MISE$method <- factor(MISE$method, levels = c("FMM2d", "GAM2d"))
  MISE$direction <- arg
  cover <- rbind(coverf1, coverf2, coverf3, coverg1, coverg2, coverg3)
  cover$method <- factor(cover$method, levels = c("FMM2d", "GAM2d"))
  cover$direction <- arg
  MIAW <- rbind(MIAWf1, MIAWf2, MIAWf3, MIAWg1, MIAWg2, MIAWg3)
  MIAW$method <- factor(MIAW$method, levels = c("FMM2d", "GAM2d"))
  MIAW$direction <- arg
  time <- rbind(timef1, timef2, timef3, timeg1, timeg2, timeg3)
  time$method <- factor(MISE$method, levels = c("FMM2d", "GAM2d"))
  time$direction <- arg

  ret <- list(MISE=MISE, cover=cover, MIAW=MIAW, time=time)
  return(ret)
}


#direction <- "Sample size (N)" 
direction <- "Number of grids (R)"
#direction <- "Number of grids (L)"
val <- c(10, 15, 20)
ret <- load_results(arg=direction, val=val, 
                    dir=c("simu1S2_5010100", "simu1S2_5015100","simu1S2_5020100"))


p1 <- ggplot(ret$MISE, aes(x=arg, y=beta, fill=method, group = paste(method,arg))) + 
  geom_boxplot() + 
  theme_bw() +
  scale_fill_discrete(name = "Method", labels = c("FMM2d", "GAM2d")) + 
  labs(y = expression(paste("MISE(", hat(beta)[1], ")")), x = direction, 
       title = expression(paste("MISE(", hat(beta)[1], ") (S2)"))) + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size=13),  axis.text.x = element_text(size=13), 
        legend.title = element_text(size=14), legend.text = element_text(size=13),
        plot.title =element_text(size=16, hjust = 0.5)) +
  scale_x_continuous(breaks=val, labels=val)

p2 <- ggplot(ret$cover, aes(x=arg, y=beta, fill = method, group = paste(method,arg))) + 
  geom_boxplot() + 
  theme_bw() +
  scale_fill_discrete(name = "Method", labels = c("FMM2d", "GAM2d")) + 
  labs(y = expression(paste("Coverage(", hat(beta)[1], ")")), x = direction, 
       title = expression(paste("Coverage(", hat(beta)[1], ") (S2)"))) + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size=13),  axis.text.x = element_text(size=13), 
        legend.title = element_text(size=14), legend.text = element_text(size=13),
        plot.title =element_text(size=16, hjust = 0.5)) +
  scale_x_continuous(breaks=val, labels=val)

p3 <- ggplot(ret$MIAW, aes(x=arg, y=beta, fill = method, group = paste(method,arg))) + 
  geom_boxplot() + 
  theme_bw() +
  scale_fill_discrete(name = "Method", labels = c("FMM2d", "GAM2d")) + 
  labs(y = expression(paste("MIAW(", hat(beta)[1], ")")), x = direction, 
       title = expression(paste("MIAW(", hat(beta)[1], ") (S2)"))) +
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size=12),  axis.text.x = element_text(size=12), 
        legend.title = element_text(size=14), legend.text = element_text(size=13),
        plot.title =element_text(size=16, hjust = 0.5)) + 
  scale_x_continuous(breaks=val, labels=val)

p4 <- ggplot(ret$time, aes(x=arg, y=elapsed/3600, fill = method, group = paste(method,arg))) + 
  geom_boxplot() + 
  theme_bw() +
  scale_fill_discrete(name = "Method", labels = c("FMM2d", "GAM2d")) + 
  labs(y = "Computing Time (hours)", x = direction, 
       title = "Computing Time (S2)") + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size=13),  axis.text.x = element_text(size=13), 
        legend.title = element_text(size=14), legend.text = element_text(size=13),
        plot.title =element_text(size=16, hjust = 0.5)) + 
  scale_x_continuous(breaks=val, labels=val)


ggarrange(p1, p2, p3, p4, ncol=4, nrow=1, align = "hv",  
          common.legend = T, legend = "right")




















