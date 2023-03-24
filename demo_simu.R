## Simulation Demo

################################################################################
## Load packages
################################################################################
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(plotly)

source("codes/simu_generate.R") 
source("codes/2DFMM.R")

################################################################################
## Set parameters
################################################################################
n <- 50 ## number of subjects
S <- 10 ##number of observations per subject
T <- 100 ## dimension of the functional domain
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

beta_true <- array(0, dim = c(T,S,2))
beta_true[,,1] <- t(intercept(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))) #intercept


################################################################################
## S1
################################################################################
type <- 1 ##beta_true type 
data <- sim_generate(family = "gaussian", n = n, S = S, T = T, type = type, phi = phi, phi1 = phi1, phi2 = NULL, rho = 0.5)   
if(type==1){
  beta_true[,,2] <- t(f1(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))) #beta1
}else{
  beta_true[,,2] <- t(f2(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))) #beta2
}

## fit the fmm2d model 
ptm <- proc.time()
fit_S1 <- fmm2d(formula=Y~X, data=data, S=S, smoother="sandwich", knots=c(S-6, min(round(T/4), 35)),
                   fpca.opt = list(dataType = 'Dense', methodSelectK = 'FVE'),  parallel = TRUE)
time_S1 <- (proc.time() - ptm)[3]

#plot the 3D surface
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
fig <- plot_ly(x = s.dir, y = t.dir, z = fit_S1$betaHat$X, showscale = FALSE) %>% 
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
zmean <- fit_S1$betaHat$X
zvar <- diag(fit_S1$betaHat.cov[,,2])

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


################################################################################
## S2
################################################################################
type <- 2 ##beta_true type 
data <- sim_generate(family = "gaussian", n = n, S = S, T = T, type = type, phi = phi, phi1 = phi1, phi2 = NULL, rho = 0.5)   
if(type==1){
  beta_true[,,2] <- t(f1(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))) #beta1
}else{
  beta_true[,,2] <- t(f2(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))) #beta2
}

## fit the fmm2d model 
ptm <- proc.time()
fit_S2 <- fmm2d(formula=Y~X, data=data, S=S, smoother="te", knots=c(max(round(S/4),4), min(round(T/4), 35)),
                fpca.opt = list(dataType = 'Dense', methodSelectK = 'FVE'),  parallel = TRUE)
time_S2 <- (proc.time() - ptm)[3]

#plot the 3D surface 
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
fig <- plot_ly(x = s.dir, y = t.dir, z = fit_S2$betaHat$X, showscale = FALSE) %>% 
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
zmean <- fit_S2$betaHat$X
zvar <- diag(fit_S2$betaHat.cov[,,2])

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




