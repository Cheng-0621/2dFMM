## Simulation study comparing FUI (implemented in the "lfosr3s" function) and FAMM (implemented in the "pffr" function)

################################################################################
## Load packages
################################################################################
rm(list = ls())
library(refund)
library(dplyr)
library(FDboost)
source("~/Documents/CityU/research/Project on Shanghai Actigraph Data/codes/simu_generate.R") ## required package: mvtnorm
source("~/Documents/CityU/research/Project on Shanghai Actigraph Data/reference/codes/FUI/lfosr3s.R") ## required packages: lme4, refund, dplyr, mgcv, progress, mvtnorm, parallel
source("~/Documents/CityU/research/Project on Shanghai Actigraph Data/reference/codes/LFPCA/LFPCA_simu.r") 

################################################################################
## Set parameters
################################################################################
I <- 50 ## number of subjects
J <- 5 ## mean number of observations per subject
L <- 50 ## dimension of the functional domain
SNR_B <- 0.5 ## relative importance of random effects
SNR_sigma <- 1 ## signal-to-noise ratio
family <- "gaussian" ## family of longitudinal functional data

## specify true fixed and random effects functions
grid <- seq(0, 1, length = L)
beta_true <- matrix(NA, 2, L)
scenario <- 1 ## indicate the scenario to use
if(scenario == 1){
  beta_true[1,] = -0.15 - 0.1*sin(2*grid*pi) - 0.1*cos(2*grid*pi)
  beta_true[2,] = dnorm(grid, .6, .15)/20
}else if(scenario == 2){
  beta_true[1,] <- 0.53 + 0.06*sin(3*grid*pi) - 0.03*cos(6.5*grid*pi)
  beta_true[2,] <- dnorm(grid, 0.2, .1)/60 + dnorm(grid, 0.35, .1)/200 - 
    dnorm(grid, 0.65, .06)/250 + dnorm(grid, 1, .07)/60
}
rownames(beta_true) <- c("Intercept", "x")

psi_true <- matrix(NA, 2, L)
psi_true[1,] <- (1.5 - sin(2*grid*pi) - cos(2*grid*pi) )
psi_true[1,] <- psi_true[1,] / sqrt(sum(psi_true[1,]^2))
psi_true[2,] <- sin(4*grid*pi)
psi_true[2,] <- psi_true[2,] / sqrt(sum(psi_true[2,]^2))

## self-defined psi1_true
psi1_true <- matrix(NA, 2, L)
psi1_true[1,] <- cos(2*pi*grid)
psi1_true[1,] <- psi1_true[1,] / sqrt(sum(psi1_true[1,]^2))
psi1_true[2,] <- sin(2*pi*grid)
psi1_true[2,] <- psi1_true[2,] / sqrt(sum(psi1_true[2,]^2))


## self-defined psi2_true
psi2_true <- matrix(NA, 2, L)
psi2_true[1,] <- sqrt(3)*(2*grid-1)
psi2_true[1,] <- psi2_true[1,] / sqrt(sum(psi2_true[1,]^2))
psi2_true[2,] <- sqrt(5)*(6*grid^2 - 6*grid + 1)
psi2_true[2,] <- psi2_true[2,] / sqrt(sum(psi2_true[2,]^2))



################################################################################
## Do simulations on a local laptop (results in the manuscript are from JHPCE)
################################################################################
nsim <- 6
sim_local <- list() ## store simulation results

for(iter in 1:nsim){
  
  ################################################################################
  ## Generate simulated data
  ################################################################################
  set.seed(iter)
  data <- sim_generate(family, I, J, L, beta_true, psi_true, psi1_true = psi1_true, psi2_true = psi2_true,
                       SNR_B = SNR_B, SNR_sigma = SNR_sigma)
  
  ################################################################################
  ## Implement different estimation methods
  ################################################################################
  ## fit the lfosr3s model
  ptm <- proc.time()
  fit_lfosr3s <- lfosr3s(formula = Y ~ X + (1 + X | ID), data = data, family = family, 
                         var = TRUE, analytic = TRUE, parallel = FALSE, silent = FALSE)
  lfosr3stime <- (proc.time() - ptm)[3]
  
  y.pred <- fit_lfosr3s$yHat
  resid <- data$Y - y.pred
  resid.test_lfosr3s <- apply(resid, 1, function(x) Box.test(x, lag = 5, type = "Ljung-Box")$p.value)
  mean(resid.test_lfosr3s < 0.05)
  
  ## fit the pffr model
  ### we use "bam" for faster computation, and change k = 15 in bs.yindex to capture curvatures in both scenarios
  ptm <- proc.time()
  fit_pffr <- pffr(formula = Y ~ X + s(ID, bs = "re") + s(X, ID, bs = "re"), data = data, family = family, algorithm = "bam",
                   bs.yindex = list(bs = "ps", k = 15, m = c(2, 1)))
  pffrtime <- (proc.time() - ptm)[3]
  y.pred <- predict(fit_pffr)
  resid <- data$Y - y.pred
  resid.test_pffr <- apply(resid, 1, function(x) Box.test(x, lag = 5, type = "Ljung-Box")$p.value)
  mean(resid.test_pffr < 0.05)

  ## fit the LFPCA model #not right
  ptm <- proc.time()
  fit_lfpca <- LFPCA(Y = data$Y,               # an n x D matrix with rows=subject-visits and columns=locations (along curve)
                     subject = data$ID,         # a vector of length n containing subject identifiers for rows of Y
                     T = data$X,               # a vector of length n containing the covariate for the random slope
                     L = 0.90,        # the pre-specified level of variance explained, determines number of components   
                     N.X = NA,        # the number of components to keep for X; N.X and N.U override L if not NA
                     N.U = NA,        # the number of components to keep for U; N.X and N.U override L if not NA
                     smooth = FALSE,  # smooth=TRUE adds smoothing of the covariance matrices not done for smooth=FALSE    
                     bf = 10          # the number of basis functions used per dimension for all smooths
  )
  lfpcatime <- (proc.time() - ptm)[3]
  
  ## fit the bayes_fosr model
  ptm <- proc.time()
  fit_bayes <- bayes_fosr(formula = Y ~ X + re(ID, by=X), data = data, Kt = 10, Kp = 10, est.method = "Gibbs", 
                          N.iter = 1000, N.burn = 500)
  bayestime <- (proc.time() - ptm)[3]
  y.pred <- fit_bayes$Yhat
  resid <- data$Y - y.pred
  resid.test_bayes <- apply(resid, 1, function(x) Box.test(x, lag = 5, type = "Ljung-Box")$p.value)
  mean(resid.test_bayes < 0.05)
  
  ## fit the FDboost
  boost <- list()
  boost$Y <- data$Y
  boost$X <- data$X
  boost$id <- as.factor(data$I)
  boost$t <- seq(1:L)
  
  ptm <- proc.time()
  fit_FDboost <- FDboost(formula = Y ~ 1 + bolsc(X, df=3) + brandomc(id, df=3) + brandomc(id, by=X, df=3),
                         data = boost, control = boost_control(mstop = 500), timeformula = ~ bbs(t, df = 4))
  y.pred <- predict(fit_FDboost)
  resid <- data$Y - y.pred
  resid.test_FDboost <- apply(resid, 1, function(x) Box.test(x, lag = 5, type = "Ljung-Box")$p.value)
  mean(resid.test_FDboost < 0.05)
  
  bootstrap <- function(data, object, num_bootstrap)
  {
    B <- num_bootstrap
    betaHat_boot <- array(NA, dim = c(2, L, B))
    ID.number <- length(unique(data$id))
    for(boots in 1:B){
      sample.ind <- sample(1:ID.number, size = ID.number, replace = TRUE)
      boost.sub <- data
      boost.sub$Y <- boost.sub$Y[boost.sub$id %in% sample.ind, ]
      boost.sub$X <- boost.sub$X[boost.sub$id %in% sample.ind]
      boost.sub$id <- boost.sub$id[boost.sub$id %in% sample.ind]
      
      fit_boot <- try(FDboost(formula = as.formula(object$formulaFDboost),
                          data = boost.sub, control = boost_control(mstop = 500), timeformula = as.formula(object$timeformula)))
      coef_boot <- coef(fit_boot, n1 = L, n2 = L)
      
      betaHat_boot[1,,boots] <- as.vector(coef_boot$smterms[[1]]$value) + coef_boot$offset$value
      betaHat_boot[2,,boots] <- predict(fit_boot, which=2, newdata=list(X=1, t=1:L))
    }
    ## obtain bootstrap variance
    betaHat.var <- array(NA, dim = c(L,L,2))
    for(r in 1:2){
      betaHat.var[,,r] <- 1.2*var(t(betaHat_boot[r,,])) ## account for within-subject correlation
    }
    return(betaHat.var)
  }
  bound <- bootstrap(data = boost, object = fit_FDboost, num_bootstrap = 10)
  FDboosttime <- (proc.time() - ptm)[3]
  

  ################################################################################
  ## Organize simulation results
  ################################################################################
  ## lfosr3s
  MISE_lfosr3s <- rowMeans((fit_lfosr3s$betaHat - beta_true)^2)
  cover_lfosr3s <- rep(NA, nrow(beta_true))
  cover_lfosr3s_pw <- matrix(FALSE, nrow(beta_true), L) 
  for(p in 1:length(cover_lfosr3s)){
    cover_upper <- which(fit_lfosr3s$betaHat[p,]+fit_lfosr3s$qn[p]*sqrt(diag(fit_lfosr3s$betaHat.var[,,p])) > beta_true[p,])
    cover_lower <- which(fit_lfosr3s$betaHat[p,]-fit_lfosr3s$qn[p]*sqrt(diag(fit_lfosr3s$betaHat.var[,,p])) < beta_true[p,])
    cover_lfosr3s[p] <- length(intersect(cover_lower, cover_upper)) == L
    cover_upper_pw <- which(fit_lfosr3s$betaHat[p,]+1.96*sqrt(diag(fit_lfosr3s$betaHat.var[,,p])) > beta_true[p,])
    cover_lower_pw <- which(fit_lfosr3s$betaHat[p,]-1.96*sqrt(diag(fit_lfosr3s$betaHat.var[,,p])) < beta_true[p,])
    cover_lfosr3s_pw[p,intersect(cover_lower_pw, cover_upper_pw)] <- TRUE
  }
  
  ## pffr
  coef_pffr <- coef(fit_pffr, n1 = L)
  betaHat_pffr <- betaHat_pffr.var <- beta_true
  betaHat_pffr[1,] <- as.vector(coef_pffr$smterms$`Intercept(yindex)`$value) + coef_pffr$pterms[1]
  betaHat_pffr[2,] <- as.vector(coef_pffr$smterms$`X(yindex)`$value)
  betaHat_pffr.var[1,] <- as.vector(coef_pffr$smterms$`Intercept(yindex)`$se^2)
  betaHat_pffr.var[2,] <- as.vector(coef_pffr$smterms$`X(yindex)`$se^2)
  MISE_pffr <- rowMeans((betaHat_pffr - beta_true)^2)
  cover_pffr <- matrix(FALSE, nrow(beta_true), L)
  for(p in 1:length(cover_lfosr3s)){
    cover_upper <- which(betaHat_pffr[p,]+1.96*sqrt(betaHat_pffr.var[p,]) > beta_true[p,])
    cover_lower <- which(betaHat_pffr[p,]-1.96*sqrt(betaHat_pffr.var[p,]) < beta_true[p,])
    cover_pffr[p,intersect(cover_lower, cover_upper)] <- TRUE
  }
  ## organize pffr results
  fit_pffr_cleaned <- list(betaHat = betaHat_pffr, betaHat.var = betaHat_pffr.var)
  
  ## lfpca
  betaHat_lfpca <- betaHat_lfpca.se <- beta_true

  #estimated time-constant mean function eta(d): apply(fit_lfpca$eta.matrix, 2, mean) 
  #equivalently: predict(fit_lfpca$eta0, newdata = data.frame(d.vec = 1:D))
  
  beta0 <- predict(fit_lfpca$eta, newdata = data.frame(d.vec = 1:L, time.vec=0), se.fit = TRUE)
  betaHat_lfpca[1,] <- beta0$fit
  beta1 <- predict(fit_lfpca$eta, newdata = data.frame(d.vec = 1:L, time.vec=1), se.fit = TRUE) 
  betaHat_lfpca[2,] <- beta1$fit - beta0$fit
  MISE_lfpca <- rowMeans((betaHat_lfpca - beta_true)^2)
  cover_lfpca <- matrix(FALSE, nrow(beta_true), L)
  betaHat_lfpca.se[1,] <- beta0$se.fit
  betaHat_lfpca.se[2,] <- beta1$se.fit
  
  for(p in 1:length(cover_lfosr3s)){
    cover_upper <- which(betaHat_lfpca[p,]+1.96*betaHat_lfpca.se[p,] > beta_true[p,])
    cover_lower <- which(betaHat_lfpca[p,]-1.96*betaHat_lfpca.se[p,] < beta_true[p,])
    cover_lfpca[p,intersect(cover_lower, cover_upper)] <- TRUE
  }
  
  ## bayes
  betaHat_bayes <- fit_bayes$beta.hat
  MISE_bayes <- rowMeans((betaHat_bayes - beta_true)^2)
  cover_bayes <- matrix(FALSE, nrow(beta_true), L)
  for(p in 1:length(cover_lfosr3s)){
    cover_upper <- which(fit_bayes$beta.UB[p,] > beta_true[p,])
    cover_lower <- which(fit_bayes$beta.LB[p,] < beta_true[p,])
    cover_bayes[p, intersect(cover_lower, cover_upper)] <- TRUE
  }
  
  ## FDboost
  betaHat_FDboost  <- array(0, dim=c(2, L))
  coef_FDboost <- coef(fit_FDboost, n1 = L, n2 = L)
  betaHat_FDboost[1,] <- as.vector(coef_FDboost$smterms[[1]]$value) + coef_FDboost$offset$value
  #betaHat_FDboost[1,] <- predict(fit_FDboost, which=1, newdata=list(X=1, t=1:L))+ coef_FDboost$offset$value
  betaHat_FDboost[2,] <- predict(fit_FDboost, which=2, newdata=list(X=1, t=1:L)) 
  #betaHat_FDboost[2,] <- apply(coef_FDboost$smterms[[2]]$value, 2, mean)
  
  FDboostbeta.LB <- FDboostbeta.UB <- array(0, dim=c(2, L))
  for (r in 1:2){
    FDboostbeta.LB[r,] <- betaHat_FDboost[r,] - 1.96*sqrt(diag(bound[,,r]))
    FDboostbeta.UB[r,] <- betaHat_FDboost[r,] + 1.96*sqrt(diag(bound[,,r]))
  }

  MISE_FDboost <- rowMeans((betaHat_FDboost - beta_true)^2)
  cover_FDboost <- matrix(FALSE, nrow(beta_true), L)
  for(p in 1:length(cover_lfosr3s)){
    cover_upper <- which(FDboostbeta.UB[p,] > beta_true[p,])
    cover_lower <- which(FDboostbeta.LB[p,] < beta_true[p,])
    cover_FDboost[p, intersect(cover_lower, cover_upper)] <- TRUE
  }
  
  
  par(mfrow=c(1,2))
  for (i in 1:2) {
    if(i == 1){
      plot(x=1:L, y=beta_true[i,], type = "l", ylim = c(-0.32, 0.05), lty = 2, lwd = 2, ylab = "beta0")
    } else {
      plot(x=1:L, y=beta_true[i,], type = "l", ylim = c(-0.1, 0.32), lty = 2, lwd = 2, ylab = "beta1")
    }
    lines(x=1:L, y=fit_lfosr3s$betaHat[i,], col="red", lwd = 2)
    lines(x=1:L, y=betaHat_pffr[i,], col="blue", lwd = 2)
    lines(x=1:L, y=betaHat_lfpca[i,], col="skyblue", lwd =2)
    lines(x=1:L, y=fit_bayes$beta.hat[i,], col="green", lwd = 2)
    lines(x=1:L, y=betaHat_FDboost[i,], col="orange", lwd = 2)
    legend("bottomright", legend=c("FAMM", "Bayes", "FDboost", "FUI"),
           col=c("blue", "green", "orange", "red"), lty=1, cex=0.8)
  }

  
  ## results of a single simulation
  sim_result <- list(MISE_lfosr3s = MISE_lfosr3s, MISE_pffr = MISE_pffr, MISE_lfpca = MISE_lfpca, MISE_bayes = MISE_bayes, MISE_FDboost = MISE_FDboost,
                     cover_lfosr3s = cover_lfosr3s, cover_lfosr3s_pw = cover_lfosr3s_pw, 
                     cover_pffr = cover_pffr, cover_lfpca = cover_lfpca, cover_bayes = cover_bayes,  cover_FDboost = cover_FDboost,
                     time_lfosr3s = lfosr3stime, time_pffr = pffrtime, time_lfpca = lfpcatime, time_bayes = bayestime, time_FDboost = FDboosttime,
                     I = I, J = J, L = L, SNR_B = SNR_B, SNR_sigma = SNR_sigma,
                     family = family, scenario = scenario)
  
  sim_local[[iter]] <- sim_result
  print(iter)
  
}



library(mgcv)
n<-200
sig <- 2
dat <- gamSim(1,n=n,scale=sig)

b<-gam(y~s(x0)+s(I(x1^2))+s(x2, x3)+offset(x3),data=dat)

newd <- data.frame(x0=(0:30)/30,x1=(0:30)/30,x2=(0:30)/30,x3=(0:30)/30)
pred <- predict.gam(b,newd)
pred0 <- predict(b,newd,exclude="s(x0)") ## prediction excluding a term
## ...and the same, but without needing to provide x0 prediction data...
newd1 <- newd;newd1$x0 <- NULL ## remove x0 from `newd1'
pred1 <- predict(b,newd1,exclude="s(x0)",newdata.guaranteed=TRUE)


################################################################################
## Obtain MISE, coverage, and computing time
################################################################################
## MISE
MISE_lfosr3s <- lapply(sim_local, '[[', 1) %>% bind_rows()
MISE_pffr <- lapply(sim_local, '[[', 2) %>% bind_rows()
MISE_lfpca <- lapply(sim_local, '[[', 3) %>% bind_rows()
MISE_bayes <- lapply(sim_local, '[[', 4) %>% bind_rows()
MISE_FDboost <- lapply(sim_local, '[[', 5) %>% bind_rows()

colMeans(MISE_lfosr3s)
colMeans(MISE_pffr)
colMeans(MISE_lfpca)
colMeans(MISE_bayes)
colMeans(MISE_FDboost)


## coverage
### joint
cover_lfosr3s <- t(lapply(sim_local, '[[', 6) %>% bind_cols())
rownames(cover_lfosr3s) <- 1:nsim
### pointwise
cover_lfosr3s_pw <- array(NA, dim = c(nsim, L, p))
for(i in 1:nsim){
  cover_lfosr3s_pw[i,,1] <- sim_local[[i]]$cover_lfosr3s_pw[1,]
  cover_lfosr3s_pw[i,,2] <- sim_local[[i]]$cover_lfosr3s_pw[2,]
}
cover_pffr_pw <- array(NA, dim = c(nsim, L, p))
for(i in 1:nsim){
  cover_pffr_pw[i,,1] <- sim_local[[i]]$cover_pffr[1,]
  cover_pffr_pw[i,,2] <- sim_local[[i]]$cover_pffr[2,]
}
cover_lfpca_pw <- array(NA, dim = c(nsim, L, p))
for(i in 1:nsim){
  cover_lfpca_pw[i,,1] <- sim_local[[i]]$cover_lfpca[1,]
  cover_lfpca_pw[i,,2] <- sim_local[[i]]$cover_lfpca[2,]
}
cover_bayes_pw <- array(NA, dim = c(nsim, L, p))
for(i in 1:nsim){
  cover_bayes_pw[i,,1] <- sim_local[[i]]$cover_bayes[1,]
  cover_bayes_pw[i,,2] <- sim_local[[i]]$cover_bayes[2,]
}
cover_FDboost_pw <- array(NA, dim = c(nsim, L, p))
for(i in 1:nsim){
  cover_FDboost_pw[i,,1] <- sim_local[[i]]$cover_FDboost[1,]
  cover_FDboost_pw[i,,2] <- sim_local[[i]]$cover_FDboost[2,]
}


### print results
print(paste0("lfosr3s coverage (joint CI): ", colMeans(cover_lfosr3s)[2]))
print(paste0("lfosr3s coverage (pointwise CI): ", apply(cover_lfosr3s_pw, 3, mean)[2]))
print(paste0("pffr coverage (pointwise CI): ", apply(cover_pffr_pw, 3, mean)[2]))
print(paste0("lfpca coverage (pointwise CI): ", apply(cover_lfpca_pw, 3, mean)[2]))
print(paste0("bayes coverage (pointwise CI): ", apply(cover_bayes_pw, 3, mean)[2]))
print(paste0("FDboost coverage (pointwise CI): ", apply(cover_FDboost_pw, 3, mean)[2]))


## time
time_lfosr3s <- lapply(sim_local, '[[', 12) %>% bind_rows()
print(apply(time_lfosr3s, 2, median))
time_pffr <- lapply(sim_local, '[[', 13) %>% bind_rows()
print(apply(time_pffr, 2, median))
time_lfpca <- lapply(sim_local, '[[', 14) %>% bind_rows()
print(apply(time_lfpca, 2, median))
time_bayes <- lapply(sim_local, '[[', 15) %>% bind_rows()
print(apply(time_bayes, 2, median))
time_FDboost <- lapply(sim_local, '[[', 16) %>% bind_rows()
print(apply(time_FDboost, 2, median))



#############################
f1 <- function(S,T, type=2){
  y <- array(0,dim=c(length(S), length(T)))
  for(s in 1:length(S)){
    for(t in 1:length(T)){
      y[s,t] <- sin(0.8*pi*(S[s]+0.5)^2)*cos(4*pi*(T[t]))
    }
  }
  block <- y[(2*length(S)/10):(5*length(S)/10), (14*length(T)/100):(38*length(T)/100)]
  f <- array(0,dim=c(length(S), length(T)))
  if(type==1){
     f[(1*length(S)/10):(4*length(S)/10),(14*length(T)/100):(38*length(T)/100)] <- 5*block
     f[(7*length(S)/10):(10*length(S)/10),(14*length(T)/100):(38*length(T)/100)] <- -5*block
     f[(1*length(S)/10):(4*length(S)/10),(62*length(T)/100):(86*length(T)/100)] <- -5*block
     f[(7*length(S)/10):(10*length(S)/10),(62*length(T)/100):(86*length(T)/100)] <- 5*block
  }
  if(type==2){
    f[(2*length(S)/10):(5*length(S)/10),(26*length(T)/100):(50*length(T)/100)] <- 5*block
    f[(6*length(S)/10):(9*length(S)/10),(26*length(T)/100):(50*length(T)/100)] <- -5*block
    f[(2*length(S)/10):(5*length(S)/10),(51*length(T)/100):(75*length(T)/100)] <- -5*block
    f[(6*length(S)/10):(9*length(S)/10),(51*length(T)/100):(75*length(T)/100)] <- 5*block
  }
  if(type==3){
    f[(1*length(S)/10):(4*length(S)/10),(14*length(T)/100):(38*length(T)/100)] <- 5*block
    f[(7*length(S)/10):(10*length(S)/10),(14*length(T)/100):(38*length(T)/100)] <- -5*block
    f[(1*length(S)/10):(4*length(S)/10),(62*length(T)/100):(86*length(T)/100)] <- -5*block
    f[(7*length(S)/10):(10*length(S)/10),(62*length(T)/100):(86*length(T)/100)] <- 5*block
    for(s in 1:length(S)){
      for(t in 1:length(T)){
        if((S[s]-0.2)^2 + (T[t]-0.25)^2 > 0.1^2 & (S[s]-0.8)^2 + (T[t]-0.25)^2 > 0.1^2 & (S[s]-0.2)^2 + (T[t]-0.75)^2 > 0.1^2 & (S[s]-0.8)^2 + (T[t]-0.75)^2 > 0.1^2) f[s,t] <- 0 
      }
    }
  }
  f
}


S.size <- 10
T.size <- 100
S.grid <- seq(0,1,length.out = S.size)
T.grid <- seq(0,1,length.out = T.size)
f.1 <- f1(S=S.grid, T=T.grid, type=3)

plot.subject <- data.frame(y = as.vector(f.1), 
                           s = rep(S.grid, T.size), 
                           t = rep(T.grid, each = S.size) )
                           
ggplot(plot.subject, aes(x=t, y=s, fill= y)) + 
  geom_tile() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 16), 
        legend.title = element_text(size = 14), legend.key.width = unit(1, 'cm'), 
        plot.title = element_text(size = 19)) + 
  labs(y = "S", x = "T", title = "") +
  scale_y_continuous(breaks = S.grid, labels = 1:S.size) +
  scale_x_continuous() +
  scale_fill_gradient2(low="blue", mid="grey98",
                       high="red", space ="Lab", name=expression(paste(beta, "(s,t)")))


psi_true <- array(0, dim=c(2, T.size))
psi_true[1,] <- (1.5 - sin(2*T.grid*pi) - cos(2*S.grid*pi) )
psi_true[1,] <- psi_true[1,] / sqrt(sum(psi_true[1,]^2))
psi_true[2,] <- sin(4*S.grid*pi)
psi_true[2,] <- psi_true[2,] / sqrt(sum(psi_true[2,]^2))

