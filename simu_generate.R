#' Simulate longitudinal functional response Y (matrix) and a scalar predictor X (vector)
#' 

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





#' @param family distribution of longitudinal functional data, including "gaussian", "binomial", "poisson".
#' @param I number of subjects.
#' @param J mean number of observations per subject.
#' @param L number of grid points on the functional domain.
#' @param beta_true true fixed effects functions.
#' @param psi_true true orthonormal functions to generate subject-specific random effects.
#' @param psi2_true true orthonormal functions to generate subject/visit-specific random effects.
#' @param SNR_B relative importance of random effects
#' @param SNR_sigma signal-to-noise ratio in Gaussian data generation
#' @export
#' @return a data frame containing generated predictors and longitudinal functional outcomes.

sim_generate <- function(family = "gaussian", I = 100, S = 10, T = 100, 
                         beta_true, psi_true, psi1_true, psi2_true,
                         SNR_B = 1, SNR_sigma = 1){
  library(mvtnorm)
  
  ## generate number of visits for each subject from poisson distribution
  J_subj <- pmax(rpois(I, J), 1)
  
  ## generate fixed effects
  n <- sum(J_subj)
  X_des <- cbind(1, rnorm(n, 0, 2))
  colnames(X_des) <- c("intercept", "x")
  fixef <- X_des %*% beta_true
  
  ## generate random effects
  subj <- as.factor(rep(1:I, J_subj))
  Z_intercept <- model.matrix( ~ 0 + subj)
  Z_slope <- Z_intercept * X_des[,"x"]

  c_true <- rmvnorm(I, mean = rep(0, 2), sigma = diag(c(3, 1.5))) ## simulate score function
  b_true <- c_true %*% psi_true
  c1_true <- rmvnorm(I, mean = rep(0, 2), sigma = diag(c(3, 1.5))) ## simulate score function
  b1_true <- c1_true %*% psi1_true
  
  ranef <- Z_intercept %*% b_true + Z_slope %*% b1_true
  ## generate linear predictors
  ranef <- sd(fixef)/sd(ranef)/SNR_B*ranef ## adjust for relative importance of random effects
  
  ## by default add subject-visit random deviation
  c2_true <- rmvnorm(n, mean = rep(0, 2), sigma = diag(c(3, 1.5)))
  eta <- c2_true %*% psi2_true
  
  
  library(ggplot2)
  library(reshape2)
  cormat <- cor(eta)
  melted_cormat <- melt(cormat)
  ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill = value)) +
    geom_tile() + 
    scale_fill_continuous(trans = 'reverse')
  
  Y_true <- fixef + ranef + eta
  
  ## generate longitudinal functional data
  Y_obs <- matrix(NA, n, L) 
  p_true <- plogis(Y_true)
  lam_true <- exp(Y_true)
  sd_signal <- sd(Y_true)
  for(i in 1:n){
    for(j in 1:L){
      if(family == "gaussian"){
        Y_obs[i, j] <- rnorm(1, mean = Y_true[i, j], sd = sd_signal/SNR_sigma)
      }else if(family == "binomial"){
        Y_obs[i, j] <- rbinom(1, 1, p_true[i, j])
      }else if(family == "poisson"){
        Y_obs[i, j] <- rpois(1, lam_true[i, j])
      }
    }
  }
  
  # combine simulated data
  visit <- rep(1:I, J_subj)
  for(i in 1:I){
    visit[which(subj == i)] <- 1:J_subj[i]
  }
  dat.sim <- data.frame(ID = subj, visit = visit, X = X_des[,"x"], Y = I(Y_obs))
  
  return(dat.sim)
}
