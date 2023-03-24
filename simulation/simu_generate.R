#' Simulate longitudinal functional response Y (matrix) and a scalar predictor X (vector)
#' @param family distribution of longitudinal functional data, including "gaussian", "binomial", "poisson".
#' @param N number of subjects.
#' @param S number of S grid points on the functional domain.
#' @param T number of T grid points on the functional domain.
#' @param psi true orthonormal functions to generate subject-specific random intercept effects.
#' @param psi1 true orthonormal functions to generate subject specific random slope effects
#' @param psi2 true orthonormal functions to generate subject/visit-specific random effects.
#' @param SNR_B relative importance of random effects
#' @param SNR_sigma signal-to-noise ratio in Gaussian data generation
#' @export
#' @return a data frame containing generated predictors and longitudinal functional outcomes.

sim_generate <- function(family = "gaussian", bivariate = T, n = 100, S = 10, T = 100, 
                         type, phi, phi1, phi2, rho = 1.5,
                         SNR_B = 1, SNR_sigma = 1){
  library(mvtnorm)
  
  subj <- as.factor(rep(1:n, each=S))
  S.grid <- seq(0,1,length.out = S)
  T.grid <- seq(0,1,length.out = T)
  beta0 <-  intercept(S=S.grid, T=T.grid)
  if(type==1){
    beta1 <- f1(S=S.grid, T=T.grid)
  }else{
    beta1 <- f2(S=S.grid, T=T.grid)
  }

  X_des <- lapply(1:n, function(i) rnorm(S, 0, 2)) #6*(S.grid-0.5)^2 + rnorm(S, 0, 2))
  if (bivariate) {
    ## generate fixed effects
    Xbeta <- lapply(1:n, function(i) X_des[[i]]*beta1)
    fixef <- lapply(1:n, function(i) beta0 + Xbeta[[i]])
  } else {
    beta0 <- apply(beta0, 2, mean)
    beta1 <- apply(beta1, 2, mean)
    fixef <- lapply(1:n, function(i) matrix(rep(beta0,S), nrow=S, byrow = T) + X_des[[i]]%*%t(beta1))
  }
  
  ## generate random effects
  a <- rmvnorm(n, mean = rep(0, 2), sigma = diag(c(3, 1.5))) ## simulate score function
  u0 <- a %*% phi
  a1 <- rmvnorm(n, mean = rep(0, 2), sigma = diag(c(3, 1.5))) ## simulate score function
  u1 <- a1 %*% phi1
  
  Z_des <- lapply(1:n, function(i) 6*(S.grid-0.5) + rnorm(S, 0, rho))
  ranef <- lapply(1:n, function(i) t(replicate(S,u0[i,])) + t(replicate(S,u1[i,]))*Z_des[[i]])
  
  if(!is.null(phi2)){ ## by default do not add subject-visit random deviation
    c2 <- rmvnorm(n, mean = rep(0, 2), sigma = diag(c(3, 1.5)))
    eta <- c2 %*% phi2
    ranef <- lapply(1:n, function(i) ranef[[i]] + t(replicate(S, eta[i,]))) 
  }
  
  ## generate linear predictors
  ranef <- lapply(1:n, function(i) sd(fixef[[i]])/sd(ranef[[i]])/SNR_B*ranef[[i]]) ## adjust for relative importance of random effects
  Y_true <- lapply(1:n, function(i) fixef[[i]] + ranef[[i]])
  
  ## generate longitudinal functional data
  Y_obs <- lapply(1:n, function(i) array(0, dim=c(S,T)))
  for(i in 1:n){
    sd_signal <- sd(Y_true[[i]]) #sd(unlist(Y_true))
    for(s in 1:S){
      for(t in 1:T){
        if(family == "gaussian"){
          Y_obs[[i]][s, t] <- rnorm(1, mean = Y_true[[i]][s, t], sd = sd_signal/SNR_sigma)
        }
      }
    }
  }
  X <- unlist(X_des) 
  Z <- unlist(Z_des)
  Y <- do.call("rbind", Y_obs)
  dat.sim <- data.frame(ID = subj, X = X, Longit = rep(S.grid, n), Y = I(Y))
  
  return(dat.sim)
}

intercept <- function(S,T){
  intercept <- array(0,dim=c(length(S), length(T)))
  for(s in 1:length(S)){
    for(t in 1:length(T)){
      intercept[s,t] <- 3*sin(pi*(S[s]+0.5)^2)*cos(1*pi*(T[t])+0.5) + 1
    }
  }
  intercept
}

f1 <- function(S,T){
  y <- array(0,dim=c(length(S), length(T)))
  for(s in 1:length(S)){
    for(t in 1:length(T)){
      y[s,t] <- sin(0.8*pi*(S[s]+0.5)^2)*cos(4*pi*(T[t]))
    }
  }
  block <- y[(2*length(S)/10):(5*length(S)/10), (14*length(T)/100):(38*length(T)/100)]
  f <- array(0,dim=c(length(S), length(T)))
  f[(1*length(S)/10):(4*length(S)/10),(14*length(T)/100):(38*length(T)/100)] <- 5*block
  f[(7*length(S)/10):(10*length(S)/10),(14*length(T)/100):(38*length(T)/100)] <- -5*block
  f[(1*length(S)/10):(4*length(S)/10),(62*length(T)/100):(86*length(T)/100)] <- -5*block
  f[(7*length(S)/10):(10*length(S)/10),(62*length(T)/100):(86*length(T)/100)] <- 5*block
  f
}

f2 <- function(S, T){
  y <- array(0,dim=c(length(S), length(T)))
  for(s in 1:length(S)){
    for(t in 1:length(T)){
      y[s,t] <- 5*sin(0.5*pi*(S[s]+0.5)^2)*cos(2*pi*(T[t])+0.5) 
    }
  }
  y
}


#plot the 2D surface
S <- 10
T <- 100
S.grid <- seq(0,1,length.out = S)
T.grid <- seq(0,1,length.out = T)
pic <- list()
for(i in 1:3){
  if(i==1){
    f <- intercept(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))
  }else if(i==2){
    f <- f1(S=seq(0,1,length.out = S), T=seq(0,1,length.out = T))
  }else{
    f <- f2(S=seq(0,1,length.out = S), T=seq(0,1,length.out = T))
  }
  plot.subject <- data.frame(y = as.vector(f), 
                             s = rep(S.grid, T), 
                             t = rep(T.grid, each = S) )
  
  pic[[i]] <- ggplot(plot.subject, aes(x=t, y=s, fill= y)) + 
    geom_tile() + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 16), 
          legend.title = element_text(size = 14), legend.key.width = unit(1, 'cm'), 
          plot.title = element_text(size = 19)) + 
    labs(y = "S", x = "T", title = "") +
    scale_y_continuous(breaks = S.grid, labels = seq(0.1, 1, length.out = S)) +
    scale_x_continuous(breaks = seq(0,1, length.out = 5), labels =seq(0,1, length.out = 5)) +
    scale_fill_gradient2(low="#4575b4", mid="#fafafa",
                         high="#d73027", space ="Lab", name=expression(paste(beta, "(s,t)")))
}
ggarrange(pic[[1]], pic[[2]], pic[[3]], ncol=3, nrow=1)


#plot the 3D surface
library(plotly)
beta <- 0
if(beta==0){
  f <- intercept(S = seq(0,1,length.out=S), T = seq(0,1, length.out=T))
}else if(beta==1){
  f <- f1(S=seq(0,1,length.out = S), T=seq(0,1,length.out = T))
}else{
  f <- f2(S=seq(0,1,length.out = S), T=seq(0,1,length.out = T))
}
t.dir <- 1:T
s.dir <- 1:S
#true beta(s,t)
fig <- plot_ly(x = s.dir, y = t.dir, z = t(f), showscale = FALSE) %>% 
  layout(font=list(size=15), scene=list(xaxis = list(title = "S"), 
                                        yaxis = list(title = "T"),
                                        zaxis = list(title = "Beta"),
                                        camera = list(eye = list(x = -1.5, y = -1.5, z = 1.5)))) %>% 
  add_surface(colorscale = list(list(0,"#4575b4"), list(0.5,"#fafafa"), list(1,"#d73027")))
fig





