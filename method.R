#' Fit longitudinal function-on-scalar regression using the proposed 3-step approach
#' 
#' @param formula two-sided formula object in lmer() format, except that the response is a matrix  
#' @param data data frame containing variables in formula
#' @param family GLM family of the response
#' @param argvals locations of observations on the functional domain
#' @param var whether to estimate variance. By default FALSE
#' @param analytic whether to use the analytic inference approach
#' @param parallel whether to run parallel computing
#' @param silent whether to show descriptions of each step
#' @export
#' @return a list containing estimated beta(s)

our.lfosr3s <- function(formula, data, family = "gaussian", argvals = NULL,
                        bs = "cr", nknots = 35, raw = FALSE,
                        var = FALSE, analytic = TRUE, parallel = FALSE, silent = FALSE){
  library(lme4) ## mixed models
  library(refund) ## fpca.face
  library(dplyr) ## organize lapply results
  library(progress) ## display progress bar
  library(mgcv) ## smoothing in step 2
  library(mvtnorm) ## joint CI
  library(parallel) ## mcapply
  
  if(family != "gaussian") analytic <- FALSE ## bootstrap inference for non-Gaussian family
  
  ## Organize the input
  model_formula <- as.character(formula)
  stopifnot(model_formula[1] == "~" & length(model_formula) == 3)
  
  
  ##########################################################################################
  ## Step 1
  ##########################################################################################
  if(silent == FALSE) print("Step 1: Massive Univariate Mixed Models")
  
  L <- ncol(data[,model_formula[2]]) ## number of observations on the functional domain
  if(is.null(argvals)) argvals <- 1:L
  
  ## function "unimm" fit univariate mixed model at location l
  unimm <- function(l){
    #print(l)
    data$Yl <- unclass(data[,model_formula[2]][,l])
    if(family == "gaussian"){
      fit_uni <- suppressMessages(lmer(formula = as.formula(paste0("Yl ~ ", model_formula[3])), 
                                       data = data, control = lmerControl(optimizer = "bobyqa")))
    }else{
      fit_uni <- suppressMessages(glmer(formula = as.formula(paste0("Yl ~ ", model_formula[3])), 
                                        data = data, family = family, control = glmerControl(optimizer = "bobyqa")))
    }
    betaTilde <- lme4::fixef(fit_uni)
    #uTilde <- lme4::ranef(fit_uni, drop=TRUE)[[1]]
    yTilde <- predict(fit_uni)
    if(analytic == TRUE){
      varcorr <- as.data.frame(VarCorr(fit_uni))
      var_random <- varcorr[,4] ## variance/covariance estimates
      ind_var <- which(is.na(varcorr[,3]) & varcorr[,1] != "Residual") ## variance of random components
      names(var_random)[ind_var] <- paste0("var.",varcorr[ind_var,1],".",varcorr[ind_var,2])
      names(var_random)[which(varcorr[,1] == "Residual")] <- "var.Residual"
      names(var_random)[which(!is.na(as.data.frame(VarCorr(fit_uni))[,3]))] <- "cov"
      return(list(betaTilde = betaTilde, yTilde = yTilde, var_random = var_random, group = varcorr[1,1]))
    }else{
      return(list(betaTilde = betaTilde, yTilde = yTilde, group = as.data.frame(VarCorr(fit_uni))[1,1]))
    }
  }
  
  ## fit massive univariate mixed models
  if(parallel == TRUE){
    massmm <- mclapply(argvals, unimm, mc.cores = detectCores() - 1)
  }else{
    massmm <- lapply(argvals, unimm)
  }
  
  ## obtain betaTilde
  betaTilde <- t(lapply(massmm, '[[', 1) %>% bind_rows())
  colnames(betaTilde) <- argvals
  ## obtain uTilde
  #group <- data.frame(massmm[[1]]$uTilde)
  #uTilde <- array(NA, dim = c(nrow(group),L,ncol(group)))
  #for (l in 1:dim(uTilde)[2]){
  #  for (r in 1:dim(uTilde)[3]){
  #    uTilde[,l,r] <- data.frame(massmm[[l]]$uTilde)[,r] 
  #  }
  #}
  ## obtain yTilde
  yTilde <- t(lapply(massmm, '[[', 2) %>% bind_rows())
  
  ## obtain variance estimates of random part when necessary
  if(analytic == TRUE){
    var_random <- t(lapply(massmm, '[[', 3) %>% bind_rows())
    sigmaesqHat <- var_random["var.Residual",,drop = FALSE]
    sigmausqHat <- var_random[which(rownames(var_random) != "var.Residual"),,drop = FALSE]
    rm(var_random)
    
    ## fit a fake model to obtain design matrix
    data$Yl <- unclass(data[,model_formula[2]][,1])
    if(family == "gaussian"){
      fit_uni <- suppressMessages(lmer(formula = as.formula(paste0("Yl ~ ", model_formula[3])), 
                                       data = data, control = lmerControl(optimizer = "bobyqa")))
    }else{
      fit_uni <- suppressMessages(glmer(formula = as.formula(paste0("Yl ~ ", model_formula[3])), 
                                        data = data, family = family, control = glmerControl(optimizer = "bobyqa")))
    }
    designmat <- model.matrix(fit_uni) ## model design matrix, necessary for analytic inference
    name_random <- as.data.frame(VarCorr(fit_uni))[which(!is.na(as.data.frame(VarCorr(fit_uni))[,3])),3]
    rm(fit_uni)
  }
  
  
  ##########################################################################################
  ## Step 2
  ##########################################################################################
  if(silent == FALSE) print("Step 2: Smoothing")
  
  nknots <- nknots #min(round(L/4), 35) ## number of knots for penalized splines smoothing
  if(analytic == TRUE){ ## we need smoothing parameter, spline basis and penalty matrix for analytic inference
    p <- nrow(betaTilde) ## number of fixed effects parameters
    betaHat <- matrix(NA, nrow = p, ncol = L)
    lambda <- rep(NA, p) ## smoothing parameter
    for(r in 1:p){
      fit_smooth <- gam(betaTilde[r,] ~ s(argvals, bs = bs, k = (nknots + 1)), method = "REML")
      betaHat[r,] <- fit_smooth$fitted.values
      lambda[r] <- fit_smooth$sp ## get smoothing parameter
    }
    #s <- dim(uTilde)[3]
    #group <- data.frame(massmm[[1]]$uTilde)
    #uHat <- array(NA, dim = c(nrow(group),L,ncol(group)))
    #for(r in 1:s){
    #  uHat[,,r] <- t(apply(uTilde[,,r], 1, function(x) gam(x ~ s(argvals, bs = "cr", k = (nknots + 1)), method = "REML")$fitted.values))
    #}
    
    sm <- smoothCon(s(argvals, bs = "cr", k = (nknots + 1)), data=data.frame(argvals=argvals), absorb.cons=TRUE)
    S<- sm[[1]]$S[[1]] ## penalty matrix
    B <- sm[[1]]$X ## basis functions
    
    rm(fit_smooth, sm)
  }else{
    betaHat <- t(apply(betaTilde, 1, function(x) gam(x ~ s(argvals, bs = "cr", k = (nknots + 1)), method = "REML")$fitted.values))
  }
  rownames(betaHat) <- rownames(betaTilde)
  colnames(betaHat) <- 1:L
  
  yHat <- designmat %*% betaHat #yHat 
  
  ##########################################################################################
  ## Step 3
  ##########################################################################################
  if(var == TRUE){
    if(analytic == TRUE){ ## analytic inference
      
      ##########################################################################################
      ## Analytic Inference
      ##########################################################################################
      if(silent == FALSE) print("Step 3: Analytic Inference")
      
      ################ Preparation
      if(silent == FALSE) print("Step 3.1: Preparation")
      ## 1. derive variance estimates of random components (\hat{H}(s), \hat{R}(s))
      ### smooth raw estimates
      HHat <- t(apply(sigmausqHat, 1, function(b) smooth.spline(x = argvals, y = b)$y))
      HHat[which(grepl("var", rownames(HHat)) == TRUE),][which(HHat[which(grepl("var", rownames(HHat)) == TRUE),] < 0)] <- 0
      RHat <- t(apply(sigmaesqHat, 1, function(b) smooth.spline(x = argvals, y = b)$y))
      RHat[which(RHat < 0)] <- 0
      if(silent == FALSE) print("Step 3.1.1: variance of random components")
      
      ## 2. derive covariance estimates of random components (\hat{G}(s1,s2))
      if(length(which(rownames(HHat) == "cov")) == 0){
        GTilde <- matrix(NA, nrow = L, ncol = L)
        for(i in 1:L){
          for(j in 1:L){ # no need to begin from 1, can begin from i
            GTilde[i,j] <- cov(data[,model_formula[2]][,i], data[,model_formula[2]][,j], use = "pairwise.complete.obs") -
              t(betaHat[,i]) %*% var(designmat) %*% betaHat[,j]
            
            #print(GTilde[i,j])
            
            #GTilde[i,j] <- (1/n) * t(data[,model_formula[2]][,i] - designmat %*% betaHat[,i]) %*%
            #  (data[,model_formula[2]][,j] - designmat %*% betaHat[,j])
            
            #print(GTilde[i,j])
            
            #print(paste(i, j, "MoM done"))
          }
        }
        diag(GTilde) <- HHat[1,]
        GHat <- fbps(GTilde)$Yhat ## fast bivariate smoother
        diag(GHat)[which(diag(GHat) < 0)] <- diag(GTilde)[which(diag(GHat) < 0)]
      }else{ ## fit a simple regression model between each of location pairs
        GTilde <- array(NA, dim = c(nrow(HHat), L, L))
        for(i in 1:L){
          for(j in i:L){ 
            data_cov <- data.frame(Y = (data[,model_formula[2]][,i] - designmat %*% matrix(betaHat[,i], ncol = 1)) *
                                     (data[,model_formula[2]][,j] - designmat %*% matrix(betaHat[,j], ncol = 1)),
                                   twoZ = 2 * data[,name_random], Zsq = data[,name_random]^2) 
            #twoZ is intercept*random variable; sqaureZ is random variable*random variable 
            fit_cov <- lm(Y ~ Zsq + twoZ, data = data_cov)
            GTilde[,i,j] <- GTilde[,j,i] <- fit_cov$coefficients
          }
        }
        rm(data_cov, fit_cov)
        GHat <- GTilde
        for(r in 1:nrow(HHat)){
          GHat[r,,] <- fbps(GTilde[r,,])$Yhat
        }
      }
      if(silent == FALSE) print("Step 3.1.2: covariance of random components")
      
      ## 3. variance estimates at each location (V(s))
      obs.ind <- list()
      group <- massmm[[1]]$group ## group name in the data
      for(id in unique(data[,group])){
        obs.ind[[as.character(id)]] <- which(data[,group] == id) ## the corresponding rows of each subject
      }
      V.inv.all <- list()
      for(s in 1:L){
        V.subj.inv <- list()
        ## we first do inverse of each block matrix then combine them
        if(!length(which(rownames(HHat) == "cov")) == 0){
          cov.raw <- matrix(c(HHat[1,s], HHat[3,s], HHat[3,s], HHat[2,s]), nrow = 2, ncol = 2)
          edcomp <- eigen(cov.raw) ## trim non-positive eigenvalues to ensure positive semidefinite
          eigen.positive <- which(edcomp$values > 0)
          if(length(eigen.positive) == 2){
            cov.trimmed <- cov.raw
          }else{
            cov.trimmed <- matrix(edcomp$vectors[,1], ncol = 1) %*% edcomp$values[1] %*% matrix(edcomp$vectors[,1], nrow = 1)
          }
        }
        for(id in unique(data[,group])){
          if(length(which(rownames(HHat) == "cov")) == 0){
            Ji <- length(obs.ind[[as.character(id)]])
            V.subj <- matrix(1, nrow = Ji, ncol = 1) %*% HHat[1,s] %*% matrix(1, nrow = 1, ncol = Ji) + diag(RHat[s], Ji)
          }else{
            subj.ind <- obs.ind[[as.character(id)]]
            V.subj <- cbind(rep(1, length(subj.ind)), data[subj.ind, name_random]) %*% cov.trimmed %*%
              t(cbind(rep(1, length(subj.ind)), data[subj.ind, name_random])) + diag(RHat[s], length(subj.ind))
          }
          V.subj.inv[[as.character(id)]] <- solve(V.subj) ## calculate inverse of each block
        }
        V.inv.all[[s]] <- V.subj.inv
      }
      if(silent == FALSE) print("Step 3.1.3: V(s)")
      
      ## 4. covariance between any two pairs on the functional domain (V(s_1, s_2)) #it is actually W(s_1, s_2)
      V.cov.all <- list()
      for(i in 1:L){
        for(j in i:L){
          V.cov.subj <- list()
          if(!length(which(rownames(HHat) == "cov")) == 0){
            cov.raw <- matrix(c(GHat[1,i,j], GHat[3,i,j], GHat[3,i,j], GHat[2,i,j]), nrow = 2, ncol = 2)
            edcomp <- eigen(cov.raw)
            eigen.positive <- which(edcomp$values > 0)
            if(length(eigen.positive) == 2){
              cov.trimmed <- cov.raw
            }else{
              cov.trimmed <- matrix(edcomp$vectors[,1], ncol = 1) %*% edcomp$values[1] %*% matrix(edcomp$vectors[,1], nrow = 1)
            }
          }
          #print("done trimmed")
          for(id in unique(data[,group])){
            if(length(which(rownames(HHat) == "cov")) == 0){
              Ji <- length(obs.ind[[as.character(id)]])
              V.cov.subj[[as.character(id)]] <- matrix(1, nrow = Ji, ncol = 1) %*% GHat[i,j] %*% matrix(1, nrow = 1, ncol = Ji)
            }else{
              subj.ind <- obs.ind[[as.character(id)]]
              V.cov.subj[[as.character(id)]] <- cbind(rep(1, length(subj.ind)), data[subj.ind, name_random]) %*% cov.trimmed %*%
                t(cbind(rep(1, length(subj.ind)), data[subj.ind, name_random]))
            }
          }
          #print(paste(i,j, "W(s1,s2) done"))
          V.cov.all[[L*(i-1)+j]] <- V.cov.subj
        }
      }
      rm(V.subj, V.subj.inv, V.cov.subj)
      if(silent == FALSE) print("Step 3.1.4: W(s1, s2)")
      
      ################ First step
      if(silent == FALSE) print("Step 3.2: First step")
      ## 1. intra-location variance (Var(\Tilde{\beta}(s)))
      var.beta.tilde.theo <- array(NA, dim = c(p,p,L))
      for(s in 1:L){
        tmp <- matrix(0, nrow = p, ncol = p)
        for(id in unique(data[,group])){
          tmp <- tmp + t(matrix(designmat[obs.ind[[as.character(id)]],], ncol = p)) %*%
            V.inv.all[[s]][[as.character(id)]] %*% matrix(designmat[obs.ind[[as.character(id)]],], ncol = p)
        }
        var.beta.tilde.theo[,,s] <- solve(tmp)
      }
      if(silent == FALSE) print("Step 3.2.1:intra-location variance (Var(Tilde{beta}(s)))")
      
      
      ## 2. inter-location covariance (Cov(\Tilde{\beta}(s_1), \Tilde{\beta}(s_2)))
      cov.beta.tilde.theo <- array(NA, dim = c(p,p,L,L))
      for(i in 1:L){
        for(j in i:L){
          tmp <- matrix(0, nrow = p, ncol = p)
          for(id in unique(data[,group])){
            tmp <- tmp + t(matrix(designmat[obs.ind[[as.character(id)]],], ncol = p)) %*% V.inv.all[[i]][[as.character(id)]] %*%
              V.cov.all[[L*(i-1)+j]][[as.character(id)]] %*% V.inv.all[[j]][[as.character(id)]] %*% matrix(designmat[obs.ind[[as.character(id)]],], ncol = p)
          }
          cov.beta.tilde.theo[,,i,j] <- var.beta.tilde.theo[,,i] %*% tmp %*% var.beta.tilde.theo[,,j]
        }
      }
      rm(V.cov.all, V.inv.all)
      if(silent == FALSE) print("Step 3.2.2:inter-location covariance (Cov(Tilde{beta}(s_1), Tilde{beta}(s_2))")
      
      
      ################ Second step
      if(silent == FALSE) print("Step 3.3: Second step")
      ## 1. intermediate step for covariance estimate
      var.beta.tilde.s <- array(NA, dim = c(L,L,p))
      for(j in 1:p){
        for(r in 1:L){
          for(t in 1:L){
            if(t == r){
              var.beta.tilde.s[r,t,j] <- var.beta.tilde.theo[j,j,r]
            }else{
              var.beta.tilde.s[r,t,j] <- cov.beta.tilde.theo[j,j,min(r,t),max(r,t)]
            }
          }
        }
      }
      
      var.betaTilde <- array(NA, dim = c(p, L))
      for (j in 1:p){
        var.betaTilde[j,] <- diag(var.beta.tilde.s[,,j])
      }
      
      ## 2. final variance estimate of betaHat (Cov(\hat{\beta}(s_1), \hat{\beta}(s_2)))
      var.beta.hat <- array(NA, dim = c(L,L,p))
      for(r in 1:p){
        M <- B %*% solve(t(B)%*%B + lambda[r]*S) %*% t(B) + matrix(1/L, nrow = L, ncol = L)
        var.raw <- M %*% var.beta.tilde.s[,,r] %*% t(M)
        ## trim eigenvalues to make final variance matrix positive-semidefinite
        edcomp <- eigen(var.raw)
        eigen.positive <- which(edcomp$values > 0)
        var.beta.hat[,,r] <- edcomp$vectors[,eigen.positive] %*% diag(edcomp$values[eigen.positive]) %*% t(edcomp$vectors[,eigen.positive])
      }
      betaHat.var <- var.beta.hat ## final variance estimate
      
      
      ## obtain qn to construct joint CI
      qn <- rep(0, length = nrow(betaHat))
      N <- 10000 ## sample size in simulation-based approach
      for(i in 1:length(qn)){
        Sigma <- betaHat.var[,,i]
        x_sample <- rmvnorm(N, mean = betaHat[i,], sigma = Sigma)
        un <- rep(NA, N) 
        for(j in 1:N){
          un[j] <- max(abs((x_sample[j,] - betaHat[i,])/sqrt(diag(Sigma))))
        }
        qn[i] <- quantile(un, 0.95)
      }
      
      if (raw) {
        return(list(betaHat = betaTilde, betaHat.var = var.betaTilde, yHat = yTilde))
      } else {
        return(list(betaHat = betaHat, betaHat.var = betaHat.var, yHat = yHat, qn = qn))
      }
      
    }else{ ## bootstrap inference
      
      ##########################################################################################
      ## Bootstrap Inference
      ##########################################################################################
      if(silent == FALSE) print("Step 3: Bootstrap Inference")
      
      B <- 100
      betaHat_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), B))
      group <- massmm[[1]]$group
      ID.number <- unique(data[,group])
      pb <- progress_bar$new(total = B)
      for(boots in 1:B){
        pb$tick()
        sample.ind <- sample(1:length(ID.number), size = length(ID.number), replace = TRUE)
        dat.ind <- c()
        for(i in 1:length(ID.number)){
          dat.ind <- c(dat.ind, which(data[,group] == ID.number[sample.ind[i]])) ## subject-level bootstrap
        }
        fit_boot <- lfosr3s(formula = formula, data = data[dat.ind,], family = family, var = FALSE, 
                            parallel = parallel, silent = TRUE)
        betaHat_boot[,,boots] <- fit_boot$betaHat
      }
      
      ## obtain bootstrap variance
      betaHat.var <- array(NA, dim = c(L,L,nrow(betaHat)))
      for(r in 1:nrow(betaHat)){
        betaHat.var[,,r] <- 1.2*var(t(betaHat_boot[r,,])) ## account for within-subject correlation
      }
      
      ## obtain qn to construct joint CI using the fast approach
      qn <- rep(0, length = nrow(betaHat))
      N <- 10000 ## sample size in simulation-based approach
      for(i in 1:length(qn)){
        est_bs <- t(betaHat_boot[i,,])
        fit_fpca <- fpca.face(est_bs)
        ## extract estimated eigenfunctions/eigenvalues
        phi <- fit_fpca$efunctions
        lambda <- fit_fpca$evalues
        K <- length(fit_fpca$evalues)
        ## simulate random coefficients 
        theta <- matrix(rnorm(N*K), nrow=N, ncol=K) # generate independent standard normals 
        theta <- theta %*% diag(sqrt(lambda)) # scale to have appropriate variance
        X_new <- theta %*% t(phi) # simulate new functions
        x_sample <- X_new + t(fit_fpca$mu %o% rep(1,N)) # add back in the mean function
        Sigma <- apply(x_sample, 2, var)
        x_mean <- colMeans(est_bs)
        un <- rep(NA, N) 
        for(j in 1:N){
          un[j] <- max(abs((x_sample[j,] - x_mean)/sqrt(Sigma)))
        }
        qn[i] <- quantile(un, 0.95)
      }
    
    return(list(betaHat = betaHat, betaHat.var = betaHat.var, qn = qn))
    
    }
  }else{
    
    return(list(betaHat = betaHat))
    
  }
}


setwd("~/Documents/CityU/research/Project on Shanghai Actigraph Data")
load("example.RData")
library(dplyr) ## organize lapply results
library(parallel) ## mcapply
library(ggplot2)
library(ggpubr)
library(fdapace) ## PACE

data <- acc.all
n <- length(unique(data$ID))
Y <- lapply(1:n, function(i) data$Y[((i-1)*7+1):(i*7),])
T <- ncol(Y[[1]]) ## number of observations on the functional domain
S <- nrow(Y[[1]])

Ysm <- list()
for (i in 1:n) {
  ori  <- as.vector(t(Y[[i]]))
  sm <- stats::filter(ori, sides=2, filter = (rep(1, 30)/30)) ## smooth accelerometer dat
  #lowess(ori, f=.08)$y
  sm[is.na(sm)] <- 0
  Ysm[[i]] <- t(array(sm, dim=c(T,S)))
}
 
plot(Y[[4]][1,], type="l")
lines(Ysm[[4]][1,], col="red")

library(scales)
plot.subject <- data.frame(y = as.vector(t(Ysm[[5]])), t=rep(1:T, 7), day = rep(1:7, each=T)) 

p1 <- ggplot(plot.subject, aes(1:(7*T), y)) + 
  geom_line() + 
  theme_bw() + 
  labs(y = "Accelerometer", x = "Time of Day") +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14),) + 
  scale_x_continuous(breaks = seq(T/2, 7*T-T/2, length.out = 7), labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) 

p2 <- ggplot(plot.subject, aes(t, day, fill= y)) + 
  geom_tile() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) + 
  labs(y = "", x = "Time of Day") +
  scale_y_continuous(breaks=1:7, labels=c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) +
  scale_x_continuous(breaks = seq(1, T, length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) + 
  scale_fill_viridis_c(name="Accelerometer",values=c(0,0.1,0.55,1))
ggarrange(p1, p2, ncol=2)


##########################################################################################
## Step 1
##########################################################################################
S.argvals <- T.argvals <- NULL
family <- "gaussian" 

formula <- as.formula("Y ~ grade + gender + income + edu + depression")
model_formula <- as.character(formula)
stopifnot(model_formula[1] == "~" & length(model_formula) == 3)
vars <- unlist(strsplit(model_formula[3], split = ".\\+."))

if(is.null(T.argvals)) T.argvals <- 1:T
if(is.null(S.argvals)) S.argvals <- 1:S

## function "unim" fit univariate model at t and s
unim <- function(s,t){
  
  Yst <- sapply(Y, function(x) x[s,t])
  Xst <- data[seq(s, (n-1)*S + s, by=S), vars]
  dat <- cbind(Xst, Yst)

  if(family == "gaussian"){
    fit_uni <- suppressMessages(lm(formula = as.formula(paste0("Yst ~ ", model_formula[3])), data=dat))
  }else{
    fit_uni <- suppressMessages(glm(formula = as.formula(paste0("Yst ~ ", model_formula[3])), data = dat, family = family))
  }
  betaTilde <- coef(fit_uni)
  designmat <- model.matrix(fit_uni)
  Covst <- summary(fit_uni)$sigma^2
  
  return(list(betaTilde = betaTilde, designmat = designmat, cov = Covst))
}
#massive univariate model
parallel = TRUE
massm <- list() 
if(parallel == TRUE){
  for(i in S.argvals) massm[[i]] <- mclapply(T.argvals, unim, s=i, mc.cores = detectCores() - 1)
}else{
  for(i in S.argvals) massm[[i]] <- lapply(T.argvals, unim, s=i)
}

betaTilde <- lapply(S.argvals, function(i) lapply(massm[[i]], '[[', 1) %>% bind_rows()) %>% bind_rows()
designmat <- lapply(S.argvals, function(i) lapply(massm[[i]], '[[', 2)) 
RTilde <- lapply(S.argvals, function(i) lapply(massm[[i]], '[[', 3) %>% unlist()) %>% do.call("cbind",.)

###############Smoothing
library(mgcv)
library(refund)
smooth="fbps"
betaHat <- array(0, dim=c(T*S, ncol(betaTilde)))
RHat <- array(0, dim=c(T,S))

### smooth raw covariance estimates (estimates to measurement error)
if (smooth=="te"){
  sm <- data.frame(y = as.vector(RTilde), t=rep(1:T, S), day = rep(1:S, each=T)) 
  RHat <- gam(y ~ te(t, day, k = c(15,5)), data = sm)$fitted.values %>% array(dim=c(T,S))
} else if (smooth=="fbps"){
  RHat <- fbps(RTilde, knots = c(15,5))$Yhat 
}
RHat[which(RHat < 0)] <- 0

heatmap(RHat, Rowv = NA, Colv = 'Rowv')

for (p in 1:ncol(betaTilde)){
  if (smooth=="te"){
    sm <- data.frame(y = betaTilde[[p]], t=rep(1:T, S), day = rep(1:S, each=T)) 
    betaHat[,p] <- gam(y ~ te(t, day, k = c(15,5)), data = sm)$fitted.values
  } else if (smooth=="fbps"){
    sm <- array(betaTilde[[p]], dim=c(T,S)) 
    betaHat[,p] <- fbps(sm, knots = c(25,3))$Yhat %>% as.vector()
  }
}

###############
#pointwise beta
pic <- list()
for (p in 1:length(vars)){
  plot.subject <- data.frame(y = betaTilde[[p+1]], t=rep(1:T, 7), day = rep(1:7, each=T)) 
  pic[[p]] <- ggplot(plot.subject, aes(t, day, fill= y)) + 
    geom_tile() + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 16), 
          legend.title = element_text(size = 14), legend.key.width = unit(1, 'cm'), 
          plot.title = element_text(size = 19)) + 
    labs(y = "S", x = "T", title = vars[p]) +
    scale_y_continuous(breaks=1:7, labels=c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) +
    scale_x_continuous(breaks = seq(1, T, length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) +
    scale_fill_gradient2(low="blue", mid="grey98",
                         high="red", space ="Lab", name=expression(paste(hat(beta), "(s,t)")))
  
    
}
ggarrange(pic[[1]],pic[[2]],pic[[3]], 
          pic[[4]],pic[[5]], ncol=3, nrow=2)
#smoothed beta
pic <- list()
for (p in 1:length(vars)){
  plot.subject <- data.frame(y = betaHat[,p+1], t=rep(1:T, 7), day = rep(1:7, each=T)) 
  pic[[p]] <- ggplot(plot.subject, aes(t, day, fill= y)) + 
    geom_tile() + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 16), 
          legend.title = element_text(size = 14), legend.key.width = unit(1, 'cm'), 
          plot.title = element_text(size = 19)) + 
    labs(y = "S", x = "T", title = vars[p]) +
    scale_y_continuous(breaks=1:7, labels=c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) +
    scale_x_continuous(breaks = seq(1, T, length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) +
    scale_fill_gradient2(low="blue", mid="grey98", #limits=c(min(betaHat[,-1]),max(betaHat[,-1])),
                         high="red", space ="Lab", name=expression(paste(hat(beta), "(s,t)")))
}
ggarrange(pic[[1]],pic[[2]],pic[[3]], 
          pic[[4]],pic[[5]], ncol=3, nrow=2)
###############

betaHat <- lapply(1:ncol(betaHat), function(p) array(betaHat[,p], dim =  c(T, S)))

##########################################################################################
## Step 2.1 Centering
##########################################################################################
## Get eta to do FPCA
get_eta <- function(s,t){
  Yst <- sapply(Y, function(x) x[s,t])
  betaHat.st <- sapply(1:(length(vars)+1), function(p) betaHat[[p]][t,s]) %>% array(dim = c(length(vars)+1, 1))
  if(family == "gaussian"){
    Ytilde <- designmat[[s]][[t]] %*% betaHat.st
    eta <- Yst - Ytilde
  }
  return(eta = eta)
}
#massive eta
massEta <- list()
if(parallel == TRUE){
  for(i in S.argvals) massEta[[i]] <- mclapply(T.argvals, get_eta, s=i, mc.cores = detectCores() - 1)
}else{
  for(i in S.argvals) massEta[[i]] <- lapply(T.argvals, get_eta, s=i)
}

eta <- NULL
for(i in 1:n) {
 eta.i <- lapply(S.argvals, function(s) lapply(T.argvals, function(t) massEta[[s]][[t]][i])) %>% unlist() #ith sample
 eta <- rbind(eta, array(eta.i, dim=c(T,S)))
} 

##########################################################################################
## Step 2.2 Computing marginal eigenfunctions psi
##########################################################################################
pc.s <- 3
dataType <- 'DenseWithMV' #'Sparse'
fpca.op <- list(dataType = dataType, methodSelectK = 'FVE')

sList <- list()
etaList <- list()

for (nt in 1:(n*T)){
  ind <- which(!is.na(eta[nt,]))
  etaList[[nt]] <- eta[nt,ind]
  sList[[nt]] <- S.argvals[ind]
}

res.psi <- suppressWarnings(FPCA(etaList, sList, optns = fpca.op))
psi <- res.psi$phi #psi(s)_hat
xiEst <- res.psi$xiEst #xi(t)_hat

par(mfrow=c(1,3))
plot(1:7, psi[,1], type="l", xlab="Day") #res.psi$cumFVE is cumulative variance
plot(1:7, psi[,2], type="l", xlab="Day")
plot(1:7, psi[,3], type="l", xlab="Day")

##########################################################################################
## Step 2.3 Computing B-spline phi
##########################################################################################
nbasis <- 30
res.b <- list()

Bt <- bs(T.argvals, df = nbasis)
bspl.fit <- function(i){
  fit <- lm(bspl[i,] ~ 0 + Bt)
  b <- coef(fit)
  return(b)
}

for (j in 1:pc.s){
  bspl <- matrix(xiEst[,j], nrow=n, ncol=T, byrow=TRUE)
  res.b[[j]] <- lapply(1:n, function(i) bspl.fit(i)) %>% bind_rows() %>% as.matrix()
}

par(mfrow=c(1,1))
plot(xiEst[1:T,1], type="l")
lines(Bt%*% res.b[[1]][1,], col="red")


kernel <- function(i,s,u){
  ker <- array(0, dim=c(ncol(Bt), ncol(Bt)))
  for (j in 1:pc.s){
    for (h in 1:pc.s){
      ker <- ker + psi[s,j]*psi[u,h]*cov(rbind(res.b[[j]][i,], res.b[[h]][i,]))
    }
  }
  edcomp <- eigen(ker)
  eigen.positive <- which(edcomp$values > 0)
  if(length(eigen.positive) == ncol(Bt)){
    ker.trimmd <- ker
  }else{
    ker.trimmed <- edcomp$vectors[, eigen.positive] %*% diag(edcomp$values[eigen.positive]) %*% t(edcomp$vectors[,eigen.positive])
  }
  ker.trimmed
}

covBeta <- function(s,u,t,v){ #set s,u,t,v
  cBeta <- 0
  cEta <- array(0, dim=c(n,n))
  for (i in 1:n){
    ker <- kernel(i, s, u)
    cEta[i,i] <- t(Bt[t,]) %*% ker %*% Bt[v,] + RHat[t,s]
    print( t(Bt[t,]) %*% ker %*% Bt[v,]) ; print(RHat[t,s])
  }  
  xst1 <- designmat[[s]][[t]]
  xuv2 <- designmat[[u]][[v]]
  cBeta <- solve(t(xst1) %*% xst1) %*% t(xst1) %*% cEta %*% xuv2 %*% solve(t(xuv2) %*% xuv2) 
  cBeta
}

###raw estimate of Var(beta(s,t))
var.beta.ts.tilde <- array(0, dim = c(length(betaHat),length(betaHat),S,T))
for(s in S.argvals){
  for (t in T.argvals){
    print(t)
    u <- s; v <- t
    var.beta.ts.tilde[,,s,t] <- covBeta(s,u,t,v)
  }
}


###final estimate of Var(beta(s,t))
var.beta.ts.hat <- array(0, dim = c(T,S,length(betaHat)))
for(p in 1:length(betaHat)){
  for(s in 1:S){
    for(t in 1:T){
        print(var.beta.ts.tilde[p,p,s,t])
        var.beta.ts.hat[t,s,p] <- var.beta.ts.tilde[p,p,s,t]
    }
  }
}

var.beta.ts.hatsm <- array(0, dim = c(T*S,length(betaHat)))
for (p in 1:length(betaTilde)){
  if (smooth=="te"){
    sm <- data.frame(y = var.beta.ts.hat[,,p], t=rep(1:T, S), day = rep(1:S, each=T)) 
    var.beta.ts.hatsm[,p] <- gam(y ~ te(t, day, k = c(15,5)), data = sm)$fitted.values
  } else if (smooth=="fbps"){
    var.beta.ts.hatsm[,p] <- fbps(var.beta.ts.hat[,,p], knots = c(15,3))$Yhat %>% as.vector()
  }
}


pic <- list()
for (p in 1:length(vars)){
  
  betaHat.up <- as.vector(betaHat[[p+1]]) + 2*sqrt(as.vector(var.beta.ts.hat[,,p+1]))
  betaHat.low <- as.vector(betaHat[[p+1]]) - 2*sqrt(as.vector(var.beta.ts.hat[,,p+1]))
  
  plot.subject <- data.frame(betaHat.low=betaHat.low, betaHat.up=betaHat.up, t=rep(1:T, 7), day = rep(1:S, each=T)) 
  plot.subject$sign <- ifelse(plot.subject$betaHat.up < 0,  plot.subject$betaHat.up, 
                              ifelse(plot.subject$betaHat.low > 0, plot.subject$betaHat.low, 0))
  
  pic[[p]] <- ggplot(plot.subject, aes(t, day, fill= sign)) + 
    geom_tile() + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 16), 
          legend.title = element_text(size = 14), legend.key.width = unit(1, 'cm'), 
          plot.title = element_text(size = 19)) + 
    labs(y = "S", x = "T", title = vars[p]) +
    scale_y_continuous(breaks=1:7, labels=c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) +
    scale_x_continuous(breaks = seq(1, T, length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) +
    scale_fill_gradient2(low="blue", mid="grey98", 
                         high="red", space ="Lab", name=expression(paste(hat(beta), "(s,t)")))

}
ggarrange(pic[[1]],pic[[2]],pic[[3]], 
          pic[[4]],pic[[5]], ncol=3, nrow=2)


pic <- list()
for (p in 1:length(vars)){
  
  betaHat.up <- as.vector(betaHat[[p+1]]) + 2*sqrt(as.vector(RHat))
  betaHat.low <- as.vector(betaHat[[p+1]]) - 2*sqrt(as.vector(RHat))
  
  plot.subject <- data.frame(betaHat.low=betaHat.low, betaHat.up=betaHat.up, t=rep(1:T, 7), day = rep(1:S, each=T)) 
  plot.subject$sign <- ifelse(plot.subject$betaHat.up < 0,  plot.subject$betaHat.up, 
                              ifelse(plot.subject$betaHat.low > 0, plot.subject$betaHat.low, 0))
  
  pic[[p]] <- ggplot(plot.subject, aes(t, day, fill= sign)) + 
    geom_tile() + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 16), 
          legend.title = element_text(size = 14), legend.key.width = unit(1, 'cm'), 
          plot.title = element_text(size = 19)) + 
    labs(y = "S", x = "T", title = vars[p]) +
    scale_y_continuous(breaks=1:7, labels=c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) +
    scale_x_continuous(breaks = seq(1, T, length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) +
    scale_fill_gradient2(low="blue", mid="grey98", 
                         high="red", space ="Lab", name=expression(paste(hat(beta), "(s,t)")))
  
}
ggarrange(pic[[1]],pic[[2]],pic[[3]], 
          pic[[4]],pic[[5]], ncol=3, nrow=2)
