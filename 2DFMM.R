#' Fit longitudinal 2D functional mixed model regression using the proposed 3-step approach
#' 
#' @param formula two-sided formula object in lm() format, except that the response is a matrix  
#' @param data dataframe containing variables in formula; 
#' dataframe contains covariates (vec or data.frame) and response (data.frame)
#' @param S number of visits
#' method="OLS" OLS estimator, or otherwise, Ridge estimator is only applied when there are at least two covariates of interest
#' @param smoother sandwich smoother (sandwich) or tensor product smoother (te)
#' @param knots=c(Sknots,Tknots) the number of knots for S and T
#' @param fpca.opt A list of options control parameters specified by list, dataType: dense and regular (Dense); Very densely and regularly 
#' observed data: empirical mean and Densely recorded but irregular design, or contaminated with error: pre-smoothing for individual 
#' curves (DenseWithMV); Sparse random design (Sparse)
#' @param parallel whether to run parallel computing (True only for Linux/Mac users)
#' @export
#' @return a list containing estimated beta(s,t), covariance matrix of beta(s,t)


fmm2d <- function(formula, data, S, smoother = "sandwich", knots = NULL, fpca.opt = list(dataType = 'DenseWithMV', methodSelectK = 'FVE'), parallel = FALSE)
{
  
  library(dplyr) ## organize lapply results
  library(parallel) ## mcapply
  #library(glmnet) ##ridge estimator
  library(mgcv) ## tensor product smoothing
  library(refund) ## sandwich smoothing
  library(fdapace) ## PACE
  library(Matrix)

  model_formula <- as.character(formula)
  stopifnot(model_formula[1] == "~" & length(model_formula) == 3)
  vars <- unlist(strsplit(model_formula[3], split = ".\\+."))
  
  n <- nrow(data[,model_formula[2]])/S
  T <- ncol(data[,model_formula[2]])
  T.argvals <- 1:T
  S.argvals <- 1:S
  Y <- lapply(1:n, function(i) data[,model_formula[2]][((i-1)*S+1):(i*S),])
  var.mat <- NULL
  for(v in vars){
    if(is.array(data[,v])) var.mat <- c(v,var.mat) #determine which covariate is bivariate function
  } 
  
  ###############Presmoothing  
  #Ysm <- list()
  #for (i in 1:n) {
  #  ori  <- as.vector(t(Y[[i]]))
  #  sm <- stats::filter(ori, sides=2, filter = (rep(1, 10)/10)) ## smooth accelerometer dat  #lowess(ori, f=.08)$y
  #  sm[is.na(sm)] <- ori[is.na(sm)]
  #  Ysm[[i]] <- t(array(sm, dim=c(T,S)))
  #}
  #Y <- Ysm

  ##########################################################################################
  ## Step 1 Bivariate pointwise estimate
  ##########################################################################################
  print("Step 1: Bivariate pointwise estimate")
  #function "bim" fit bivariate model at t and s
  bim <- function(s,t){
    
    Yst <- sapply(Y, function(x) x[s,t])
    Xst <- data[seq(s, (n-1)*S + s, by=S), vars]
    
    if(length(vars) == 1){ #if there is only one covariate
      if(!is.null(var.mat)){ #if the covariate is bivariate function
        Xst <- Xst[,t]
        Xst <- array(Xst, dim=c(n, 1)) 
        colnames(Xst) <- vars
      }else{
        Xst <- array(Xst, dim=c(n, 1)) 
        colnames(Xst) <- vars
      }
    } else{ #if more than one covariate
      if(!is.null(var.mat)){ #if at least one covariates have bivariate functions
        Xst.sub <- as.matrix(Xst[,var.mat])[,seq(t, (length(var.mat)-1)*T + t, by=T)] #subset covariates which are bivariate functions
        colnames(Xst.sub) <- var.mat
        Xst <- cbind(Xst[,!names(Xst) %in% var.mat], Xst.sub)
        Xst <- Xst[, vars]
      }
    }
    dat <- data.frame(cbind(Xst, Yst))
      
    fit_bi <- suppressMessages(lm(formula = as.formula(paste0("Yst ~ ", model_formula[3])), data=dat))
    betaTilde <- coef(fit_bi)
    designmat <- model.matrix(fit_bi)
    Covst <- summary(fit_bi)$sigma^2
    
    #if(method == "Ridge"){
    #  fit_bi <- suppressMessages(cv.glmnet(x=matrix(dat[,model_formula[3]]), y=dat[,"Yst"], alpha=0))
    #  betaTilde <- as.numeric(coef(fit_bi))
    #  designmat <- model.matrix(~ dat[,model_formula[3]])
    #  pred <- predict(fit_bi, newx = dat[,model_formula[3]])
    #  Covst <- sum((pred - dat[,"Yst"])^2)/(n-length(vars)-1)
    #}
    
    return(list(betaTilde = betaTilde, designmat = designmat, cov = Covst))
  }
  
  #massive bivariate model
  massm <- list() 
  if(parallel == TRUE){
    for(i in S.argvals) massm[[i]] <- mclapply(T.argvals, bim, s=i, mc.cores = detectCores() - 1)
  }else{
    for(i in S.argvals) massm[[i]] <- lapply(T.argvals, bim, s=i)
  }
  
  betaTilde <- lapply(S.argvals, function(i) lapply(massm[[i]], '[[', 1) %>% bind_rows()) %>% bind_rows()
  designmat <- lapply(S.argvals, function(i) lapply(massm[[i]], '[[', 2)) 
  RTilde <- lapply(S.argvals, function(i) lapply(massm[[i]], '[[', 3) %>% unlist()) %>% do.call("cbind",.)
  
  ##########################################################################################
  ## Step 2 Smoothing
  ##########################################################################################
  print("Step 2: Bivariate smoothing")
  betaHat <- array(0, dim=c(T*S, ncol(betaTilde))) 
  RHat <- array(0, dim=c(T,S)) #R is the transponse of R from the paper
  lambda <- array(0, dim=c(2, ncol(betaTilde)))
  S2kS1 <- S2 <- S1 <- list()
  
  ### smooth raw covariance estimates (estimates to measurement error)
  if(is.null(knots)) knots <- c(min(round(S/1.5), 35), min(round(T/4), 35))
  Sknots <- knots[1]
  Tknots <- knots[2]
  
  #obtain S1 and S2 from Sandwich smoother
  sL <- refund:::pspline.setting(x=S.argvals, knots=refund:::select_knots(S.argvals,Sknots), p=3, m=2, periodicity=FALSE, weight=NULL)
  tL <- refund:::pspline.setting(x=T.argvals, knots=refund:::select_knots(T.argvals,Tknots), p=3, m=2, periodicity=FALSE, weight=NULL)
  
  B1 <- tL$B
  Bt1 <- Matrix(t(as.matrix(B1))) #B1^T
  P1 <- tL$P #D1^T %*% D1
  
  B2 <- sL$B 
  Bt2 <- Matrix(t(as.matrix(B2))) #B2^T
  P2 <- sL$P #D1^T %*% D1
  
  #Sandwich/tensor product smoother 
  for (p in 1:ncol(betaTilde)){ 
    sm <- array(betaTilde[[p]], dim=c(T,S)) 
    y <- fbps(sm, knots = c(Tknots,Sknots))
    if(smoother == "te"){
      te.sm <- data.frame(y = betaTilde[[p]], t=rep(1:T, S), s = rep(1:S, each=T)) 
      betaHat[,p] <- gam(y ~ te(t, s, k = c(Tknots,Sknots)), data = te.sm)$fitted.values
    }else if(smoother == "sandwich"){
      betaHat[,p] <- y$Yhat %>% as.vector()
    }else{ 
      stop("Smoother must be either sanwich smoother or tensor product smooths!")
    }
    lambda[,p] <- y$lambda
    S1[[p]] <- B1 %*% solve(Bt1 %*% B1 + lambda[1,p]*P1) %*% Bt1
    S2[[p]] <- B2 %*% solve(Bt2 %*% B2 + lambda[2,p]*P2) %*% Bt2
    S2kS1[[p]] <- kronecker(S2[[p]], S1[[p]]) #S2 \otimes S1
  }
  betaHat <- lapply(1:ncol(betaHat), function(p) array(betaHat[,p], dim =  c(T, S))) #betaHat is the transpose of betaHat from paper
  names(betaHat) <- c("intercept", vars)
  colnames(lambda) <- c("intercept", vars)
  RHat <- fbps(RTilde, knots = c(Tknots,Sknots))$Yhat 
  RHat[which(RHat < 0)] <- 0
  
  rm(betaTilde, RTilde)

  ##########################################################################################
  ## Step 3.1 Centering
  ##########################################################################################
  print("Step 3: Computing variance of bivariate estimates")
  
  ## Get eta to do FPCA
  get_eta <- function(s,t){
    Yst <- sapply(Y, function(x) x[s,t])
    betaHat.st <- sapply(1:(length(vars)+1), function(p) betaHat[[p]][t,s]) %>% array(dim = c(length(vars)+1, 1))
    Ytilde <- designmat[[s]][[t]] %*% betaHat.st
    eta <- Yst - Ytilde
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
  ## Step 3.2 Computing marginal covariance w.r.t. S and eigenfunctions psi
  ##########################################################################################
  sList <- list()
  etaList <- list()
  
  for (nt in 1:(n*T)){
    ind <- which(!is.na(eta[nt,]))
    etaList[[nt]] <- eta[nt,ind]
    sList[[nt]] <- S.argvals[ind]
  }
  
  fpca <- suppressWarnings(FPCA(etaList, sList, optns = fpca.opt)) #fpca$cumFVE is cumulative variance
  psi <- fpca$phi 
  xiEst <- fpca$xiEst 
  #tau <- fpca$lambda
  pc.s <- ncol(xiEst)
  s.cov <- array(0, dim=c(pc.s, S, S))
  for(j in 1:pc.s){
      s.cov[j,,] <- psi[, j] %*% t(psi[,j]) #fpca$fittedCov
  }
  
  ##########################################################################################
  ## Step 3.3 Computing B-spline smoothing xi(t) and marginal covariance w.r.t. T 
  ##########################################################################################
  Bt <- as.matrix(B1)
  b <- array(0, dim=c(n,ncol(Bt),pc.s))
  
  bspl.fit <- function(i){
    fit <- lm(bspl[i,] ~ 0 + Bt)
    b <- coef(fit)
    return(b)
  }
  
  for (j in 1:pc.s){
    bspl <- matrix(xiEst[,j], nrow=n, ncol=T, byrow=TRUE)
    bspl.result <- lapply(1:n, function(i) bspl.fit(i))
    b[,,j] <- bspl.result %>% bind_rows() %>% as.matrix()
  }
     
  #store marginal covariance of T 
  t.cov <- array(0, dim=c(pc.s, T, T))
  for(j in 1:pc.s){
    ker <- array(0, dim=c(ncol(Bt), ncol(Bt)))
    for (i in 1:n){
      ker <- ker + kronecker(t(b[i,,j]), b[i,,j])
    }
    t.cov[j,,] <- (1/n) * Bt %*% ker %*% t(Bt)
  }

  ##########################################################################################
  ## Step 3.4 Computing four-dim covariance operator w.r.t. beta(s,t) 
  ##########################################################################################
    
  covBeta <- function(s,u,t,v){ #set s,u,t,v
    if(s==u & t==v){
      cEta <- RHat[t,s] 
    }else{
      cEta <- Reduce("+", lapply(1:pc.s,function(j) s.cov[j,s,u]*t.cov[j,t,v]))
    }
    xst1 <- designmat[[s]][[t]]
    xuv2 <- designmat[[u]][[v]]
    cBeta <- cEta * solve(t(xst1) %*% xst1) %*% t(xst1) %*% xuv2 %*% solve(t(xuv2) %*% xuv2) 
    cBeta
  }
  
  ## raw estimate of Cov(beta(s,t), beta(s,t))
  cov.beta.ts.tilde <-  array(0, dim = c(S*T, S*T, length(betaHat)))
  for(s in 1:S){
    for(u in s:S){
      for(t in 1:T){
        for(v in t:T){
          tmp <- covBeta(s,u,t,v)
          for(p in 1:length(betaHat)){
            cov.beta.ts.tilde[(u-1)*T+v, (s-1)*T+t, p] <- 
              cov.beta.ts.tilde[(u-1)*T+t, (s-1)*T+v, p] <- 
              cov.beta.ts.tilde[(s-1)*T+v, (u-1)*T+t, p] <- 
              cov.beta.ts.tilde[(s-1)*T+t, (u-1)*T+v, p] <- tmp[p,p]
          }
        }
      }
    }
  }
  
  ## refined smoothing estimates covariance 
  cov.beta.ts.hat <- array(0, dim = c(T*S, T*S, length(betaHat)))
  for(p in 1:length(betaHat)){
    edcomp <- eigen(cov.beta.ts.tilde[,,p])
    eigen.positive <- which(edcomp$values > 0)
    if(length(eigen.positive) == S*T){
      ker.trimmed <- cov.beta.ts.tilde[,,p]
    }else{
      ker.trimmed <- edcomp$vectors[,eigen.positive] %*% diag(edcomp$values[eigen.positive]) %*% t(edcomp$vectors[,eigen.positive])
    }
    cov.beta.ts.hat[,,p] <- matrix(S2kS1[[p]] %*% ker.trimmed %*% t(S2kS1[[p]]))
  }
  
  rm(cov.beta.ts.tilde)

  return(list(betaHat = betaHat, betaHat.cov = cov.beta.ts.hat))
}
