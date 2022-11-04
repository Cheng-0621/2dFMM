#' Fit longitudinal 2D functional mixed model regression using the proposed 3-step approach
#' 
#' @param formula two-sided formula object in lm() format, except that the response is a matrix  
#' @param data data frame containing variables in formula
#'@param S number of visits
#' @param family GLM family of the response
#' @param smoother tensor product ("te")/sandwich smoothing ("fbps")
#' @param knots=c(Sknots,Tknots) the number of knots for S and T
#' @param fpca.opt A list of options control parameters specified by list
#' @param parallel whether to run parallel computing
#' @export
#' @return a list containing estimated beta(s)


fmm2d <- function(formula, data, S, family = "gaussian", smoother, knots = NULL, fpca.opt = list(dataType = 'DenseWithMV', methodSelectK = 'FVE'), parallel = TRUE){
  
  library(dplyr) ## organize lapply results
  library(parallel) ## mcapply
  library(mgcv) ## tensor product smoothing
  library(refund) ## sandwich smoothing
  library(fdapace) ## PACE

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
  Ysm <- list()
  for (i in 1:n) {
    ori  <- as.vector(t(Y[[i]]))
    sm <- stats::filter(ori, sides=2, filter = (rep(1, 10)/10)) ## smooth accelerometer dat  #lowess(ori, f=.08)$y
    sm[is.na(sm)] <- ori[is.na(sm)]
    Ysm[[i]] <- t(array(sm, dim=c(T,S)))
  }

  ##########################################################################################
  ## Step 1 Bivariate pointwise estimate
  ##########################################################################################
  print("Step 1: Bivariate pointwise estimate")
  #function "bim" fit bivariate model at t and s
  bim <- function(s,t){
    
    Yst <- sapply(Ysm, function(x) x[s,t])
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
      
    
    if(family == "gaussian"){
      fit_bi <- suppressMessages(lm(formula = as.formula(paste0("Yst ~ ", model_formula[3])), data=dat))
    }else{
      fit_bi <- suppressMessages(glm(formula = as.formula(paste0("Yst ~ ", model_formula[3])), data = dat, family = family))
    }
    betaTilde <- coef(fit_bi)
    designmat <- model.matrix(fit_bi)
    Covst <- summary(fit_bi)$sigma^2
    
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
  print("Step 2: Smoothing")
  betaHat <- array(0, dim=c(T*S, ncol(betaTilde)))
  RHat <- array(0, dim=c(T,S))
  
  ### smooth raw covariance estimates (estimates to measurement error)
  if(is.null(knots)) knots <- c(min(round(S/1.5), 35), min(round(T/4), 35))
  Sknots <- knots[1]
  Tknots <- knots[2]

  if (smoother=="te"){
    sm <- data.frame(y = as.vector(RTilde), t=rep(1:T, S), s = rep(1:S, each=T)) 
    RHat <- gam(y ~ te(t, s, k = c(Tknots,Sknots)), data = sm)$fitted.values %>% array(dim=c(T,S))
  } else if (smoother=="fbps"){
    RHat <- fbps(RTilde, knots = c(Tknots,Sknots))$Yhat 
  }
  RHat[which(RHat < 0)] <- 0
  
  for (p in 1:ncol(betaTilde)){
    if (smoother=="te"){
      sm <- data.frame(y = betaTilde[[p]], t=rep(1:T, S), s = rep(1:S, each=T)) 
      betaHat[,p] <- gam(y ~ te(t, s, k = c(Tknots,Sknots)), data = sm)$fitted.values
    } else if (smoother=="fbps"){
      sm <- array(betaTilde[[p]], dim=c(T,S)) 
      betaHat[,p] <- fbps(sm, knots = c(Tknots,Sknots))$Yhat %>% as.vector()
    }
  }
  betaHat <- lapply(1:ncol(betaHat), function(p) array(betaHat[,p], dim =  c(T, S)))
  names(betaHat) <- c("intercept", vars)
  
  ##########################################################################################
  ## Step 3.1 Centering
  ##########################################################################################
  print("Step 3: Computing variance of bivariate estimates")
  
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
  ## Step 3.2 Computing marginal eigenfunctions psi
  ##########################################################################################
  sList <- list()
  etaList <- list()
  
  for (nt in 1:(n*T)){
    ind <- which(!is.na(eta[nt,]))
    etaList[[nt]] <- eta[nt,ind]
    sList[[nt]] <- S.argvals[ind]
  }
  
  res.psi <- suppressWarnings(FPCA(etaList, sList, optns = fpca.opt)) #res.psi$cumFVE is cumulative variance
  psi <- res.psi$phi #psi(s)_hat
  xiEst <- res.psi$xiEst #xi(t)_hat
  pc.s <- ncol(xiEst)
  
  ##########################################################################################
  ## Step 3.3 Computing B-spline phi
  ##########################################################################################
  nbasis <- 30
  res.b <- list()
  
  Bt <- bs(T.argvals, df = nbasis)
  bspl.fit <- function(i){
    fit <- lm(bspl[i,] ~ 0 + Bt)
    b <- coef(fit)
    return(b)
  }
  
  for (j in 1:ncol(xiEst)){
    bspl <- matrix(xiEst[,j], nrow=n, ncol=T, byrow=TRUE)
    res.b[[j]] <- lapply(1:n, function(i) bspl.fit(i)) %>% bind_rows() %>% as.matrix()
  }
  
  ##########################################################################################
  ## Step 3.4 Estimating covariance of beta(s,t)
  ##########################################################################################
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
      #print( t(Bt[t,]) %*% ker %*% Bt[v,]) ; print(RHat[t,s])
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
      u <- s; v <- t
      var.beta.ts.tilde[,,s,t] <- covBeta(s,u,t,v)
    }
  }
  
  ###final estimate of Var(beta(s,t))
  var.beta.ts.hat <- array(0, dim = c(T,S,length(betaHat)))
  for(p in 1:length(betaHat)){
    for(s in 1:S){
      for(t in 1:T){
        var.beta.ts.hat[t,s,p] <- var.beta.ts.tilde[p,p,s,t]
      }
    }
  }
  
  ###smoothing final estimate of Var(beta(s,t))
  var.beta.ts.hatsm <- array(0, dim = c(T*S,length(betaHat)))
  for (p in 1:length(betaTilde)){
    if (smoother=="te"){
      sm <- data.frame(y = var.beta.ts.hat[,,p], t=rep(1:T, S), day = rep(1:S, each=T)) 
      var.beta.ts.hatsm[,p] <- gam(y ~ te(t, day, k = c(Tknots,Sknots)), data = sm)$fitted.values
    } else if (smoother=="fbps"){
      var.beta.ts.hatsm[,p] <- fbps(var.beta.ts.hat[,,p], knots = c(Tknots,Sknots))$Yhat %>% as.vector()
    }
  }
  
  return(list(betaHat = betaHat, betaHat.var = var.beta.ts.hatsm))
  
}
