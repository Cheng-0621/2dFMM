#' An individual file to obtain SCB by bootstrap method. It is decoupled with PCB and gives more stable way. 

bootstrap <- function(formula, data, S, T, knots, betaHat, betaHatcov, B, M, smoother){
  T.argvals <- 1:T
  S.argvals <- 1:S
  n <- nrow(data)/S
  Sknots <- knots[1]
  Tknots <- knots[2]
  sL <- refund:::pspline.setting(x=S.argvals, knots=refund:::select_knots(S.argvals,Sknots), p=3, m=2, periodicity=FALSE, weight=NULL)
  tL <- refund:::pspline.setting(x=T.argvals, knots=refund:::select_knots(T.argvals,Tknots), p=3, m=2, periodicity=FALSE, weight=NULL)
  
  B1 <- tL$B
  B2 <- sL$B
  Bt1 <- Matrix(t(as.matrix(B1))) #B1^T
  P1 <- tL$P #D1^T %*% D1
  
  qn <- rep(0, length = length(betaHat))
  B <- B #bootstrap times
  
  Bs <- as.matrix(B2)
  Bt <- as.matrix(B1)
  betaHat.boot <- NULL
  
  fpca.opt <- list(dataType = 'Dense', methodSelectK = 'FVE')
  parallel <- FALSE
  
  for (b in 1:B){
    print(b)
    sample.ind <- sample(unique(data$ID), size = length(unique(data$ID)), replace = TRUE) #bootstrap id with replacement
    row.ind <- NULL
    for (id in 1:length(sample.ind)){
      row.ind <- c(row.ind, which(data$ID == sample.ind[id]))
    }
    data.boot <- data[row.ind,] #b_th dataset
    fmm2d_boot <- fmm2d(formula, data = data.boot, S, smoother, knots, fpca.opt, parallel,
                        pcb = FALSE, scb = FALSE, silence = TRUE)
    betaHat.boot[[b]] <- fmm2d_boot$betaHat
  }
  #computing marginal decomposition on FPCA for beta.boot
  psi.boot <- list()
  ker.boot <- list()
  betaHat.boot.avg <- list()
  for (p in 1:length(qn)){
    betaHat.boot.p <- lapply(1:B, function(b) betaHat.boot[[b]][[p]])
    betaHat.boot.avg[[p]] <- Reduce("+", betaHat.boot.p)/B
    betaHat.boot.p.demean <- lapply(1:B, function(b) betaHat.boot.p[[b]] - betaHat.boot.avg[[p]]) #demean each betaHat.boot.ps
    eigenS.boot <- s.fpca(eta = do.call("rbind", betaHat.boot.p.demean), S = S, fpca.opt = fpca.opt)
    psi.boot[[p]] <- eigenS.boot$psi
    bsplineT.boot <- b.lm(Bt = Bt, pc.s = eigenS.boot$pc.s, xiEst = eigenS.boot$xiEst) 
    ker.boot[[p]] <- bsplineT.boot$ker
  }
  rm(betaHat.boot.p, betaHat.boot.p.demean)
  
  M <- M
  for(p in 1:length(qn)){
    sp <- sqrt(diag(betaHatcov[,,p]))
    qm <- lapply(1:M, function(m) array(0, dim=c(T,S)))
    JB <- dim(ker.boot[[p]])[1]
    for(j in 1:JB){
      Sig <- (1/n)*ker.boot[[p]][j,,]  
      Sig.r <- Sig[-c(1,ncol(Bt)), -c(1,ncol(Bt))] #remove tail problems
      Bt.r <- Bt[,-c(1,ncol(Bt))] #remove tail problems
      umat <- rmvnorm(M, mean = rep(0, ncol(Bt.r)), sigma = Sig.r)
      Btu <- apply(umat, 1, function(x) Bt.r %*% x)  
      qm_j <- lapply(1:M, function(m) kronecker(Btu[,m], t(psi.boot[[p]][,j])))
      qm <- mapply("+", qm_j, qm, SIMPLIFY = FALSE)
      betaHat.boot.p <- lapply(qm, function(q) q + betaHat.boot.avg[[p]]) 
    }
    qmdvar <- lapply(1:M, function(m) abs(betaHat.boot.p[[M]] - betaHat[[p]])/array(sp, dim=c(T,S)))
    qmstar <- lapply(1:M, function(m) max(qmdvar[[m]]))
    qn[p] <- quantile(unlist(qmstar), 0.80) 
    
  }
  return(qn)
}