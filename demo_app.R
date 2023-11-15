rm(list = ls())

load("example.RData")
source("2DFMM.R")

library(dplyr) ## organize lapply results
library(ggplot2)
library(ggpubr)

data <- phyact
N <- length(unique(data$ID))
Y <- lapply(1:N, function(i) data$Y[((i-1)*7+1):(i*7),])
T <- ncol(Y[[1]]) #the number of functional grids
S <- nrow(Y[[1]]) #the number of longitudinal grids


#####################
###baseline model#### 
#####################
fit_baseline <- fmm2d(formula=Y ~ gender + grade + income + edu, data=data, S=7, smoother="te", 
                   knots=c(4, 35), fpca.opt = list(dataType = 'Dense', methodSelectK = 'FVE'),  parallel = TRUE)


vars <- c("intercept", "gender", "grade", "income", "edu")
betaHat <- fit_baseline$betaHat
betaHatcov <- fit_baseline$betaHat.cov

pic <- list()
for (p in 1:length(betaHat)){
  
  betaHat.up <- as.vector(betaHat[[p]]) + 2*sqrt(diag(betaHatcov[,,p]))
  betaHat.low <- as.vector(betaHat[[p]]) - 2*sqrt(diag(betaHatcov[,,p]))
  
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
    scale_fill_gradient2(low="Steel Blue 3", mid="grey98", 
                         high="red", space ="Lab", name=expression(paste(hat(beta), "(s,t)")))
  
}
ggarrange(pic[[1]],pic[[2]],pic[[3]], pic[[4]],pic[[5]], ncol=3, nrow=2)


#####################
###baseline+var model
#####################


pic <- list()
vars <-  c("stress", "anxiety", "depression","acips","erq_cr", "erq_es", "sf12v2_mcs", "sf12v2_pcs")
for (p in 1:length(vars)){
  
  formula <- paste("Y ~ gender + grade + income + edu + ", vars[p])
  fit_fmm2d <- fmm2d(formula=as.formula(formula), data=data, S=S, smoother="te", 
                     knots=c(4, 35), fpca.opt = list(dataType = 'Dense', methodSelectK = 'FVE'),
                     parallel = FALSE)
  
  varBetaEst <- fit_fmm2d$betaHat[[6]]
  varBetaEstcov <- fit_fmm2d$betaHat.cov[,,6]
  
  betaHat.up <- as.vector(varBetaEst) + 2*sqrt(diag(varBetaEstcov))
  betaHat.low <- as.vector(varBetaEst) - 2*sqrt(diag(varBetaEstcov))
  
  plot.subject <- data.frame(betaHat.low=betaHat.low, betaHat.up=betaHat.up, t=rep(1:T, S), day = rep(1:S, each=T)) 
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
    scale_fill_gradient2(low="Steel Blue 3", mid="grey98", 
                         high="red", space ="Lab", name=expression(paste(hat(beta), "(s,t)")))
}
ggarrange(pic[[1]],pic[[2]],pic[[3]],pic[[4]],
          pic[[5]],pic[[6]],pic[[7]],pic[[8]], ncol=4, nrow=2)

