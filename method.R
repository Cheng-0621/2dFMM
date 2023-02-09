rm(list = ls())

setwd("~/Documents/CityU/research/Project on Shanghai Actigraph Data")
load("data/example.RData")
source("codes/2DFMM.R")

library(dplyr) ## organize lapply results
library(parallel) ## mcapply
library(ggplot2)
library(ggpubr)
library(fdapace) ## PACE

data <- phyact
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
plot.subject <- data.frame(y = as.vector(t(Y[[22]])), t=rep(1:T, 7), day = rep(1:7, each=T)) 

p1 <- ggplot(plot.subject, aes(1:(7*T), y)) + 
  geom_line() + 
  theme_bw() + 
  labs(y = "Accelerometer", x = "Time of Day") +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14),) + 
  scale_x_continuous(breaks = seq(T/2, 7*T-T/2, length.out = 7), labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) 

p2 <- ggplot(plot.subject, aes(t, day, fill= is.na(y))) + 
  geom_tile() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) + 
  labs(y = "", x = "Time of Day") +
  scale_y_continuous(breaks=1:7, labels=c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) +
  scale_x_continuous(breaks = seq(1, T, length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) + 
  scale_fill_viridis_c(name="Accelerometer",values=c(0,0.1,0.55,1))
ggarrange(p1, p2, ncol=2)



ptm <- proc.time()
fit_fmm2d <- fmm2d(formula=Y~ gender + grade + income + edu +
                   stress, data=data, S=7, method = "OLS", knots=c(4, 35),
                   fpca.opt = list(dataType = 'Dense', methodSelectK = 'FVE'),  parallel = TRUE)
time_fmm2d <- (proc.time() - ptm)[3]

vars <- c("intercept", "gender", "grade", "income", "edu", 
          "stress", "anxiety", "depression",
          "acips","erq_cr","erq_es", "sf12v2_mcs", "sf12v2_pcs" )
betaHat <- fit_fmm2d$betaHat
betaHatcov <- fit_fmm2d$betaHat.cov


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
ggarrange(pic[[1]],pic[[2]],pic[[3]], 
          pic[[4]],pic[[5]],pic[[6]],
          ncol=3, nrow=2)


pic.var[[5]] <- pic[[6]]
ggarrange(pic.var[[1]],pic.var[[2]], pic.var[[3]], pic.var[[4]],
          pic.var[[5]],
          ncol=3, nrow=2)
                

#section 3.5
betaHat.sub <- lapply(1:3, function(i) array(0, dim=c(200,7)))
betaHat.cov.sub <- lapply(1:3, function(i) rep(0, 200*7))
for(c in unique(data$class)){
  print(c)
  subdata <- data[data$class == c,]
  if(length(unique(subdata$gender)) == 2 &  length(unique(subdata$stress)) >= 2 ){
    fit_fmm2d <- fmm2d(formula= Y~ gender +
                         stress, data=subdata, S=7, method = "OLS", knots=c(4, 35),
                       fpca.opt = list(dataType = 'Dense', methodSelectK = 'FVE'),  parallel = FALSE)
    for(p in 1:3){
      betaHat.sub[[p]] <- fit_fmm2d$betaHat[[p]]*mean(data$class == c) + betaHat.sub[[p]]
    }
    for(p in 1:3){
      betaHat.cov.sub[[p]] <- diag(fit_fmm2d$betaHat.cov[,,p])*mean(data$class == c)^2 +  betaHat.cov.sub[[p]]
    }
  }else{
    betaHat.sub <- betaHat.sub
    betaHat.cov.sub <- betaHat.cov.sub
  }
  
  
}


vars <- c("intercept", "gender", "stress")
pic <- list()
for (p in 1:length(betaHat.sub)){
  
  betaHat.up <- as.vector(betaHat.sub[[p]]) + 2*sqrt(betaHat.cov.sub[[p]])
  betaHat.low <- as.vector(betaHat.sub[[p]]) - 2*sqrt(betaHat.cov.sub[[p]])
  
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
ggarrange(pic[[1]],pic[[2]],pic[[3]], 
          ncol=3, nrow=1)





### smooth raw covariance estimates (estimates to measurement error)
if (smooth=="te"){
  sm <- data.frame(y = as.vector(RTilde), t=rep(1:T, S), day = rep(1:S, each=T)) 
  RHat <- gam(y ~ te(t, day, k = c(15,5)), data = sm)$fitted.values %>% array(dim=c(T,S))
} else if (smooth=="fbps"){
  RHat <- fbps(RTilde, knots = c(15,5))$Yhat 
}
RHat[which(RHat < 0)] <- 0


#draw pointwise covariance matrix
sub.avg <- data.frame(avg=as.vector(fit_fmm2d$betaHat.cov[1:1000,1:1000,2]), x=rep(1:1000, 1000), y=rep(1:1000, each=1000))

ggplot(sub.avg, aes(x, y, fill= avg)) + 
  geom_tile() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 13), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 13), axis.title.y = element_blank(), 
        legend.title = element_text(size = 14), legend.key.width = unit(1, 'cm'), 
        plot.title = element_text(size = 19)) + 
  scale_fill_gradient2(low="Steel Blue 3", mid="grey97", midpoint = 0,
                       high="red", name = "Cov") + 
  scale_y_reverse() 


library(RColorBrewer)
col<- colorRampPalette(c("blue", "white", "orange", "orangered", "red"))(10)
heatmap(cov.beta.ts.tilde[1:200,1:200,2], Rowv = NA, Colv = "Rowv", col = col)

edcomp <- eigen(cov.beta.ts.tilde[,,2])

