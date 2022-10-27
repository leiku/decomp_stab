
library(wsyn)
source("Fn_decomposition.R")


######## Survey ###########

A <- readRDS("Data/Data_Survey.RDS")


Scomip.surr <- vector("list", length(A))  # to store all stab.comip.surr
names(Scomip.surr) <- names(A) 

n.surr <- 1000 #number of surrogates
res.surr <- data.frame(site=rep(NA,1000), plot=rep(NA,1000), stab.comm=rep(NA,1000), 
                       stab.pop=rep(NA,1000), CPE=rep(NA,1000), SAE=rep(NA,1000), 
                       CPE.lb=rep(NA,1000), CPE.ub=rep(NA,1000),
                       SAE.lb=rep(NA,1000), SAE.ub=rep(NA,1000)) #store 95% CI
k <- 1
for(i in 1:9){
  A1 <- A[[i]]
  scomip.surr <- matrix(NA, n.surr, length(A1)) # a matrix to store all SAE
  colnames(scomip.surr) <- names(A1)
  for(j in 1:length(A1)){
    res.surr$site[k] <- names(A)[i]
    res.surr$plot[k] <- names(A1)[j]
    
    tmp <- decomposition(A1[[j]])
    res.surr[k, 3:6] <- tmp[c(1:2,4:5)]
    
    #### surrogate
    X <- A1[[j]]
    X <- X[, colSums(X)>0]
    mu.tot <- mean(rowSums(X))
    
    Y <- cleandat(t(X), times=1:nrow(X), clev=1)   #demean (which should not affect on sd)
    mat.surr <- surrog(Y$cdat, nsurrogs = n.surr, surrtype = "aaft", syncpres = F)
    sd.tot.surr <- sapply(mat.surr, function(xx) sd(colSums(xx)))
    stab.comip.surr <- mu.tot / sd.tot.surr
    scomip.surr[,j] <- stab.comip.surr
    sae.surr <- stab.comip.surr / tmp[2]   #SAE = S.comip / Spop
    cpe.surr <- tmp[1] / stab.comip.surr   #CPE = Scom / S.comip
    
    cpe.CI <-  quantile(cpe.surr, c(0.025,0.975))
    sae.CI <-  quantile(sae.surr, c(0.025,0.975))
    
    res.surr[k, 7:10] <- c(cpe.CI, sae.CI)
    
    k <- k+1
  }
  Scomip.surr[[i]] <- scomip.surr
  print(i)
}
res.surr <- res.surr[!is.na(res.surr$site),]


## make a table to show significance of CPE
sig <- rep("NS", nrow(res.surr))
sig[res.surr$CPE.lb>1] <- "compensatory" 
sig[res.surr$CPE.ub<1] <- "synchronous" 
res.surr$sig <- sig
print(round(prop.table(table(res.surr$sig, res.surr$site), margin=2),3))






####### Cedar: AAFT surrogates to get CIs ###########
ds <- readRDS("Data/Data_Cedar_clean.Rds")

for(i in 1:2){
  A <- ds[[i]]
  n.sp <- unlist(lapply(A, ncol))
  table(n.sp)
  
  A.mix <- A[n.sp>1]
  ds[[i]] <- A.mix
}



Scomip.surr <- vector("list", length(ds))  # to store all stab.comip.surr
names(Scomip.surr) <- names(ds) 

n.surr <- 1000 #number of surrogates
res.surr <- data.frame(site=rep(NA,1000), plot=rep(NA,1000), stab.comm=rep(NA,1000), 
                       stab.pop=rep(NA,1000), CPE=rep(NA,1000), SAE=rep(NA,1000), 
                       CPE.lb=rep(NA,1000), CPE.ub=rep(NA,1000),
                       SAE.lb=rep(NA,1000), SAE.ub=rep(NA,1000)) #store 95% CI
k <- 1
for(i in 1:2){
  ds1 <- ds[[i]]
  scomip.surr <- matrix(NA, n.surr, length(ds1)) # a matrix to store all SAE
  colnames(scomip.surr) <- names(ds1)
  for(j in 1:length(ds1)){
    res.surr$site[k] <- names(ds)[i]
    res.surr$plot[k] <- names(ds1)[j]
    
    
    tmp <- decomposition(ds1[[j]])
    res.surr[k, 3:6] <- tmp[c(1:2,4:5)]
    
    #### surrogate
    X <- ds1[[j]]
    X <- X[, colSums(X, na.rm=T)>0]
    X <- X[rowSums(X,na.rm=T)>0,]
    mu.tot <- mean(rowSums(X))
    
    Y <- cleandat(t(X), times=1:nrow(X), clev=1)   #demean (which should not affect on sd)
    mat.surr <- surrog(Y$cdat, nsurrogs = n.surr, surrtype = "aaft", syncpres = F)
    sd.tot.surr <- sapply(mat.surr, function(xx) sd(colSums(xx)))
    stab.comip.surr <- mu.tot / sd.tot.surr
    scomip.surr[,j] <- stab.comip.surr
    sae.surr <- stab.comip.surr / tmp[2]   #SAE = S.comip / Spop
    cpe.surr <- tmp[1] / stab.comip.surr   #CPE = Scom / S.comip
    
    cpe.CI <-  quantile(cpe.surr, c(0.025,0.975))
    sae.CI <-  quantile(sae.surr, c(0.025,0.975))
    
    res.surr[k, 7:10] <- c(cpe.CI, sae.CI)
    
    k <- k+1
  }
  Scomip.surr[[i]] <- scomip.surr
  print(i)
}
res.surr <- res.surr[!is.na(res.surr$site),]


## make a table to show significance of CPE
sig <- rep("NS", nrow(res.surr))
sig[res.surr$CPE.lb>1] <- "compensatory" 
sig[res.surr$CPE.ub<1] <- "synchronous" 
res.surr$sig <- sig
print(t(round(prop.table(table(res.surr$sig, res.surr$site), margin=2),3)))



