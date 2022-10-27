##### X should be a list containing several plots.
##### In each plot, there is a matrix of density arranged by year * species
## evenness should be alpha.shannon / log(alpha.richness)
## Notice that shannon index is calculated based on the sum across all years,
## and richness is the number of sp. exist in at least one year
## Evenness is also calculated based on the sum of each sp. across all years

decomposition.surrogate.covariance <- function(X){
  X1 <- lapply(X, function(Y) Y[,colSums(Y, na.rm=T)>0])  #remove sp not showing
  n.sp <- unlist(lapply(X1, ncol))  #check number of species for each plot
  
  #make a matrix to show if a species shows in a plot
  all.sp <- sort(unique(unlist(lapply(X1, colnames))))
  sp.plot <- matrix(0, length(X1), length(all.sp))
  rownames(sp.plot) <- names(X1)
  colnames(sp.plot) <- all.sp
  for(i in 1:length(X1)){
    sp.plot[i, all.sp %in% colnames(X1[[i]])] <- 1
  }
  
  
  
  res <- matrix(NA, length(X1), 14)
  colnames(res) <- c("stab.com", "stab.pop", 
                     "asyn", "SAE", 
                     "CPE.env", "CPE.int", 
                     "richness", "shannon",
                     "inv.simpson", "evenness", 
                     "CPE.env.LB","CPE.env.UB",
                     "CPE.int.LB", "CPE.int.UB")
  for(i in 1:length(X1)){
    X2 <- X1[[i]]
    X2 <- X2[,order(colnames(X2))] #reorder columns
    
    #if all sp. in this year are NAs, delete this year from tot 
    i.not.na <- !(apply(X2, 1, function(xx) all(is.na(xx))))
    X3 <- X2[i.not.na,]
    
    tot <- rowSums(X3, na.rm=T)
    utot <- mean(tot)
    stab.com <- utot / sd(tot)       #the reverse of CVcom
    p <- colSums(X3, na.rm=T) / sum(tot)
    alpha.inv.simpson <- 1/sum(p^2)      #Simpson
    alpha.shannon <- -sum(p * log(p))  #Shannon
    nn <- apply(X3, 1, function(xx) sum(xx>0))
    richness <- mean(nn)  #species richness
    evenness.shannon <- alpha.shannon / log(length(p))  #evenness of the distribution of density across sp.
    
    stab.sp <- utot / sum(apply(X3, 2, sd))   #the reverse of CVpop
    asyn <- stab.com / stab.sp      #based on LdM
    
    stab.com.ip <- utot / sqrt(sum(apply(X3, 2, var)))
    SAE <- stab.com.ip / stab.sp
    
    #with 100 replicates
    CPE.env <- rep(NA,100)
    CPE.int <- rep(NA,100)
    for(rr in 1:100){
      mat.cov <- cov(X3)
      for(i1 in 1:ncol(mat.cov)){
        abun1 <- X3[, i1]  #timeseries of the first species
        for(j1 in 1:ncol(mat.cov)){
          if(j1 != i1){
            name.sp <- colnames(X2)[j1]  #name of sp to be chosen
            #show which plots this sp exist in
            tmp <- sp.plot[, colnames(sp.plot) == name.sp]
            tmp <- tmp[-i]
            #randomly choose a plot if there is more than 0
            ind <- 1:length(X1)
            ind <- ind[-i]
            if(sum(tmp)==0){
              mat.cov[i1,j1] <- 0  #if this sp only exist in the current plot
            }else{
              if(sum(tmp)==1){
                i.chose <- ind[which(tmp==1)]
              }else(i.chose <- ind[sample(which(tmp==1),1)])
              X.chose <- X1[[i.chose]]  # the population matrix chosen
              #get the timeseries of the chosen sp. with the same years of the target commnity
              ts2 <- X3[, j1];
              abun2 <- X.chose[i.not.na, colnames(X.chose)==name.sp] 
              abun2 <- (abun2 - mean(abun2, na.rm=T)) / sd(abun2, na.rm=T) * sd(ts2) + mean(ts2)
              mat.cov[i1,j1] <- cov(abun1, abun2, use="pairwise.complete.obs") 
            }
          }
        }
      }
      stab.com.sur <- utot / sqrt(sum(mat.cov, na.rm=T)) 
      CPE.env[rr] <- stab.com.sur/ stab.com.ip    # CPE from response to environmental fluctuation
      CPE.int[rr] <- stab.com / stab.com.sur      # CPE from interaction
    }
    
    
    
    res[i,] <- c(stab.com=stab.com,  stab.pop=stab.sp, asyn=asyn, SAE=SAE, 
      CPE.env=mean(CPE.env, na.rm=T), CPE.int=mean(CPE.int, na.rm=T),
      richness=richness, shannon=alpha.shannon,
      inv.simpson=alpha.inv.simpson, evenness=evenness.shannon,
      CPE.env.LB=quantile(CPE.env,0.025, na.rm=T), CPE.env.UB=quantile(CPE.env,0.975, na.rm=T),
      CPE.int.LB=quantile(CPE.int, 0.025, na.rm=T), CPE.int.UB=quantile(CPE.int, 0.975, na.rm=T))
  }
  
  return(as.data.frame(res))
}


