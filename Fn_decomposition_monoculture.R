##### X should be a population matrix (year * sp)
##### comm.list should be a list of pool population matrix, based on which a fake community would be assembled.
##### In each plot, there is a matrix of density arranged by year * species

##### Notice that each matrix in comm.list should have the same nrow with X

## evenness should be alpha.shannon / log(alpha.richness)
## Notice that Shannon index is calculated based on the sum across all years,
## and richness is the average number of sp. across years
## Evenness is also calculated based on the sum of each sp. across all years

decomposition.fake.comm <- function(X, comm.list){
  X1 <- X[,colSums(X, na.rm=T)>0]  #remove sp not showing
  
  ### make a big pool population matrix, containing all sp. from comm.list
  X.pool <- comm.list[[1]]
  for(i in 2:length(comm.list)){
    X.pool <- cbind(X.pool, comm.list[[i]])
  }
  sp.pool <- colnames(X.pool)
  
  
  
  #### get results

  #if all sp. in this year are NAs, delete this year from tot 
  i.not.na <- !(apply(X1, 1, function(xx) all(is.na(xx))))
  X2 <- X1[i.not.na,]
  
  tot <- rowSums(X2, na.rm=T)
  utot <- mean(tot)
  stab.com <- utot / sd(tot, na.rm=T)       #the reverse of CVcom
  p <- colSums(X2, na.rm=T) / sum(tot)
  alpha.inv.simpson <- 1/sum(p^2)      #Simpson
  alpha.shannon <- -sum(p * log(p))  #Shannon
  nn <- apply(X2, 1, function(xx) sum(xx>0))
  richness <- mean(nn)  #species richness
  evenness.shannon <- alpha.shannon / log(length(p))  #evenness of the distribution of density across sp.
  
  
  stab.sp <- utot / sum(apply(X2, 2, sd))   #the reverse of CVpop
  asyn <- stab.com / stab.sp      #based on LdM
  
  stab.com.ip <- utot / sqrt(sum(apply(X2, 2, var)))   #stable of com_ip, assuming all cov = 0
  SAE <- stab.com.ip / stab.sp
  
  
  
  #with 100 replicates
  CPE.env <- rep(NA,100)
  CPE.int <- rep(NA,100)
  for(rr in 1:100){
    #### make a fake community
    Y <- X1    #to store the fake community
    for(i in 1:ncol(X1)){
      tsp <- colnames(X1)[i]   #name of target species
      ind.tsp <- which(sp.pool == tsp)  #index of target sp in pool
      if(length(ind.tsp)==0){Error(paste0("No target sp ",tsp," in pool"))}
      i.sel <- sample(ind.tsp, 1)   #randomly choose a plot
      ts <- X.pool[,i.sel]
      Y[,i] <- (ts- mean(ts, na.rm=T)) / sd(ts, na.rm=T) * sd(X1[,i], na.rm=T) + mean(X1[,i], na.rm=T)
    }
    
    
    Y <- Y[i.not.na,]
    ## stab.com.fake, the community stability of fake community
    tot.fake <- rowSums(Y, na.rm=T)
    #stab.com.fake <- mean(tot.fake) / sd(tot.fake, na.rm=T) 
    stab.com.fake <- utot / sd(tot.fake, na.rm=T)       
    
    CPE.env[rr] <- stab.com.fake/ stab.com.ip    # CPE from response to environmental fluctuation
    CPE.int[rr] <- stab.com / stab.com.fake      # CPE from interaction
  }
  
  
  
    
  res <- data.frame(stab.com=stab.com,  stab.pop=stab.sp, asyn=asyn, SAE=SAE, 
                 CPE.env=mean(CPE.env, na.rm=T), CPE.int=mean(CPE.int, na.rm=T),
                 richness=richness, shannon=alpha.shannon,
                 inv.simpson=alpha.inv.simpson, evenness=evenness.shannon,
                 CPE.env.LB=quantile(CPE.env,0.025), CPE.env.UB=quantile(CPE.env,0.975),
                 CPE.int.LB=quantile(CPE.int, 0.025), CPE.int.UB=quantile(CPE.int, 0.975))

  return(res)
}


