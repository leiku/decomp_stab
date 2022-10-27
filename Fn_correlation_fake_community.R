##### X should be a population matrix (year * sp)
##### comm.list should be a list of pool population matrix, based on which a fake community would be assembled.
##### In each plot, there is a matrix of density arranged by year * species

##### Notice that each matrix in comm.list should have the same nrow with X

## Gross K, Cardinale BJ, and Fox JW, et al. 2014. Species Richness and the Temporal Stability of Biomass Production: 
## A New Analysis of Recent Biodiversity Experiments. Am Nat 183: 1-12.


correlation.fake.comm <- function(X, comm.list){
  X1 <- X[,colSums(X, na.rm=T)>0]  #remove sp not showing
  
  ### make a big pool population matrix, containing all sp. from comm.list
  X.pool <- comm.list[[1]]
  for(i in 2:length(comm.list)){
    X.pool <- cbind(X.pool, comm.list[[i]])
  }
  sp.pool <- colnames(X.pool)
  
  
  
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
  
  
  #### get results
  
  #if all sp. in this year are NAs, delete this year from tot 
  i.not.na <- !(apply(X1, 1, function(xx) all(is.na(xx))))
  X2 <- X1[i.not.na,]
  Y <- Y[i.not.na,]
  
  
  Gross.eta <- function(M){
    eta <- rep(NA, ncol(M))
    for(i in 1:ncol(M)){
      tot <- rowSums(M, na.rm=T)
      eta[i] <- cor(M[,i], tot-M[,i], use="pairwise.complete.obs")
    }
    return(mean(eta, na.rm=T))
  }
  
  
  return(data.frame(eta.real=Gross.eta(X2), eta.fake=Gross.eta(Y)))
}


