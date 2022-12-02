

source("Fn_decomposition_surrogate_covariance.R")



##### get data #####
set.seed(1)

A <- readRDS("Data/Data_Survey.RDS")


res <- data.frame(site=NA, plot=NA, stab.com=NA, 
                  stab.pop=NA, asyn=NA, 
                  SAE=NA, CPE.env=NA, CPE.int=NA, 
                  richness=NA, shannon=NA, 
                  inv.simpson=NA,evenness=NA,
                  CPE.env.LB=NA, CPE.env.UB=NA,
                  CPE.int.LB=NA, CPE.int.UB=NA)
for(i in 1:9){
  A1 <- A[[i]]
  res1 <- decomposition.surrogate.covariance(A1)
  res1$site <- names(A)[i]
  res1$plot <- names(A1)
  res1 <- res1[, c(15,16,1:14)]
  
  res <- rbind(res, res1)
  print(i)
}
res <- res[-1,]


write.csv(res, "Results/result_decomposition_Survey_4parts.csv", row.names = F)


## uncertainty test
sig <- rep("NS", nrow(res))
sig[res$CPE.env.LB>1] <- "compensatory" 
sig[res$CPE.env.UB<1] <- "synchronous" 
res$sig <- sig
print("Significance of CPE.env:")
print(t(round(prop.table(table(res$sig, res$site), margin=2),3)))


sig <- rep("NS", nrow(res))
sig[res$CPE.int.LB>1] <- "compensatory" 
sig[res$CPE.int.UB<1] <- "synchronous" 
res$sig <- sig
print("Significance of CPE.int:")
print(t(round(prop.table(table(res$sig, res$site), margin=2),3)))


####### make figures ########
res <- read.csv("Results/result_decomposition_Survey_4parts.csv")

res <- res[!is.na(res$CPE.env),]
## make log-transformation (not variables about shannon because shannon has been transformed)
res[,c(3:9,11)] <- apply(res[,c(3:9, 11)], 2, log10)
res$CPE <- res$CPE.env + res$CPE.int
res$portfolio <- res$SAE + res$CPE.env

source("Fn_myplot.R")


## decomposition ##
tiff("Figs/Fig_SI_Survey_4parts.tif", width=7, height=6, 
     units="in", res=300, compression = "lzw")
myplot.avg(res[,c("CPE.env","richness","site")], fig.xlab="Richness", 
              fig.ylab=expression(CPE[env]), mycex=1.3) +
  myplot.avg(res[,c("CPE.int","richness","site")], fig.xlab="Richness", 
              fig.ylab=expression(CPE[int]), mycex=1.3) +
  myplot.avg(res[,c("stab.com","portfolio","site")], fig.xlab=expression(SAE %*% CPE[env]), 
             fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  myplot.avg(res[,c("portfolio","richness","site")], fig.xlab="Richness", 
             fig.ylab=expression(SAE%*%CPE[env]), mycex=1.3) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect', ncol=2)&
  theme(legend.position="bottom")
dev.off()



