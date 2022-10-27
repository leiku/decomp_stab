

source("Fn_decomposition_monoculture.R")


####### BigBio
ds <- readRDS("Data/Data_Cedar_clean.Rds")

A <- ds[[1]]
n.sp <- unlist(lapply(A, ncol))
table(n.sp)

A.mono <- A[n.sp==1]
A.mix <- A[n.sp>1]


### get results
set.seed(10)

res.bigbio <- data.frame(site=NA, plot=NA, stab.com=NA, 
                         stab.pop=NA, asyn=NA, 
                         SAE=NA, CPE.env=NA, CPE.int=NA, 
                         richness=NA, shannon=NA, 
                         inv.simpson=NA,evenness=NA,
                         CPE.env.LB=NA, CPE.env.UB=NA,
                         CPE.int.LB=NA, CPE.int.UB=NA)
for(i in 1:length(A.mix)){
  A1 <- A.mix[[i]]
  res1 <- decomposition.fake.comm(A1, A.mono)
  res1$site <- "BigBio"
  res1$plot <- names(A.mix)[i]
  res1 <- res1[, c(15,16,1:14)]
  
  res.bigbio <- rbind(res.bigbio, res1)
  
}
res.bigbio <- res.bigbio[-1,]


###### BioCON  #####

B <- ds[[2]]
n.sp <- unlist(lapply(B, ncol))
table(n.sp)

B.mono <- B[n.sp==1]
B.mix <- B[n.sp>1]


### get results
set.seed(10)

res.biocon <- data.frame(site=NA, plot=NA, stab.com=NA, 
                         stab.pop=NA, asyn=NA, 
                         SAE=NA, CPE.env=NA, CPE.int=NA, 
                         richness=NA, shannon=NA, 
                         inv.simpson=NA,evenness=NA,
                         CPE.env.LB=NA, CPE.env.UB=NA,
                         CPE.int.LB=NA, CPE.int.UB=NA)
for(i in 1:length(B.mix)){
  B1 <- B.mix[[i]]
  res1 <- decomposition.fake.comm(B1, B.mono)
  res1$site <- "BioCON"
  res1$plot <- names(B.mix)[i]
  res1 <- res1[, c(15,16,1:14)]
  
  res.biocon <- rbind(res.biocon, res1)
  
}
res.biocon <- res.biocon[-1,]

res <- rbind(res.bigbio, res.biocon)

write.csv(res, "Results/result_decomposition_Cedar_monoculture.csv", row.names = F)



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
res <- read.csv("Results/result_decomposition_Cedar_monoculture.csv")

## make log-transformation (not variables about shannon because shannon has been transformed)
res[,c(3:9,11, 13:16)] <- apply(res[,c(3:9, 11, 13:16)], 2, log10)
res$CPE <- res$CPE.env + res$CPE.int

res$portfolio <- res$SAE + res$CPE.env

source("Fn_myplot.R")


tiff("Figs/Fig_SI_Cedar_4parts.tif", width=10, height=5, 
     units="in", res=300, compression = "lzw")

myplot.each(res[,c("CPE.env","richness","site")], fig.xlab="Richness", 
              fig.ylab=expression(CPE[env]), mycex=.9) +
  myplot.each(res[,c("CPE.int","richness","site")], fig.xlab="Richness", 
              fig.ylab=expression(CPE[int]), mycex=.9) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect', ncol=2)&
  theme(legend.position="bottom")
dev.off()

