
res1 <- read.csv("Results/result_decomposition_Survey_surrogate.csv")

res1 <- res1[!is.na(res1$CPE.env),]
## make log-transformation (not variables about shannon because shannon has been transformed)
res1[,c(3:9,11)] <- apply(res1[,c(3:9, 11)], 2, log10)
res1$portfolio <- res1$SAE + res1$CPE.env


res2 <- read.csv("Results/result_decomposition_Cedar_monoculture.csv")

res2 <- res2[!is.na(res2$CPE.env),]
## make log-transformation (not variables about shannon because shannon has been transformed)
res2[,c(3:9,11)] <- apply(res2[,c(3:9, 11)], 2, log10)
res2$portfolio <- res2$SAE + res2$CPE.env


source("Fn_myplot.R")


## decomposition ##
tiff("Figs/Fig_SI_SAE_prod_CPEenv.tif", width=10, height=8, 
     units="in", res=300, compression = "lzw")
myplot.avg(res1[,c("stab.com","portfolio","site")], fig.xlab=expression(SAE %*% CPE[env]), 
           fig.ylab=expression(Community~stability~","~S[com]), mycex=.9) +
  myplot.avg(res1[,c("portfolio","richness","site")], fig.xlab="Richness", 
             fig.ylab=expression(SAE%*%CPE[env]), mycex=.9) +
  myplot.each(res2[,c("stab.com","portfolio","site")], fig.xlab=expression(SAE%*%CPE[env]), 
             fig.ylab=expression(Community~stability~","~S[com]), mycex=.9) +
  myplot.each(res2[,c("portfolio","richness","site")], fig.xlab="Richness", 
             fig.ylab=expression(SAE%*%CPE[env]), mycex=.9) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect', ncol=2)&
  theme(legend.position="bottom")
dev.off()



