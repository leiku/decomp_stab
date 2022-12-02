
source("Fn_decomposition.R")



##### get data #####


A <- readRDS("Data/Data_Survey.RDS")

res <- data.frame(site=rep(NA,1000), plot=rep(NA,1000), stab.comm=rep(NA,1000), 
                  stab.pop=rep(NA,1000), asynchrony=rep(NA,1000), compensatory=rep(NA,1000), statistical=rep(NA,1000),
                  richness=rep(NA, 1000), shannon=rep(NA,1000), inv.simpson=rep(NA,1000),evenness=rep(NA,1000))
k <- 1
for(i in 1:9){
  A1 <- A[[i]]
  for(j in 1:length(A1)){
    res$site[k] <- names(A)[i]
    res$plot[k] <- names(A1)[j]
    
    tmp <- decomposition(A1[[j]])
    res[k, 3:11] <- tmp
    
    k <- k+1
  }
}
res <- res[!is.na(res$site),]


write.csv(res, "Results/result_decomposition_Survey_3parts.csv", row.names = F)

####### make figures ########
res <- read.csv("Results/result_decomposition_Survey_3parts.csv")
## make log-transformation (not variables about shannon because shannon has been transformed)
res[,c(3:8,10)] <- apply(res[,c(3:8, 10)], 2, log10)

source("Fn_myplot.R")


## decomposition ##

pdf("Figs/Fig_Survey_decomposition.pdf", width=7, height=6)
myplot.avg(res[,c("stab.comm","stab.pop","site")], fig.xlab=expression(Population~stability~","~S[pop]), 
           fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  myplot.avg(res[,c("stab.comm","asynchrony","site")], fig.xlab=expression(Asynchrony~","~Phi), 
             fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  myplot.avg(res[,c("stab.comm","compensatory","site")], fig.xlab=expression(Compensatory~effect~","~CPE), 
             fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  myplot.avg(res[,c("stab.comm","statistical","site")], fig.xlab=expression(Statistical-averaging~effect~","~SAE), 
             fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect')
dev.off()

designs <- c(area(1,2,1,3), area(2,1,2,2), area(2,3,2,4), area(3,1,3,2), area(3,3,3,4))


## average richness across years ##
pdf("Figs/Fig_Survey_richness.pdf", width=7, height=9)
myplot.avg(res[,c("stab.comm","richness","site")], fig.xlab="Richness", fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  myplot.avg(res[,c("stab.pop","richness","site")], fig.xlab="Richness", fig.ylab=expression(Population~stability~","~S[pop]), mycex=1.3) +
  myplot.avg(res[,c("asynchrony","richness","site")], fig.xlab="Richness", fig.ylab=expression(Asynchrony~","~Phi), mycex=1.3) +
  myplot.avg(res[,c("compensatory","richness","site")], fig.xlab="Richness", fig.ylab=expression(Compensatory~effect~","~CPE), mycex=1.3) +
  myplot.avg(res[,c("statistical","richness","site")], fig.xlab="Richness", fig.ylab=expression(Statistical-averaging~effect~","~SAE), mycex=1.3) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(design=designs, widths=c(3,3,3,4.5), guides = 'collect')
dev.off()



## shannon ##
tiff("Figs/Fig_Survey_shannon.tif", width=7, height=9, 
     units="in", res=300, compression = "lzw")
myplot.avg(res[,c("stab.comm","shannon","site")], fig.xlab="Shannon index", fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  myplot.avg(res[,c("stab.pop","shannon","site")], fig.xlab="Shannon index", fig.ylab=expression(Population~stability~","~S[pop]), mycex=1.3) +
  myplot.avg(res[,c("asynchrony","shannon","site")], fig.xlab="Shannon index", fig.ylab=expression(Asynchrony~","~Phi), mycex=1.3) +
  myplot.avg(res[,c("compensatory","shannon","site")], fig.xlab="Shannon index", fig.ylab=expression(Compensatory~effect~","~CPE), mycex=1.3) +
  myplot.avg(res[,c("statistical","shannon","site")], fig.xlab="Shannon index", fig.ylab=expression(Statistical-averaging~effect~","~SAE), mycex=1.3) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(design=designs, widths=c(3,3,3,4.5), guides = 'collect')
dev.off()



## inv.simpson ##
tiff("Figs/Fig_Survey_inv.simpson.tif", width=7, height=9, 
     units="in", res=300, compression = "lzw")
myplot.avg(res[,c("stab.comm","inv.simpson","site")], fig.xlab="Inverse Simpson", fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  myplot.avg(res[,c("stab.pop","inv.simpson","site")], fig.xlab="Inverse Simpson", fig.ylab=expression(Population~stability~","~S[pop]), mycex=1.3) +
  myplot.avg(res[,c("asynchrony","inv.simpson","site")], fig.xlab="Inverse Simpson", fig.ylab=expression(Asynchrony~","~Phi), mycex=1.3) +
  myplot.avg(res[,c("compensatory","inv.simpson","site")], fig.xlab="Inverse Simpson", fig.ylab=expression(Compensatory~effect~","~CPE), mycex=1.3) +
  myplot.avg(res[,c("statistical","inv.simpson","site")], fig.xlab="Inverse Simpson", fig.ylab=expression(Statistical-averaging~effect~","~SAE), mycex=1.3) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(design=designs, widths=c(3,3,3,4.5), guides = 'collect')
dev.off()



