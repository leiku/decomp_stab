
source("Fn_decomposition.R")

ds <- readRDS("Data/Data_Cedar_clean.Rds")

for(i in 1:2){
  A <- ds[[i]]
  n.sp <- unlist(lapply(A, ncol))
  table(n.sp)
  
  A.mix <- A[n.sp>1]
  ds[[i]] <- A.mix
}



res <- data.frame(site=rep(NA,1000), plot=rep(NA,1000), stab.comm=rep(NA,1000), 
                  stab.pop=rep(NA,1000), asynchrony=rep(NA,1000), compensatory=rep(NA,1000), statistical=rep(NA,1000),
                  richness=rep(NA, 1000), shannon=rep(NA,1000), inv.simpson=rep(NA,1000),evenness=rep(NA,1000))
k <- 1
for(i in 1:2){
  ds1 <- ds[[i]]
  for(j in 1:length(ds1)){
    res$site[k] <- names(ds)[i]
    res$plot[k] <- names(ds1)[j]
    
    tmp <- decomposition(ds1[[j]])
    res[k, 3:11] <- tmp
    
    k <- k+1
  }
}
res <- res[!is.na(res$site),]

write.csv(res, "Results/result_decomposition_Cedar.csv", row.names = F)


######### plot  ########
res <- read.csv("Results/result_decomposition_Cedar.csv")

## make log-transformation (not variables about shannon because shannon has been transformed)
res[,c(3:8, 10)] <- apply(res[,c(3:8, 10)], 2, log10)


source("Fn_myplot.R")


##decomposition ##
pdf("Figs/Fig_Cedar_decomposition.pdf", width=7, height=6)
myplot.each(res[,c("stab.comm","stab.pop","site")], fig.xlab=expression(Population~stability~","~S[pop]), 
            fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  myplot.each(res[,c("stab.comm","asynchrony","site")], fig.xlab=expression(Asynchrony~","~Phi), 
              fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  myplot.each(res[,c("stab.comm","compensatory","site")], fig.xlab=expression(Compensatory~effect~","~CPE), 
              fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  myplot.each(res[,c("stab.comm","statistical","site")], fig.xlab=expression(Statistical-averaging~effect~","~SAE), 
              fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect')&
  theme(legend.position="bottom")
dev.off()


designs <- c(area(1,2,1,3), area(2,1,2,2), area(2,3,2,4), area(3,1,3,2), area(3,3,3,4))

## richness ##
pdf("Figs/Fig_Cedar_richness.pdf", width=7, height=9)
myplot.each(res[,c("stab.comm","richness","site")], fig.xlab="Richness", 
            fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  myplot.each(res[,c("stab.pop","richness","site")], fig.xlab="Richness", 
              fig.ylab=expression(Population~stability~","~S[pop]), mycex=1.3) +
  myplot.each(res[,c("asynchrony","richness","site")], fig.xlab="Richness", 
              fig.ylab=expression(Asynchrony~","~Phi), mycex=1.3) +
  myplot.each(res[,c("compensatory","richness","site")], fig.xlab="Richness", 
              fig.ylab=expression(Compensatory~effect~","~CPE), mycex=1.3) +
  myplot.each(res[,c("statistical","richness","site")], fig.xlab="Richness", 
              fig.ylab=expression(Statistical-averaging~effect~","~SAE), mycex=1.3) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(design=designs, widths=c(3,3,3,4.5), guides = 'collect')&
  theme(legend.position="bottom")
dev.off()


## shannon ##

tiff("Figs/Fig_Cedar_shannon.tif", width=7, height=9, 
     units="in", res=300, compression = "lzw")
myplot.each(res[,c("stab.comm","shannon","site")], fig.xlab="Shannon index", 
            fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  myplot.each(res[,c("stab.pop","shannon","site")], fig.xlab="Shannon index", 
              fig.ylab=expression(Population~stability~","~S[pop]), mycex=1.3) +
  myplot.each(res[,c("asynchrony","shannon","site")], fig.xlab="Shannon index", 
              fig.ylab=expression(Asynchrony~","~Phi), mycex=1.3) +
  myplot.each(res[,c("compensatory","shannon","site")], fig.xlab="Shannon index",
              fig.ylab=expression(Compensatory~effect~","~CPE), mycex=1.3) +
  myplot.each(res[,c("statistical","shannon","site")], fig.xlab="Shannon index", 
              fig.ylab=expression(Statistical-averaging~effect~","~SAE), mycex=1.3) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(design=designs, widths=c(3,3,3,4.5), guides = 'collect')&
  theme(legend.position="bottom")
dev.off()



## inv.simpson ##
designs <- c(area(1,2,1,3), area(2,1,2,2), area(2,3,2,4), area(3,1,3,2), area(3,3,3,4))

tiff("Figs/Fig_Cedar_inv.simpson.tif", width=7, height=9, 
     units="in", res=300, compression = "lzw")
myplot.each(res[,c("stab.comm","inv.simpson","site")], fig.xlab="Inverse Simpson", 
            fig.ylab=expression(Community~stability~","~S[com]), mycex=1.3) +
  myplot.each(res[,c("stab.pop","inv.simpson","site")], fig.xlab="Inverse Simpson", 
              fig.ylab=expression(Population~stability~","~S[pop]), mycex=1.3) +
  myplot.each(res[,c("asynchrony","inv.simpson","site")], fig.xlab="Inverse Simpson", 
              fig.ylab=expression(Asynchrony~","~Phi), mycex=1.3) +
  myplot.each(res[,c("compensatory","inv.simpson","site")], fig.xlab="Inverse Simpson", 
              fig.ylab=expression(Compensatory~effect~","~CPE), mycex=1.3) +
  myplot.each(res[,c("statistical","inv.simpson","site")], fig.xlab="Inverse Simpson", 
              fig.ylab=expression(Statistical-averaging~effect~","~SAE), mycex=1.3) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(design=designs, widths=c(3,3,3,4.5), guides = 'collect')&
  theme(legend.position="bottom")
dev.off()



