
fn_getstar <- function(p){
  if(p<0.001){
    star <- " ***"
  }else if(p<0.01){
    star <- " **"
  }else if(p<0.05){
    star <- " *"
  }else{star <- " "}
  return(star)
}


##### each point indicates the average of values across plots within one site 

myplot.avg <- function(D, fig.xlab="your.x", fig.ylab="your.y", mycex=1){
  colnames(D) <- c("y", "x", "Site")
  D$Site <- as.factor(D$Site)
  x <- tapply(D$x, D$Site, mean)
  x.sem <- tapply(D$x, D$Site, sd) / sqrt(table(D$Site))
  y <- tapply(D$y, D$Site, mean)
  y.sem <- tapply(D$y, D$Site, sd) / sqrt(table(D$Site))
  
  D.avg <- data.frame(x=x, y=y, x.sem=as.numeric(x.sem), y.sem=as.numeric(y.sem), Site=names(x))
  
  z <- summary(lm(y~x, D.avg))
  zz <- z$coef
  pvalue <- zz[2,4]
  if(pvalue>0.05){my.lty="dashed"}else{my.lty="solid"}
  mytext <- paste0("slope = ",round(zz[2,1], 3), fn_getstar(pvalue))
  mytext2 <-paste("R^{2}==", round(z$r.squared,2))
  
  #get position of text
  xrange <- range(D.avg$x-D.avg$x.sem)
  x1 <- xrange[1] + 0.0*diff(xrange)
  
  g1 <- ggplot(D.avg, aes(x=x, y=y)) + geom_point(size=1*mycex) + 
    geom_smooth(data=D.avg,method="lm", formula=y~x, se=T, col="black", linetype=my.lty)+
    geom_errorbar(aes(ymin=y-y.sem, ymax=y+y.sem))+
    geom_errorbarh(aes(xmin=x-x.sem, xmax=x+x.sem))+
    annotate("text", x = x1, y = Inf, hjust=0, vjust=1.3, label = mytext, parse=F, size=3*mycex, col="blue")+
    annotate("text", x = x1, y = Inf, hjust=0, vjust=2.6, label = mytext2, parse=T, size=3*mycex, col="blue")+
    geom_text_repel(label=D.avg$Site, size=2.5*mycex)+  #from "ggrepel", add labels for each point
    labs(x=fig.xlab, y=fig.ylab)+
    guides(color=guide_legend(nrow=1))+
    #xlim(min(x-x.sem), max(x+x.sem)) + 
    #ylim(min(y-y.sem), max(y+y.sem))+
    theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
          legend.position="bottom", text = element_text(size = 8*mycex))
  #stat_poly_eq(formula = y~x,parse=TRUE, aes(label=stat(eq.label),group=1),label.y = 0.95,label.x = 0.05)+
    #stat_poly_eq(formula = y~x,parse=TRUE,label.y = 0.85,label.x = 0.05,aes(label=stat(p.value.label),group=1))
  return(g1)
}




##### regress each site separately (for Cedar Greek)  ######
## only suitable for two sites

myplot.each <- function(D, fig.xlab="your.x", fig.ylab="your.y", mycex=1){
  colnames(D) <- c("y", "x", "Study")
  D$Study <- as.factor(D$Study)
  
  z1 <- summary(lm(y~x, subset(D, Study=="BigBio")))
  z2 <- summary(lm(y~x, subset(D, Study=="BioCON")))
  
  
  mytext1 <- paste0("slope = ",round(z1$coef[2,1], 3), fn_getstar(z1$coef[2,4]))
  mytext2 <- paste0("slope = ",round(z2$coef[2,1], 3), fn_getstar(z2$coef[2,4]))
  mytext1.R2 <- paste0("R^{2}==", round(z1$r.squared,2))
  mytext2.R2 <- paste0("R^{2}==", round(z2$r.squared,2))
  
  
  #mycolors<-c("#6495ED","#FFA500")
  mycolors<-c("blue","chocolate")
  
  #get position of text
  xrange <- range(D$x)
  yrange <- range(D$y)
  y.ub <- yrange[2] + 0.2*diff(yrange)
  x1 <- xrange[1] + 0.0*diff(xrange)
  x2 <- xrange[1] + 0.55*diff(xrange)
  y1 <- y.ub - 0.02 * diff(yrange)
  y2 <- y.ub - 0.15 * diff(yrange)
  
  g2 <- ggplot(data=D, aes(x=x, y=y,color=Study)) + 
    geom_point(alpha=0.6)+scale_color_manual(values=mycolors)+
    geom_smooth(method = "lm",formula = y ~ x,se=T,aes(group=Study))+
    annotate("text", x = x1, y = y1, hjust = 0, label = mytext1, parse=F, col=mycolors[1], size=3*mycex)+
    annotate("text", x = x1, y = y2, hjust = 0, label = mytext2, parse=F, col=mycolors[2], size=3*mycex)+
    annotate("text", x = x2, y = y1, hjust = 0, label = mytext1.R2, parse=T, col=mycolors[1], size=3*mycex)+
    annotate("text", x = x2, y = y2, hjust = 0, label = mytext2.R2, parse=T, col=mycolors[2], size=3*mycex)+
    labs(x=fig.xlab,y=fig.ylab)+
    ylim(min(D$y), y.ub)+
    theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
          legend.position="bottom", text = element_text(size = 8*mycex))+
    guides(color=guide_legend(nrow=1))
  
  return(g2)
}




##### regress a simple linear line  ######
## only suitable for one sites (i.e. no group)

myplot.simple <- function(D, fig.xlab="your.x", fig.ylab="your.y", mycex=1){
  colnames(D) <- c("y", "x")
  
  z1 <- summary(lm(y~x, D))
 
  mytext1 <- paste0("slope=",round(z1$coef[2,1], 3), fn_getstar(z1$coef[2,4]))
  mytext1.R2 <- paste0("R^{2}==", round(z1$r.squared,2))
  
  
  #mycolors<-c("#6495ED","#FFA500")
  mycolors<-c("blue")
  g1 <- ggplot(data=D, aes(x=x, y=y)) + 
    geom_point()+scale_color_manual(values=mycolors)+
    geom_smooth(method = "lm",formula = y ~ x,se=T)+
    annotate("text", x = -Inf, y = Inf, hjust=0, vjust=1.8, label = mytext1, parse=F, col=mycolors[1], size=5*mycex)+
    annotate("text", x = -Inf, y = Inf, hjust=-2, vjust=1.2, label = mytext1.R2, parse=T, col=mycolors[1], size=5*mycex)+
    labs(x=fig.xlab,y=fig.ylab)+
    ylim(min(D$y), max(D$y)+0.2)+
    theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
          legend.position="bottom", text = element_text(size = 15*mycex))+
    guides(color=guide_legend(nrow=1))
  
  return(g1)
}



