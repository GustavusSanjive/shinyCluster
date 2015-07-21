library("matrixStats", lib.loc="D:/R-3.1.0/library")
library("mosaic", lib.loc="D:/R-3.1.0/library")
library("dendextend", lib.loc="D:/R-3.1.0/library")
library("scales", lib.loc="D:/R-3.1.0/library")
library("colorspace", lib.loc="D:/R-3.1.0/library")
library("gplots", lib.loc="D:/R-3.1.0/library")
library("dplyr", lib.loc="D:/R-3.1.0/library")
library("reshape2", lib.loc="D:/R-3.1.0/library")
library("grid", lib.loc="D:/R-3.1.0/library")
library("gridExtra", lib.loc="D:/R-3.1.0/library")
library("RColorBrewer", lib.loc="D:/R-3.1.0/library")
library("marray", lib.loc="D:/R-3.1.0/library")
library("d3heatmap", lib.loc="D:/R-3.1.0/library")


#================ The function that generates heatmap ===================== 
# The function takes a gene dataset as an input and generates a heatmap
 

makeD3HeatMap<-function (cluster1) {
  rbg <- maPalette(low="darkblue", high="red4", mid="grey", k=200)
  bicluster<-d3heatmap(as.matrix(cluster1), color=rbg, dendrogram="both", scale="row", k_row=3, k_col=3)
  bicluster 
}


meltCluster<-function(cluster1){
  cluster1<-mutate(cluster1, probes=rownames(cluster1))
  melt(cluster1)}

ANOVA_probes<-function(cluster1){
  mod<-lm(value ~ x*Norm_HR, data=cluster1)
  mplot(mod, which=c(1,2,7), multiplot = TRUE, ncol=3)}

ANOVA_TukeyHSD<-function(cluster1){
  mod<-lm(value ~ x*Norm_HR, data=cluster1)
  mplot(TukeyHSD(mod), system="ggplot2")}

ANOVA_Summary<-function(cluster1){
  mod<-lm(value ~ x*Norm_HR, data=cluster1)
  summary(mod)}

ANOVA_pval<-function(cluster1){
  mod<-lm(value ~ x*Norm_HR, data=cluster1)
  anova(mod)}

 
 
 #================ The function that makes power test   ================== 
 #' The user may choose to calculate power of a hypothesis test. There are sliders in the manipulate
 #' box to interact with mu, sigma and sample size.
 
 makePower<- function(mu=mu, sigma=sigma, n=n)
 { #set.seed(seed)
   samp <- rnorm(n,mean=mu,sd=sigma)
   mu0 <-5
   x.bar <- mean(samp)
   margin <- qt(0.975,df=n-1)*sd(samp)/sqrt(n)
   low <- x.bar-margin
   high <- x.bar+margin
   t.stat <- (x.bar-mu0)*sqrt(n)/sd(samp)
   p.value <- round(2*(1-pt(abs(t.stat),df=n-1)),3)
   hist=histogram(samp,ylab="Density", type="density",
                  xlim=c(0,10),xaxt="n",xlab="",
                  main="Power of a Hypothesis (t) Test",
                  sub=paste("P-value = ",p.value),
                  col=rgb(0,.5,.9,.7),
                  panel=function(x,...){
                    panel.histogram(x,...)
                    #panel.densityplot(x, col="orange")
                    llines(x=c(low,high), y=c(0,0), col="red", lwd=10)
                    lpoints(x=x.bar, y=0, col="green", cex=1.3, pch=19)
                    grid.text(label=expression(mu), x=unit(x.bar,"native"), y=unit(.1, "npc"))
                    lpoints(x=mu0, y=0, col="black", cex=1.3, pch=19)
                    grid.text(label=expression(mu[0]==5), x=unit(mu0,"native"), y=unit(.03, "npc"))
                  })
   print(hist)
 }

 
 