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
 
makeHeatMap <- function (cluster1) {
  pairs.breaks <- seq(from=-2,to=2,length.out=257)
  red_green_bicluster<-heatmap.2(as.matrix(cluster1), col=greenred(256), dendrogram="both", 
                                 trace="none", breaks=pairs.breaks)
  red_green_bicluster
}

makeD3HeatMap<-function (cluster1) {
  rbg <- maPalette(low="darkblue", high="red4", mid="grey", k=200)
red_green_bicluster<-d3heatmap(as.matrix(cluster1), color=rbg, dendrogram="both", scale="row", k_row=3, k_col=3)
  red_green_bicluster 
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

 
 #================ The function that makes BayesBinom  ==================  
 # not working properly!!
 ##Beta distribution in biomaker discovery
 ##BayesBinom for each gene in the six geneLists(picker to choose different genes) 

 nthetapts=1000  #Resolution in theta for thetapts
 
 #Log Likelihood function:
 LL = function(theta,ns,Ntrials){
   ns*log(theta)+(Ntrials-ns)*log(1-theta)
 }
 posterior.poly = rgb(1,.5,0,.2)
 posterior.text = rgb(1,.5,0,1)
 
 
#  makeBayesBinom<- function(`pick.prior`=`pick.prior`, min=min, max=max, 
#                   `samp.size`=`samp.size`, `log.yaxis`=`log.yaxis`, dat=dat){
#    theta.pts = seq(0, 1, length = nthetapts)
#    uniform =function(theta){
#      if(min>max) {cat("Don't cross the streams!")}
#      dunif(theta, min = min, max=max)}
#    beta = function(theta){
#      min = min*10
#      max=max*10
#      dbeta(theta, min, max)
#    }
#    prior=list(uniform, beta)
#    prior.function = prior[[pick.prior]]
#    prior.pts = do.call( prior.function, list(theta.pts) )
#    #   prior.pts = prior[[pick.prior]](theta.pts)
#    
#    
#    baseCat <- sort(unique(dat))[1]
#    #Replicate data so you can sample more than exists! 
#    dat = rep(dat, times = 2)
#    nsuccesses=sum( dat[1:samp.size] == baseCat ) 
#    LL.pts = LL(theta.pts,nsuccesses,samp.size)
#    likeli.pts = exp(LL.pts)
#    posterior.pts = prior.pts*likeli.pts
#    posterior.pts = posterior.pts/mean(posterior.pts)
#    mypanel = function(x,y){
#      lpolygon(c(0,x,1),y=c(0,y,0), col=rgb(0,0,1,.1), border = FALSE)
#      lpolygon(c(0,x,1), c(0,posterior.pts,0), col = posterior.poly, border=FALSE)  
#      scaled.likeli = likeli.pts/max(likeli.pts)*max(prior.pts)
#      llines(x, scaled.likeli, col = "red", lwd = 2)
#      #Likelihood text
#      x.like.max = x[which(scaled.likeli==max(scaled.likeli))]
#      grid.text(x= unit(x.like.max,"native"), y = unit(max(scaled.likeli)*0.7,"native"), label = "Likelihood", rot = 270, gp = gpar(col = "red"))
#      ltext(x=x.like.max, y = max(scaled.likeli)*1.2, 
#            label= paste("theta =", signif(x.like.max, 3),
#                         "\nProb =", signif( max(likeli.pts), 3)),
#            col = "red")
#      #Prior Text:
#      succex = which(prior.pts != 0)
#      succex = x[succex]
#      if(pick.prior==1){
#        x.prior = succex[10]
#        ltext(x= x.prior, y = .9*prior[[pick.prior]](x.prior), col = "blue",
#              label="Prior", pos = 4)
#      }
#      if(pick.prior==2){
#        x.prior=x[which(prior.pts==max(prior.pts))]
#        ltext(x=x.prior, y = .9*prior[[pick.prior]](x.prior), col = "blue",
#              label="Prior")
#      }
#      #Posterior Text:
#      x.post.max = x[which(posterior.pts ==max(posterior.pts))]
#      y.post.max = posterior.pts[which(posterior.pts ==max(posterior.pts))]
#      ltext(x=x.post.max, y = 1*y.post.max, col = posterior.text, label = paste("Posterior:\n",
#                                                                                expression(theta),"=",signif(x.post.max,3),
#                                                                                "\nProb =", signif(y.post.max,3)))
#      
#    }
#    
#    graph<-xyplot(prior.pts~theta.pts, panel = mypanel, ylim = c(0,1.1*max(posterior.pts)), 
#                  xlab = expression(theta),ylab = "Probability Density",
#                  scales = list(x = list(log = FALSE), log=log.yaxis))
#    
#    graph
#    #The scales argument should accept a list, and x and y can be separate lists.
#    #On the internet I see many instances of y=list(log=TRUE) being the correct syntax. Not sure what's up.
#  }
#  
#  
 
