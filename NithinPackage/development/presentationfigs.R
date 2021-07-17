source('SparseregTE_split.R')
library(grf)
library(KRLS2)
library(sparsereg)

set.seed(1)
n<-1000
p <- 5
treat <- sort(rnorm(n))#seq(-1,1,length=n)
X<- matrix(rnorm(n*p), nrow=n)
colnames(X)<-paste("X",1:p,sep="_")

y.true <- treat^2#*treat#*sign(X[,1])
deriv.true <- (2*treat)#*sign(X[,1])
y <- y.true+rnorm(n,sd=1)

load('presrun.Rda')

#st1 <-sparseregTE(y,treat,X, splits=200)
g1 <- causal_forest(X,y,treat)

#save(st1,file="presrun.Rda")
pos <- !is.na(X[,1])

plot(treat[pos],st1$te.ave[pos],ylim=range(st1$CIs[pos,]),type="n")
segments(x0=treat[pos],x1=treat[pos],y0=st1$CIs[pos,1],y1=st1$CIs[pos,2],col=gray(.8))
points(treat[pos],st1$te.ave[pos],pch=19,cex=.2)
lines(treat[pos],deriv.true[pos])

k1<-krls(y=y,X=cbind(treat,X), epsilon = 0.001)
sk1<-inference.krls2(k1)
k1.est<-(sk1$derivatives[,1])
k1.se<-sk1$var.derivatives[,1]^.5
k1.ci<-k1.est+cbind(-1.645*k1.se, 1.645*k1.se)
plot(treat[pos],k1.est[pos],type="n",ylim=range(st1$CIs[pos,]))
segments(x0=treat[pos],x1=treat[pos],y0=k1.ci[pos,1],y1=k1.ci[pos,2],col=gray(.8))
points(treat[pos],k1.est[pos],pch=19,cex=.2)
lines(treat[pos],deriv.true[pos])

mean(apply(k1.ci-deriv.true,1,prod)<0) 
mean(apply(st1$CIs-deriv.true,1,prod)<0) 

cor(deriv.true[pos],k1.est[pos])
cor(deriv.true[pos],st1$te.ave[pos])

plot(treat,predict(g1)[,1],ylim=range(st1$CIs),type="n")
#segments(x0=treat[pos],x1=treat[pos],y0=predict(g1)[pos,2],y1=predict(g1)[pos,3],col=gray(.8))
points(treat,predict(g1)[,1],pch=19,cex=0.2)
lines(treat,deriv.true)

### Plots for pres

## The data
pdf("presfig1.pdf",h=6,w=6)
par("mar"=c(3,3,.2,.2))
plot(treat, treat^2, col=gray(.5),type="n",lwd=3,xlab="",ylab="",ylim=range(st1$CI))
points(treat,y,pch=19,cex=.25,col=gray(.8))
lines(treat,treat^2,col=gray(.2),lwd=3,lty=2)
lines(treat,2*treat,lwd=3)
mtext(outer=F,line=2,side=1,"T");mtext(outer=F,line=2,side=2,"Y")
legend("topleft",c("Conditional Mean","Slope"),lty=c(2,1),lwd=3,bty="n")
dev.off()


## The target
pdf("presfig2.pdf",h=6,w=6)
par("mar"=c(3,3,.2,.2))
plot(treat, treat^2, col=gray(.5),type="n",lwd=3,xlab="",ylab="",ylim=range(st1$CI))
#points(treat,y,pch=19,cex=.25,col=gray(.8))
#lines(treat,treat^2,col=gray(.2),lwd=3,lty=2)
lines(treat,2*treat,lwd=3)
mtext(outer=F,line=2,side=1,"T");mtext(outer=F,line=2,side=2,"Y")
dev.off()


## Fig 3: GRF Estimates
pdf("presfig3.pdf",h=6,w=6)
par("mar"=c(3,3,.2,.2))
plot(treat, treat^2, col=gray(.5),type="n",lwd=3,xlab="",ylab="",ylim=range(st1$CI))
#points(treat,y,pch=19,cex=.25,col=gray(.8))
#lines(treat,treat^2,col=gray(.2),lwd=3,lty=2)
lines(treat,2*treat,lwd=3,col=gray(.2))
points(treat,predict(g1)[,1],pch=19,cex=0.2)
mtext(outer=F,line=2,side=1,"T");mtext(outer=F,line=2,side=2,"Y")
dev.off()


## Fig 4: KRLS Estimates
pdf("presfig4.pdf",h=6,w=6)
par("mar"=c(3,3,.2,.2))
plot(treat, treat^2, col=gray(.5),type="n",lwd=3,xlab="",ylab="",ylim=range(st1$CI))
#points(treat,y,pch=19,cex=.25,col=gray(.8))
#lines(treat,treat^2,col=gray(.2),lwd=3,lty=2)
lines(treat,2*treat,lwd=3,col=gray(.2))
points(treat,k1.est,pch=19,cex=0.2)
mtext(outer=F,line=2,side=1,"T");mtext(outer=F,line=2,side=2,"Y")
dev.off()



## Fig 5: KRLS Intervals
pdf("presfig5.pdf",h=6,w=6)
par("mar"=c(3,3,.2,.2))
plot(treat, treat^2, col=gray(.5),type="n",lwd=3,xlab="",ylab="",ylim=range(st1$CI))
#points(treat,y,pch=19,cex=.25,col=gray(.8))
#lines(treat,treat^2,col=gray(.2),lwd=3,lty=2)
lines(treat,2*treat,lwd=3,col=gray(.2))
segments(x0=treat[pos],x1=treat[pos],y0=k1.ci[pos,1],y1=k1.ci[pos,2],col=gray(.8))
points(treat,k1.est,pch=19,cex=0.2)
mtext(outer=F,line=2,side=1,"T");mtext(outer=F,line=2,side=2,"Y")
#text(0,-10,"90% Intervals Cover the True Curve at 44% of Points")
dev.off()



## Fig 6: MDEI
pdf("presfig6.pdf",h=6,w=6)
par("mar"=c(3,3,.2,.2))
plot(treat, treat^2, col=gray(.5),type="n",lwd=3,xlab="",ylab="",ylim=range(st1$CI))
lines(treat,2*treat,lwd=3,col=gray(.2))
segments(x0=treat,x1=treat,y0=st1$CIs[,1],y1=st1$CIs[,2],col=gray(.8))
points(treat,st1$te.ave,pch=19,cex=0.2)
mtext(outer=F,line=2,side=1,"T");mtext(outer=F,line=2,side=2,"Y")
#text(0,-10,"90% Intervals Cover the True Curve at 44% of Points")
dev.off()
