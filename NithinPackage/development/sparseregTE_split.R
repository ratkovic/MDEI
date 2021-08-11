library(splines2)
library(ranger)
library(microbenchmark)
library(sparsereg)
library(tictoc)
if(F){
  Rcpp::sourceCpp('~/Dropbox/Github/warmup/NithinPackage/src/splinesinter_short.cpp')
  Rcpp::sourceCpp('~/Dropbox/Github/warmup/NithinPackage/src/bayesLassoRevised.cpp')
  source('~/Dropbox/Github/warmup/NithinPackage/R/sparseregTE_auxiliary.R')
  n<-1000
  p<-5
  
  X <- matrix(rnorm(n*p),nrow=n)
  X<-apply(X,2,rank)
  X<- apply(X,2,scale)
  colnames(X) <- paste("X",1:ncol(X),sep="_")
  treat<- rnorm(n)
  theta.true <-treat^2# sign(X[,1]+X[,2]) * (treat)#treat^2*X[,2]^2#+treat
  y <-theta.true +rnorm(n)
  tau.true = 2*treat#sign(X[,1]+X[,2]) #0*treat#2*treat*X[,2]^2#treat*2*X[,2]+1# +X[,2]

  
  }



## First, turn covariates and treatment into spline bases, save these ----
tic()
sp1<-sparseregTE(y,treat,X,splits=10)
toc()

table(sign(apply(sp1$CIs.theta-theta.true,1,prod)))
table(sign(apply(sp1$CIs.tau-tau.true,1,prod)))


plot(treat,tau.true,type="l",ylim=range(sp1$CIs.tau))
segments(x0=treat,y0=sp1$CIs.tau[,1],y1=sp1$CIs.tau[,2],col=gray(.75))
points(treat,sp1$tau.est, pch=19,cex=0.5)


plot(treat,0*treat,type="l",ylim=range(sp1$CIs.tau))
segments(x0=treat,y0=sp1$CIs.tau[,1]-tau.true,y1=sp1$CIs.tau[,2]-tau.true,col=gray(.75))
points(treat,sp1$tau.est-tau.true, pch=19)

mean(sp1$tau.est)
mean(tau.true)

################

plot(treat[replaceme==1],ste.EM$fitted+mean(y))
treat.seq <- seq(min(treat),max(treat),length=1000)
lines(treat.seq,treat.seq^2)

cor(cbind(0,bases.obj$MConstructDerivative[replaceme==2,])%*%beta.sp,tau.true[replaceme==2])

tau.est<-cbind(0,bases.obj$MConstructDerivative[replaceme==2,])%*%beta.sp
plot(treat,treat*2,type="l",ylim=range(tau.est))
points(treat[replaceme==2],tau.est)
lines(treat.seq,2*treat.seq)

