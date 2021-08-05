library(splines2)
library(ranger)
library(microbenchmark)
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
  
  y <- treat^2 + X[,3]+rnorm(n)
  tau.true = 2*treat# +X[,2]

  
  }



## First, turn covariates and treatment into spline bases, save these ----
n<-length(treat)
treatmat <- bs.me(treat)
Xmat<-cbind(1,matrix(apply(X,2,bs.me),nrow=n))


colnames.X <- paste("X",1:ncol(Xmat),sep="")

## Now, start split sample here ----

replaceme <- rep(1,n)
replaceme[1:floor(n/2)] <- 2
replaceme <- sample(replaceme)

## Partial out X's ----

y.partial <- partialOut(y, X, replaceme)
treat.partial <- partialOut(treat, X, replaceme)
treatmat.theta <- bs.me(treat.partial)
treatmat.tau <- dbs.me(treat.partial)
colnames.treat <- paste("treat",1:ncol(treatmat.theta),sep="")

## Calculate correlations

n1 <- sum(replaceme==1)
bases.obj <- namesAndCorrs(
  XSubsamp =  Xmat[replaceme==1,],
  # Xnames = colnames.X,
  ySubsamp = rank(y.partial[replaceme==1]),
  treatSubsamp = treatmat.theta[replaceme==1,],
  XConstruct = Xmat[replaceme==2,],
  treatConstruct = treatmat.theta[replaceme==2,],
  XConstructDerivative = Xmat[replaceme==2,],
  treatConstructDerivative = treatmat.tau[replaceme==2,],
  # treatNames = colnames.treat,
  a = ceiling(25*(1+n1^.2))
)

nc1<-ncol(bases.obj$Msubsamp)
bases.obj$Msubsamp <- bases.obj$Msubsamp[,nc1:1]
bases.obj$MConstruct <- bases.obj$MConstruct[,nc1:1]
bases.obj$MConstructDerivative <- bases.obj$MConstructDerivative[,nc1:1]

cormat <- t(matrix(unlist(bases.obj$cors),nrow=4))
cormat <- cormat[(nrow(cormat)):1,]
keeps <- which(as.vector(checkcor(apply(bases.obj$Msubsamp,2,rank), .9))==1)

bases.obj$Msubsamp <- bases.obj$Msubsamp[,keeps]
bases.obj$MConstruct <- bases.obj$MConstruct[,keeps]
bases.obj$MConstructDerivative <- bases.obj$MConstructDerivative[,keeps]


sp1<- sparsereg::sparsereg(y.partial[replaceme==1]+mean(y[replaceme==1]),bases.obj$Msubsamp,EM=T, verbose=F, iter.initialize=0)


beta.try <- sp1$coef[1,]
plot(treat[replaceme==1],sp1$fitted)
treat.seq <- seq(min(treat),max(treat),length=1000)
lines(treat.seq,treat.seq^2)

cor(cbind(0,bases.obj$MConstructDerivative)%*%beta.try,tau.true[replaceme==2])

plot(treat,treat*2,type="l")
points(treat[replaceme==2],cbind(0,bases.obj$MConstructDerivative)%*%beta.try)
lines(treat.seq,2*treat.seq)

