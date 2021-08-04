library(splines2)
library(ranger)
library(microbenchmark)
if(F){
  Rcpp::sourceCpp('~/Dropbox/Github/warmup/NithinPackage/src/splinesinter_short.cpp')
  Rcpp::sourceCpp('~/Dropbox/Github/warmup/NithinPackage/src/bayesLassoRevised.cpp')
  source('~/Dropbox/Github/warmup/NithinPackage/R/sparseregTE_auxiliary.R')
  n<-2000
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

# #splineCorrs(arma::mat XSubsamp, arma::vec ySubsamp, arma::mat treatSubsamp, arma::mat XConstruct, arma::mat treatConstruct,  long long unsigned int a)
# splineCorrs.temp <- splineCorrs(
#   XSubsamp = Xmat[replaceme==1,],
#   ySubsamp = y.partial[replaceme==1],
#   treatSubsamp = treatmat.partial[replaceme==1,],
#   XConstruct = Xmat[replaceme==2,],
#   treatConstruct = treatmat.partial[replaceme==2,],
#   a = 50
# )
# 
# splinebases.temp <- splineCorrs.temp$M

#List namesAndCorrs(arma::mat XSubsamp,
# std::vector<std::string> Xnames, 
# arma::vec ySubsamp, 
# std::vector<int> colSizes, 
# arma::mat treatSubsamp, 
# arma::mat XConstruct, 
# arma::mat treatConstruct, 
# arma::mat XConstructDerivative, 
# arma::mat treatConstructDerivative, 
# std::vector<std::string> treatNames, 
# long long unsigned int a) 
  

bases.obj <- namesAndCorrs(
  XSubsamp =  Xmat[replaceme==1,],
  # Xnames = colnames.X,
  ySubsamp = y.partial[replaceme==1],
  treatSubsamp = treatmat.theta[replaceme==1,],
  XConstruct = Xmat[replaceme==2,],
  treatConstruct = treatmat.theta[replaceme==2,],
  XConstructDerivative = Xmat[replaceme==2,],
  treatConstructDerivative = treatmat.tau[replaceme==2,],
  # treatNames = colnames.treat,
  a = 100
)

keeps <- as.vector(checkcor(bases.obj$Msubamp, 0.95))

bases.obj$Msubamp <- bases.obj$Msubamp[,keeps]
bases.obj$MConstruct <- bases.obj$MConstruct[,keeps]
bases.obj$MConstructDerivative <- bases.obj$MConstructDerivative[,keeps]

maxalpha <- n*log(ncol(bases.obj$Msubamp))*5
minalpha <- p
alpha.seq <- sequence((maxalpha),(minalpha),length=10)
g1<-GCV(y.partial[replaceme==1],cbind(1,bases.obj$Msubamp), alpha.seq, tol=sd(y)*1e-6)

beta.try <- sparsereg::sparsereg(y.partial[replaceme==1],bases.obj$Msubamp,EM=T, verbose=F, iter.initialize=20)

cor(cbind(0,bases.obj$MConstructDerivative)%*%g1$beta,tau.true[replaceme==2])
cor(cbind(0,bases.obj$MConstructDerivative)%*%beta.try$coefficients[1,],tau.true[replaceme==2])

plot(cbind(0,bases.obj$MConstructDerivative)%*%g1$beta,tau.true[replaceme==2])

plot(cbind(0,bases.obj$MConstructDerivative)%*%g1$beta,tau.true[replaceme==2])

