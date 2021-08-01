library(splines2)
library(ranger)
library(microbenchmark)
if(F){
  Rcpp::sourceCpp('Dropbox/Github/warmup/NithinPackage/src/splinesinter_short.cpp')
  Rcpp::sourceCpp('Dropbox/Github/warmup/NithinPackage/src/bayesLassoRevised.cpp')
  
  n<-200
  p<-5
  
  X <- matrix(rnorm(n*p),nrow=n)
  X<-apply(X,2,rank)
  X<- apply(X,2,scale)
  colnames(X) <- paste("X",1:ncol(X),sep="_")
  treat<- rnorm(n)
  
  y <- treat^2 + X[,3]+treat*X[,2]+rnorm(n)
  

  
  }



## First, turn covariates and treatment into spline bases, save these ----
n<-length(treat)
treatmat <- bs.me(treat)
Xmat<-matrix(apply(X,2,bs.me),nrow=n)

colnames.treat <- paste("treat",1:ncol(treatmat),sep="")
colnames.X <- paste("X",1:ncol(Xmat),sep="")

## Now, start split sample here ----

replaceme <- rep(1,n)
replaceme[1:floor(n/2)] <- 2
replaceme <- sample(replaceme)



## Partial out X's ----

y.partial <- partialOut(y, X, replaceme)
treat.partial <- partialOut(treat, X, replaceme)
treatmat.partial <- bs.me(treat.partial)
## Calculate correlations

#splineCorrs(arma::mat XSubsamp, arma::vec ySubsamp, arma::mat treatSubsamp, arma::mat XConstruct, arma::mat treatConstruct,  long long unsigned int a)
splineCorrs.temp <- splineCorrs(
  XSubsamp = Xmat[replaceme==1,],
  ySubsamp = y.partial[replaceme==1],
  treatSubsamp = treatmat.partial[replaceme==1,],
  XConstruct = Xmat[replaceme==2,],
  treatConstruct = treatmat.partial[replaceme==2,],
  a = 100
)


####################################
### Auxiliary functions -------


## B spline functions -----
##  Make a bspline matrix from a vector 
#' @export
bs.me<-function(x){
  m2<-cbind(x,bSpline2(x,df=3),bSpline2(x,df=4))#,bSpline2(x,df=5),bSpline2(x,df=6))
  return(m2)
  
}

##  Derivative of bspline from bs.me 
#' @export
dbs.me<-function(x){
  m1<-cbind(1,dbs2(x,df=3),dbs2(x,df=4))#,dbs2(x,df=5),dbs2(x,df=6))
  return(m1)
  
}

##  Making a single set of bases and derivative

bSpline2<-function(x,...){
  mx <- median(x)
  b1<-bSpline(x,knots=mx,...)
  b2<-bSpline(-x,knots=mx,...)
  cbind(b2[,ncol(b2)],b1)
}

dbs2<-function(x,...){
  mx <- median(x)
  b1<-dbs(x,knots=mx,...)
  b2<-dbs(-x,knots=mx,...)
  cbind(-b2[,ncol(b2)],b1)
}

## Partial out X
partialOut <- function(y, X, replaceme){
  mod1 <- ranger(y~., data=data.frame(X), case.weights = 1*(replaceme==1), num.trees=1000,
                 write.forest = F)
  y.partialout <- y-mod1$predictions
  y.partialout - mean(y.partialout[replaceme==1])

}

