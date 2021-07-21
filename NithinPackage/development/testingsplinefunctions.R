library(microbenchmark)
library(sparsereg)
library(splines2)
n<-200
p<-5

X <- matrix(rnorm(n*p),nrow=n)
X<-apply(X,2,rank)
X<- apply(X,2,scale)

treat<- rnorm(n)

y<- treat^2+treat*X[,1]+rnorm(n,sd=2)

###Spline functions!

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


##  Make a bspline matrix from a vector ----
#' @export
bs.me<-function(x){
  m2<-cbind(x,bSpline2(x,df=3),bSpline2(x,df=4))#,bSpline2(x,df=5),bSpline2(x,df=6))
  return(m2)
  
}

##  Derivative of bspline from bs.me ----
#' @export
dbs.me<-function(x){
  m1<-cbind(1,dbs2(x,df=3),dbs2(x,df=4))#,dbs2(x,df=5),dbs2(x,df=6))
  return(m1)
  
}

treatmat <- bs.me(treat)


Xmat<-matrix(apply(X,2,bs.me),nrow=n)

colnames.treat <- paste("treat",1:ncol(treatmat),sep="")
colnames.X <- paste("X",1:ncol(Xmat),sep="")
# splineBasesAndCorrs(arma::mat XSubsamp, 
# std::vector<std::string> Xname, 
# arma::vec ySubsamp, 
# std::vector<int> colSizes, 
#arma::mat treatSubsamp, 
# arma::mat XConstruct, arma::mat treatConstruct, 
#std::string treatName, long long unsigned int a) {
 spline1 <- splineCorrs(XSubsamp=as.matrix(Xmat),
                                #Xname=colnames.X,
                                ySubsamp=rank(y),
                                #colSizes = cumsum(rep(23,ncol(X))),
                                treatSubsamp = as.matrix(treatmat),
                                XConstruct=as.matrix(Xmat),
                                treatConstruct=as.matrix(treatmat),
                                #treatName=colnames.treat,
                                a=100) 
 
 cormat <- matrix(unlist(spline1$cors),ncol=4,byrow=T)


microbenchmark( splineCorrs(XSubsamp=as.matrix(Xmat),
                                        #Xname=colnames.X,
                                        ySubsamp=y,
                                        #colSizes = cumsum(rep(23,ncol(X))),
                                        treatSubsamp = as.matrix(treatmat),
                                        XConstruct=as.matrix(Xmat),
                                        treatConstruct=as.matrix(treatmat),
                                        #treatName=colnames.treat,
                                        a=100) , times=10
)


sc1 <- splineCorrs(XSubsamp=as.matrix(Xmat),
                       #Xname=colnames.X,
                       ySubsamp=y,
                       #colSizes = cumsum(rep(23,ncol(X))),
                       treatSubsamp = as.matrix(treatmat),
                       XConstruct=as.matrix(Xmat),
                       treatConstruct=as.matrix(treatmat),
                       #treatName=colnames.treat,
                       a=100) 
  
p <- ncol(sc1$M)
alpha.median <- (n*log(p)*(p/n)+(p)*(n/p))/(p/n+n/p)*sd(y)

alpha.seq<-seq(alpha.median*10,alpha.median/10,length=10)


microbenchmark(GCV(y,cbind(1,sc1$M),alpha.seq,1e-5),times=10)

 