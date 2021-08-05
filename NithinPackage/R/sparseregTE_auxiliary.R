####################################
### Auxiliary functions -------


## B spline functions -----
##  Make a bspline matrix from a vector 
#' @export
bs.me<-function(x){
  x <- scale(rank(x))
  m2<-cbind(x,bSpline2(x,df=3),bSpline2(x,df=4))#,bSpline2(x,df=5),bSpline2(x,df=6))
  # m2<-cbind(x,bSpline2(x,df=3),bSpline2(x,df=4),
  #           bSpline2(-x,df=3),bSpline2(-x,df=4))#,bSpline2(x,df=5),bSpline2(x,df=6))
  # 
  return(m2)
  
}

##  Derivative of bspline from bs.me 
#' @export
dbs.me<-function(x){
  x<-scale(rank(x))
  m1<-cbind(1,dbs2(x,df=3),dbs2(x,df=4))#,dbs2(x,df=5),dbs2(x,df=6))
  # m1<-cbind(x,dbs(x,df=3),dbs(x,df=4),
  #           -dbs(-x,df=3),-dbs(-x,df=4))
  return(m1)
  
}

##  Making a single set of bases and derivative

bSpline2<-function(x,df,...){
  mx <- median(x)
  b1<-bSpline(x,knots=mx,degree=df,...)
  b2<-bSpline(-x,knots=mx,degree=df,...)
  cbind(b2[,ncol(b2)],b1)
}

dbs2<-function(x,df,...){
  mx <- median(x)
  b1<-dbs(x,knots=mx,degree=df,...)
  b2<-dbs(-x,knots=mx,degree=df,...)
  cbind(-b2[,ncol(b2)],b1)
}

## Partial out X
partialOut <- function(y, X, replaceme){
  mod1 <- ranger(y~., data=data.frame(X), case.weights = 1*(replaceme==1), num.trees=1000,
                 write.forest = F)
  y.partialout <- y-mod1$predictions
  y.partialout - mean(y.partialout[replaceme==1])

}

