####################################
### Auxiliary functions -------


## B spline functions -----
##  Make a bspline matrix from a vector
bs.me <- function(x, varname) {
  x <- x - mean(x)
  if (length(unique(x)) <= 2){
    m1 <- as.matrix(x)
    colnames(m1)<-paste(varname,"linear",sep="_")
    return(m1)
  }
  if (length(unique(x)) <= 3){
    m1<-as.matrix(cbind(x, x ^ 2))
    colnames(m1)<-paste(varname,c("linear","quadratic"),sep="_")
    return(m1)
  }
  if (length(unique(x)) <= 4){
    m1<-cbind(x, x ^ 2, x ^ 3)
    colnames(m1)<-paste(varname,c("linear","quadratic","cubic"),sep="_")
    return(m1)
  }
  
  b1 <- bSpline2(x, df = 3)
  colnames(b1) <- paste(varname,"spline",1:ncol(b1),sep="_")
  m1 <-
    cbind(x, b1)#,bSpline2(x,df=4),bSpline2(x,df=5),bSpline2(x,df=6))
  colnames(m1)[1]<-paste(varname,"linear",sep="_")
  return(m1)
}

##  Derivative of bspline from bs.me
dbs.me <- function(x) {
  x <- x - mean(x)
  if (length(unique(x)) <= 2)
    return(as.matrix(1))
  if (length(unique(x)) <= 3)
    return(cbind(1, 2*x))
  m1 <- cbind(1, dbs2(x, df = 3))#,dbs2(x,df=4),dbs2(x,df=5),dbs2(x,df=6))
  return(m1)
  
}

##  Making a single set of bases and derivative

bSpline2 <- function(x, df, ...) {
  x <- x - mean(x)
  mx <- median(x)
  b1 <- bSpline(x, knots = mx, degree = df, ...)
  b2 <- bSpline(-x, knots = mx, degree = df, ...)
  knots2 <- quantile(x, seq(1:3) / 4)
  b3 <- bSpline(x, knots = knots2, degree = df, ...)
  b4 <- bSpline(-x, knots = knots2, degree = df, ...)
  cbind(b2[, ncol(b2)], b1, b4[, ncol(b4)], b3)
}

dbs2 <- function(x, df, ...) {
  x <- x - mean(x)
  
  mx <- median(x)
  b1 <- dbs(x, knots = mx, degree = df, ...)
  b2 <- -dbs(-x, knots = mx, degree = df, ...)
  #return(cbind(-b2[,ncol(b2)],b1))
  
  knots2 <- quantile(x, seq(1:3) / 4)
  b3 <- dbs(x, knots = knots2, degree = df, ...)
  b4 <- -dbs(-x, knots = knots2, degree = df, ...)
  cbind(b2[, ncol(b2)], b1, b4[, ncol(b4)], b3)
}

## Hodges Lehmann mean

hl.mean <- function (x)
{
  x <- x[!is.na(x)]
  diff1 <- as.vector(outer(x, x, "+") / 2)
  median(c(diff1, x))
}

hl.var <- function (x) {
  hl.mean(x ^ 2) - (hl.mean(x)) ^ 2
}

## Partial out X
partialOut <- function(y, X, replaceme) {
  mod1 <-
    ranger(
      y ~ .,
      data = data.frame(X),
      case.weights = length(y) * (replaceme == 3),
      num.trees = 1000,
      write.forest = F
    )
  
    preds1 <- mod1$predictions
    y - preds1
  
}

createBases <-
  function(replaceme,
           Xmat,
           y.partial,
           treatmat.theta,
           treatmat.tau,
           ratio) {
    
    n1 <- sum(replaceme == 1)
    bases.obj <- namesAndCorrs(
      XSubsamp =  Xmat[replaceme == 3, ],
      ySubsamp = rank(lm(y.partial[replaceme == 3]~treatmat.theta[replaceme == 3,1])$res),
      treatSubsamp = treatmat.theta[replaceme == 3, ],
      XConstruct = Xmat,
      treatConstruct = treatmat.theta,
      XConstructDerivative = Xmat,
      treatConstructDerivative = treatmat.tau,
      a = ceiling(min(ratio * (1 + n1 ^ .2), n1/4))
    )
    
    ## Put treatment vector and intercept in  ----
    bases.obj$Msubsamp <- cbind(treatmat.theta[replaceme == 3, 1],bases.obj$Msubsamp)
    bases.obj$MConstruct <- cbind(treatmat.theta[, 1],bases.obj$MConstruct)
    bases.obj$MConstructDerivative <-
      cbind(1,bases.obj$MConstructDerivative)
    
    ## Eliminate correlated rows ----
    cormat <- (matrix(unlist(bases.obj$cors), nrow = 4))
    keeps <- which(as.vector(checkcor(bases.obj$Msubsamp, .9)) == 1)
    if (length(keeps) < 5)
      keeps <- which(as.vector(checkcor(bases.obj$Msubsamp, .95)) == 1)
    if (length(keeps) < 5)
      keeps <- which(as.vector(checkcor(bases.obj$Msubsamp, .99)) == 1)
    
    bases.obj$cormat <- cbind(0,cormat)[, keeps]
    bases.obj$Msubsamp <- bases.obj$Msubsamp[, keeps]
    bases.obj$MConstruct <- bases.obj$MConstruct[, keeps]
    bases.obj$MConstructDerivative <-
      bases.obj$MConstructDerivative[, keeps]
    
    return(bases.obj)

  }


