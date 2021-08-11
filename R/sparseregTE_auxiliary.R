####################################
### Auxiliary functions -------


## B spline functions -----
##  Make a bspline matrix from a vector 
bs.me<-function(x, xvar=TRUE){
  x <- x-mean(x)
  if(length(unique(x)) <= 2) return(as.matrix(x))
  if(length(unique(x)) <= 3) return(as.matrix(cbind(x,x^2)))
  if(length(unique(x)) <= 4) return(as.matrix(cbind(x,x^2,x^3)))
  
  m1<-cbind(x,bSpline2(x,df=3))#,bSpline2(x,df=4),bSpline2(x,df=5),bSpline2(x,df=6))
  return(m1)
}

##  Derivative of bspline from bs.me 
dbs.me<-function(x, xvar=TRUE){
  x<-x-mean(x)
  if(length(unique(x)) <= 2) return(x)
  if(length(unique(x)) <= 3) return(cbind(x,x^2))
  m1<-cbind(1,dbs2(x,df=3))#,dbs2(x,df=4),dbs2(x,df=5),dbs2(x,df=6))
  return(m1)

}

##  Making a single set of bases and derivative

bSpline2<-function(x,df,...){
  x<-x-mean(x)
  mx <- median(x)
  b1<-bSpline(x,knots=mx,degree=df,...)
  b2<- bSpline(-x,knots=mx,degree=df,...)
  knots2 <- quantile(x, seq(1:3)/4)
  b3<-bSpline(x,knots=knots2,degree=df,...)
  b4<-bSpline(-x,knots=knots2,degree=df,...)
  cbind(b2[,ncol(b2)],b1,b4[,ncol(b4)],b3)
}

dbs2<-function(x,df,...){
  x <- x-mean(x)
  
  mx <- median(x)
  b1<-dbs(x,knots=mx,degree=df,...)
  b2<- -dbs(-x,knots=mx,degree=df,...)
  #return(cbind(-b2[,ncol(b2)],b1))
  
  knots2 <- quantile(x, seq(1:3)/4)
  b3<-dbs(x,knots=knots2,degree=df,...)
  b4<- -dbs(-x,knots=knots2,degree=df,...)
  cbind(b2[,ncol(b2)],b1,b4[,ncol(b4)],b3)
}

## Hodges Lehmann mean

hl.mean <-function (x) 
{
  x <- x[!is.na(x)]
  diff1 <- as.vector(outer(x, x, "+")/2)
  median(c(diff1, x))
}

hl.var <-function (x) {
  hl.mean(x^2)-(hl.mean(x))^2
}

## Partial out X
partialOut <- function(y, X, replaceme){
  mod1 <- ranger(y~., data=data.frame(X), case.weights = 1*(replaceme==1), num.trees=1000,
                 write.forest = F)
  y.partialout <- y-mod1$predictions
  y.partialout - mean(y.partialout[replaceme==1])

}

createBases <- function(replaceme, Xmat, y.partial, treatmat.theta, treatmat.tau, ratio){
  
  n1 <- sum(replaceme==1)
  bases.obj <- namesAndCorrs(
    XSubsamp =  Xmat[replaceme==1,],
    # Xnames = colnames.X,
    ySubsamp = rank(y.partial[replaceme==1]),
    treatSubsamp = treatmat.theta[replaceme==1,],
    XConstruct = Xmat,
    treatConstruct = treatmat.theta,
    XConstructDerivative = Xmat,
    treatConstructDerivative = treatmat.tau,
    # treatNames = colnames.treat,
    a = ceiling(ratio*(1+n1^.2))
  )
  
  
  cormat <- (matrix(unlist(bases.obj$cors),nrow=4))
  keeps <- which(as.vector(checkcor(bases.obj$Msubsamp, .9))==1)
  if(length(keeps)<5)   keeps <- which(as.vector(checkcor(bases.obj$Msubsamp, .95))==1)
  if(length(keeps)<5)   keeps <- which(as.vector(checkcor(bases.obj$Msubsamp, .99))==1)
  
  bases.obj$cormat <- cormat[,keeps]
  
  bases.obj$Msubsamp <- bases.obj$Msubsamp[,keeps]
  bases.obj$MConstruct <- bases.obj$MConstruct[,keeps]
  bases.obj$MConstructDerivative <- bases.obj$MConstructDerivative[,keeps]
  
  bases.obj
  
}

## Primary fitting function
fit.singlesubsample <- function(y0, treat0, X0, replaceme0, Xmat0){
  ## Partial out X's ----
  y<-y0
  treat<-treat0
  X<-X0
  replaceme<-replaceme0
  Xmat<-Xmat0
  
  y.partial <- partialOut(y, X, replaceme)
  Ey.x <- y-y.partial
  treat.partial <- partialOut(treat, X, replaceme)
  treatmat.theta <- cbind(treat.partial,bs.me(treat.partial, xvar=FALSE))
  treatmat.tau <- cbind(1,dbs.me(treat.partial, xvar=FALSE))
  
  keeps <- which(as.vector(checkcor(treatmat.theta, .9)==1))
  treatmat.theta <- treatmat.theta[,keeps]
  treatmat.tau <- treatmat.tau[,keeps]
  
  colnames.treat <- paste("treat",1:ncol(treatmat.theta),sep="")
  
  ## Calculate correlations
  bases.obj <- createBases(replaceme, Xmat, y.partial, treatmat.theta, treatmat.tau,ratio=50)
  # repeat.SIS <- F
  # if(repeat.SIS){
  #  ste.EM0 <- sparsereg::sparsereg(y.partial[replaceme==1],bases.obj$Msubsamp,EM=T, verbose=F, iter.initialize=0, alpha.prior="balanced", use.sparseregweights=T, thresh=1e-4)
  # y.partial2 <- y.partial
  # y.partial2[replaceme==1] <- y.partial[replaceme==1]-ste.EM0$fitted.values
  # bases.obj2 <- createBases(replaceme, Xmat, y.partial2, treatmat.theta, treatmat.tau,ratio=200)
  # 
  # # for(i.join in 2:4) bases.obj[[i.join]] <- cbind(bases.obj[[i.join]][,(ste.EM0$coef[-1])!=0],bases.obj2[[i.join]])
  # for(i.join in 2:4) bases.obj[[i.join]] <- cbind(bases.obj[[i.join]], bases.obj2[[i.join]])
  # 
  # }
  
  #ste.EM <- sparsereg::sparsereg(y.partial[replaceme==1],bases.obj$Msubsamp,EM=T, verbose=F, iter.initialize=0, use.sparseregweights=F, thresh=1e-4)
  #beta.sp <- ste.EM$coef[1,]

  n.a <-sum(replaceme==1)
  p.a <-ncol(bases.obj$Msubsamp)
  
  # m1 <- mget(ls())
  # save(file="diagnose.rda",m1)
  alpha.seq <-seq(max(n.a*log(p.a),10*p.a), p.a, length=10)
  g1 <- GCV(y.partial[replaceme==1],cbind(1,bases.obj$Msubsamp), alphas=alpha.seq,tol=1e-2)
  beta.sp <- as.vector(g1$beta)
  
  
  te.curr<-cbind(0, bases.obj$MConstructDerivative)%*%beta.sp
  fits.curr<-cbind(1, bases.obj$MConstruct)%*%beta.sp
  
  ## Variance calculations ----
  
  # Variance of fitted value
  
  res.sq <- (y.partial-fits.curr)
  res.sq <- (res.sq-mean(res.sq[replaceme==1]))^2
  treat2<-treat
  var1<-ranger(res.sq~.,data=data.frame(treat2,X),case.weights = 1*(replaceme==1))
  
  numvar <- 25
  var.treatperm <- 0
  for(i in 1:numvar) {
    treat2<-sample(treat)
    var.treatperm <- var.treatperm + predict(var1, data=data.frame(treat2,X))[[1]]/numvar
  }
  # var.tau <- abs(var.treatperm[replaceme==2]-var1$predictions[replaceme==2])
  var.tau <- pmax(var.treatperm[replaceme==2]-var1$predictions[replaceme==2],0)
  output <- list("theta.pred" = fits.curr[replaceme==2], "tau.pred" = te.curr[replaceme==2],
                 "var.theta" = var1$predictions[replaceme==2], "var.tau"=var.tau,
                  "y.partial" = y.partial[replaceme==2], 
                 "Ey.x" = Ey.x[replaceme==2]
                 )
  
          
  return(output)
}


#' MDEI function
#' 
#' Implements the Method of Direct Estimation and Inference
#' @param y The outcome variable, a vector.
#' @param treat The treatment variable, a vector.
#' @param X A matrix of covariates.
#' @param splits Number of repeated cross-fitting steps to implement.
#' @param alpha The desired level of the confidence band.
#' @export


MDEI<-function(y,treat,X, splits=10, alpha = .9){
  n<-length(treat)
  treatmat <- cbind(treat,bs.me(treat))
  #Xmat<-cbind(1,apply(X,2,rank),matrix(apply(X,2,bs.me),nrow=n))
  
  Xmat.spline <- matrix(NA,nrow=n,ncol = ncol(X)*ncol(bs.me(rnorm(20)))+10 )
  for(i.X in 1:ncol(X)){
    col.start <- which(is.na(Xmat.spline[1,]))[1]
    bmat <- bs.me(X[,i.X])
    col.stop <- ncol(bmat)+col.start-1
    Xmat.spline[,col.start:col.stop]<-bmat
    }
  
  Xmat.spline <- Xmat.spline[,is.finite(Xmat.spline[1,])]
  Xmat <- cbind(1,X,Xmat.spline)
  keeps <- which(as.vector(checkcor(Xmat, .9))==1)
  Xmat <- cbind(1,Xmat[,keeps])
  
  colnames.X <- paste("X",1:ncol(Xmat),sep="")
  
  ## Containers ----
  
  Ey.x.run <- y.partial.run <- theta.run <- tau.run <- thetavar.run <- tauvar.run <- matrix(NA,nrow=n, ncol=splits )
  ## Now, start split sample here ----
  
  for(i.runs in 1:splits){
    replaceme <- rep(1,n)
    replaceme[1:floor(n/2)] <- 2
    replaceme <- sample(replaceme)
    
    singlefit.2 <- fit.singlesubsample(y, treat, X, replaceme, Xmat)
    singlefit.1 <- fit.singlesubsample(y, treat, X, 3-replaceme, Xmat)
    # 
    theta.run[replaceme==1,i.runs] <- singlefit.1$theta.pred
    tau.run[replaceme==1,i.runs] <- singlefit.1$tau.pred
    thetavar.run[replaceme==1,i.runs] <- singlefit.1$var.theta
    tauvar.run[replaceme==1,i.runs] <- singlefit.1$var.tau
    y.partial.run[replaceme==1,i.runs] <- singlefit.1$y.partial
    Ey.x.run[replaceme==1,i.runs] <- singlefit.1$Ey.x
    
    theta.run[replaceme==2,i.runs] <- singlefit.2$theta.pred
    tau.run[replaceme==2,i.runs] <- singlefit.2$tau.pred
    thetavar.run[replaceme==2,i.runs] <- singlefit.2$var.theta
    tauvar.run[replaceme==2,i.runs] <- singlefit.2$var.tau
    y.partial.run[replaceme==2,i.runs] <- singlefit.2$y.partial
    Ey.x.run[replaceme==2, i.runs] <- singlefit.2$Ey.x
    
    
  }
  
  se.theta <- (apply(thetavar.run,1,hl.mean)+apply(theta.run,1,hl.var))^.5
  ts.theta <- (y.partial.run-theta.run)/se.theta
  critical.value.theta <- quantile(abs(ts.theta),alpha)
  
  CIs.theta <- apply(y.partial.run,1,hl.mean)+critical.value.theta*cbind(-se.theta, se.theta)
  
  se.tau <-(apply(tauvar.run,1,hl.mean)+apply(tau.run,1,hl.var))^.5
  critical.value.tau <- (critical.value.theta^2+1)^.5
  CIs.tau <- rowMeans(tau.run)+critical.value.tau*cbind(-se.tau, se.tau)
  
  output <- list("tau.est"=apply(tau.run,1,hl.mean),
                 "CIs.tau"=CIs.tau,
                 "theta.est"=apply(theta.run,1,hl.mean),
                 "CIs.theta"=CIs.theta,
                 "critical.values" = list("theta"=critical.value.theta, "tau"=critical.value.tau),
                 "Ey.x"=rowMeans(Ey.x.run)
  )
  return(output)
}
