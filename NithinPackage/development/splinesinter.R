library(splines2)
## Source functions
source('C:\\Users\\think\\OneDrive\\Documents\\NithinPackage\\R\\sparseregTE_split.R')

## Will need splines2, randomForest, coop, RcppEigen, Rfast

### Set up data
n<-1000 #n obs
p<-5 #n covs

## Treatment variable
## This is the one for which we want to estimate the first derivative.

treat<-rnorm(n)

## Covariates
X<-matrix(rnorm(n*p),n,p)

## Outcome
## Note that there is no function of only the 
## X's in here; we handle that elsewhere
y <- treat*X[,1]^2+rnorm(n)


## OK--your work starts here
## First, we are going toconvert
## The treatment vector and X's into separate matrices
## We are going to construct functions of the form
## E(y|t,x) = \phi_i(t) \phi_j(x_k) \phi_l(x_m)
## And keep track of its partial derivative
##  \phi^\prime_i(t) \phi_j(x_k) \phi_l(x_m)

## Make treatment bases and derivatives
treat.bs <- bs.me(treat)
treat.dbs <- dbs.me(treat)

par(mfcol=c(1,2),mar=c(1,1,1,1))
plot(treat,treat.bs[,4],main="spline")
plot(treat,treat.dbs[,4],main="derivative")


## Make matrix of spline bases

X.bs <- NULL
for(i in 1:p){
  X.bs <- cbind(X.bs,make.bs(X[,i]))
}

## Find the most correlated interactions.
## Note that this code has lots of duplications, 
## please don't do that in your code.  This code is meant to be simple.
cor.output <- NULL

## Fit the correlation off half the data.
subsamp <- sample(1:n, n/2, FALSE)

for(i.t in 1:ncol(treat.bs)){
  for(i.1 in 1:ncol(X.bs)){
    for(i.2 in i.1:ncol(X.bs)){
      index.curr<-c(i.t,i.1,i.2)
      inter.temp <- treat.bs[subsamp,i.t]*X.bs[subsamp,i.1]*X.bs[subsamp,i.2]
      cor.temp <- cor(y[subsamp],inter.temp)
      output.temp <- c(index.curr,cor.temp)
      cor.output <- rbind(cor.output,output.temp)
    }
  }
}

## We are going to work with the most correlated bases.
## We will get to this when I get back, but we are going
## To estimate the most correlated bases off

## One additional function to think about if time, gram schmidt.
## It is going to take an outcome y and a matrix X.
## It finds the most correlated column of X with y.
## places this column in the first column of the output matrix then uses
## least squares to partial out of the remaining columns of X.
## It then does the same thing to the remaining columns, in term.
## Code (messy!) is below.  Don't get lost in this!
gramschmidt<-function(y,X,X.cor=NULL, weights.lm){
  ## X, X.cor are the data
  ## X2, X2.cor will be populated
  
  if(is.null(X.cor)) X.cor<-X
  n<-length(y)
  
  ## Strip out covariats with no covariatnce
  keeps<-apply(X,2,sd)*(apply(X.cor,2,sd))!=0
  keeps[is.na(keeps)]<-F
  X.cor<-X.cor[,keeps]
  X<-X[,keeps]
  X.cor<-apply(X.cor,2,scale)
  X<-apply(X,2,scale)
  y<-as.vector(scale(y))
  X2.cor<-X2<-X*NA
  ng<-min(ncol(X),floor(n*.9))
  for(i in 1:ng)
  {
    cors.curr<-abs(t(X.cor)%*%y)
    ind.curr<-which(cors.curr==max(cors.curr))[1]
    x.curr<-X[,ind.curr]
    x.curr.cor<-X.cor[,ind.curr]
    if(i>1) {
      x.curr<-scale(lm(x.curr~X2[,1:(i-1)])$res)
      x.curr[is.na(x.curr)]<-0
      X<-apply(X,2,FUN=function(x) scale(lm(x~x.curr)$res))
      
      x.curr.cor<-scale(lm(x.curr.cor~X2.cor[,1:(i-1)])$res)
      X.cor<-apply(X.cor,2,FUN=function(x) scale(lm(x~x.curr.cor)$res))
      
      X[is.na(X)]<-0; X[!is.finite(X)]<-0
      X.cor[is.na(X.cor)]<-0; X.cor[!is.finite(X.cor)]<-0
    }
    X.cor[,ind.curr]<-X[,ind.curr]<-0
    
    X2[,i]<-x.curr
    X2.cor[,i]<-x.curr.cor
    
    if(i>1) y<-scale(lm(y~X2.cor[,1:(i-1)])$res)
    
  }
  X2<-X2[,colMeans(is.na(X2))==0]
  X2
}



#####  Reflating bases
##  Make treatment effect bases w/o checking correlations ----
make.justbaseste<-function(rest,resy,X,inter.schedule,tspline,tspline.del,Xte,ste, quants,lms,X.lin=NULL){
  n<-length(rest)
  ## B-spline of treatment ----
  tspline<-cbind(rest,bs.me(rest))
  tspline.del<-cbind(rest*0+1,dbs.me(rest))
  nc<-ncol(tspline)-1
  colnames(tspline.del)<-colnames(tspline)<-c("treatlin",paste("treatspline",1:nc,"_"))
  #paste("treatbs3",1:3,sep= "_"),paste("treatbs5",1:5,sep= "_"))
  
  baseste<-cbind(as.vector(rest),apply(inter.schedule,1,construct.cort,tspline0=tspline,Xt0=Xte))
  baseste.del<-cbind(1,apply(inter.schedule,1,construct.cort,tspline0=tspline.del,Xt0=Xte))
  colnames(baseste.del)<-colnames(baseste)<-c(
    "treat",#paste("treat",colnames(X2),sep=":"),
    apply(inter.schedule,1,names.cort,tspline0=tspline,Xt0=Xte)   
  )
  
  keeps<-check.cor(baseste,thresh=1e-5,nruns=5)$k
  output<-list("baseste.del"=baseste.del[,keeps],"baseste"=baseste[,keeps])
  return(output)
}

##  For removing correlated columns ----
check.cor<-function(X,thresh=0,nruns=3){
  if(length(thresh)==0) thresh<-ifelse(ncol(X)>100,0.001,0.0001)
  drops.run<-rep(0,ncol(X))
  for(j in 1:nruns){
    drops<-rep(0,ncol(X))
    cor1<-as.vector(abs(cor(X,rnorm(nrow(X)))))
    for(i in 1:(ncol(X)-1)){
      cor.temp<-cor1[(i+1):ncol(X)]
      drops[(i+1):ncol(X)][abs(cor.temp-cor1[i])<thresh]<- -1
    }
    drops.run<-drops+drops.run
  }
  keeps<-drops.run!=(-nruns)
  out1<-list("keeps"=keeps)
  return(out1)	
}