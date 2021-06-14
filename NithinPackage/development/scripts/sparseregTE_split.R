### Required packages  -----
library(splines2)
library(randomForest)
library(MASS)
library(coop)
library(RcppEigen)
library(Rfast)

### Some Rcpp functions, faster than in R ----
if(!exists("arma_cor")){
  arma_cor_code <- 
    "arma::vec arma_mm(const arma::vec& v1, const arma::vec& v2 ) {
       return arma::cor(v1 , v2 );
   };"
  arma_cor = cppFunction(code = arma_cor_code, depends = "RcppArmadillo")
  
  arma_var_code <- 
    "arma::vec arma_mm(const arma::vec& v1) {
       return arma::cov(v1);
   };"
  arma_var = cppFunction(code = arma_var_code, depends = "RcppArmadillo")
  
  arma_sd_code <- 
    "arma::vec arma_mm(const arma::vec& v1) {
       return arma::pow(arma::cov(v1),.5);
   };"
  arma_sd = cppFunction(code = arma_sd_code, depends = "RcppArmadillo")
  
  
  arma_med_code <- 
    "arma::vec arma_mm(const arma::mat& v1) {
       return arma::median(v1);
   };"
  arma_med = cppFunction(code = arma_med_code, depends = "RcppArmadillo")
  
  arma_invert_code <- 
    "arma::mat arma_mm(const arma::mat& v1) {
       return arma::inv(v1);
   };"
  arma_invert = cppFunction(code = arma_invert_code, depends = "RcppArmadillo")
  
  arma_fastcor_code <- 
    "arma::mat arma_mm(const arma::mat& indices, const arma::mat& treatbases, const arma::mat& xbases, 
    const arma::mat& res, const arma::vec& output) {
      
       arma::mat tspline = treatbases.col(indices(0,0))%
                            xbases.col(indices(1,0))%
                             xbases.col(indices(2,0));

      arma::mat cor_temp = arma::cor(tspline,res);
      output(0) = cor_temp;
      ;
   };"
#arma_fastcor = cppFunction(code = arma_fastcor_code, depends = "RcppArmadillo")
  
 
}

### Main sparseregTE function ----
#' update function
#' @export
sparseregTE <- function(y, treat, X, X.lin=NULL, id=NULL,weights = 1, splitsamp=T, 
                        interpretable=TRUE, gs=FALSE,splits=20, CI=.9, 
                        partial.out.treat = TRUE,critical.value=NULL) {
  
  time.keep<-list()
  time.keep[[1]]<-proc.time()
  
  ## Format matrices ----
  X<-cleanX(apply(X,2,rank), interpretable)
  n<-length(y)
  #sd.treat<-sd(treat)
  #treat<-as.vector(scale(treat))
  sd.treat<-diff(range(treat))
  treat<-(treat-min(treat))/sd.treat
  
  ## Format bases  ----
  X2<-X
  y2<-y
  treat2<-treat
  if(length(X.lin)>0){
    X2<-apply(X,2,FUN=function(x) lm(x~X.lin)$res)
    y2<-lm(y~X.lin)$res
    treat2<-lm(treat~X.lin)$res
  }
  
  rf1<-randomForest(y2,x=X2)
  resy.full<-resy<-y2-predict(rf1)
  
  if(partial.out.treat){
    rf1<-randomForest(treat2,x=X2)
    rest.full<-rest<-treat-predict(rf1)
  } else{
    rest.full<-rest<-treat
  }
  
  rsqt.run<-lms<-quants<-NULL
  te.bases.obj<-isis.te(rest,X2,resy,quants,lms)
  b1<-te.bases.obj$b1
  baseste<-b1$baseste
  baseste.del<-b1$baseste.del
  keeps<-check.cor(baseste,thresh = 1e-4)$k
  baseste.full<-baseste<-baseste[,keeps]
  baseste.del.full<-baseste.del<-baseste.del[,keeps]
  means.te<-te.bases.obj$means.te
  sds.te<-te.bases.obj$sds.te
  ste<-b1$ste
  Xte<-te.bases.obj$Xte
  full.inter.schedule<-te.bases.obj$full.inter.schedule
  
  tspline<-te.bases.obj$tspline
  tspline.del<-te.bases.obj$tspline.del
  
  b1<-make.justbaseste(rest,resy,X2,inter.schedule=full.inter.schedule,tspline,tspline.del,Xte,ste,quants,lms,X.lin=NULL)
  
  
  
  message("Starting split sample")
  ## Format containers ----
  replaceme.run<-treat.fits.all.run<-resy.run<-y.fits.run<-se.y.run<-t.fits.run<-te.fits.run<-se.run<-te.run<-matrix(NA,nrow=n,ncol=splits)
  betas.temp<-ste$coef/sd.treat
  coef.run<-matrix(0,nrow=1+ncol(b1$baseste),ncol=splits)
  rownames(coef.run)<-c("(Intercept)",colnames(b1$baseste))
  beta.init<-NULL
  
  time.keep[[2]]<-proc.time()
  
  ## Loop ----
  for(i.split in 1:splits){
    if(i.split%%2 == 1) replaceme<-make.replaceme(y,id) else replaceme<-3-replaceme
    
    ## Partial out X.lin
    y2<-y;X2<-X;treat2<-treat
    if(length(X.lin)>0){
      X2<-apply(X,2,FUN=function(x) lm(x~X.lin,weights=1*(replaceme==1))$res)
      y2<-lm(y~X.lin,weights=1*(replaceme==1))$res
      treat2<-lm(treat~X.lin,weights=1*(replaceme==1))$res
    }
    
    partialoutRF<-function(y,X=X2,replaceme0=replaceme,nt=2000){
      rf1<-randomForest(y,x=X,subset=(replaceme0==1),ntrees=nt)
      y-predict(rf1,newdata=X)
    }
    
    
    if(splitsamp) resy<-partialoutRF(y2) else resy<-partialoutRF(y2,replaceme0 = replaceme*0+1)
    if(partial.out.treat) {
      if(splitsamp) rest<-partialoutRF(treat2) else rest<-partialoutRF(treat2,replaceme0 = replaceme*0+1)
    } else {
      rest<-treat2
    }
    
    if(!splitsamp) replaceme<-replaceme*0+2
    rsqt.run<-c(rsqt.run,cov((treat-rest)[replaceme==2],treat[replaceme==2])/var(treat[replaceme==2]))
    
    tspline<-as.matrix(rest)
    tspline.del<-0*tspline+1
    
    b1<-make.justbaseste(rest,resy,X2,inter.schedule=full.inter.schedule,tspline,tspline.del,Xte,ste,quants,lms,X.lin=NULL)
    baseste<-b1$baseste
    baseste.del<-b1$baseste.del
    
    sds.run<-1
    
    for(i.r in 1:2){
      if(!splitsamp)  replaceme<-0*replaceme+i.r
      means<-colMeans(b1$baseste[replaceme==i.r,])
      sds<-apply(b1$baseste[replaceme==i.r,],2,sd)
      sds.run<-sds.run*sds
      for(i.scale in 1:ncol(baseste)){
        baseste[replaceme==i.r,i.scale]<-(baseste[replaceme==i.r,i.scale]-means[i.scale])/sds[i.scale]
        baseste.del[replaceme==i.r,i.scale]<-(baseste.del[replaceme==i.r,i.scale]-means[i.scale])/sds[i.scale]
      }
    }
    drops<-(sds.run==0)|is.na(sds.run)
    baseste<-baseste[,!drops]
    baseste.del<-baseste.del[,!drops]
    
    
    if(!splitsamp) replaceme<-replaceme*0+1
    screen1<-sparsereg(resy[replaceme==1],baseste[replaceme==1,],EM=T,tol=0.001,verbose=FALSE,unpen=2, beta.init=NULL,alpha.prior = "oracle",sparseregweights = TRUE, iter.initialize = 20)
    keep.cols<-colnames(baseste)%in%colnames(screen1$coefficients)
    resy2<-resy-cbind(1,baseste)%*%screen1$coefficients[1,]
    lm.adj<-lm(resy[replaceme==1]~screen1$fitted.values)
    beta.init<-screen1$coef
    screen2<-screen1
    
    keeps<-(abs(screen1$coef[1,][-1])+abs(screen1$coef[1,][-1]))>0.001*sd(resy)
    if(sum(keeps)<=1) {
      cor.ind<-sort(abs(cor(resy[replaceme==1],baseste[replaceme==1,-1])),dec=T,ind=T)$ix[1]+1
      keeps[cor.ind]<-TRUE
    }
    keeps[1]<-TRUE
    #errs.stage1<-resy2-screen2$fitted.values
    baseste<-as.matrix(baseste[,keeps])
    baseste.del<-as.matrix(baseste.del[,keeps])
    #baseste<-apply(baseste,2,partialoutRF,nt=50)
    
    if(!splitsamp)  replaceme<-0*replaceme+2
    ste.EM<-sparsereg(resy[replaceme==2],baseste[replaceme==2,],EM=T,verbose=F, tol = 0.001,id=NULL, unpen=2, sparseregweights = TRUE, iter.initialize = 20)
    ste.EM.all<-as.vector(cbind(1,baseste)%*%ste.EM$coef[1,])
    
    te.curr<-cbind(0, baseste.del)%*%(ste.EM$coef[1,])
    fits.curr<-cbind(1, baseste)%*%(ste.EM$coef[1,])
    
    ## Calculate ses ----
    XprimeX<-crossprod(ste.EM$X)
    errs.te<-(resy-ste.EM.all)#[replaceme==2]
    X.te<-cbind(0,baseste.del)#ste.EM$X#[replaceme==2,]
    
    
    if(!splitsamp)  replaceme<-0*replaceme+1
    boot.samp<-which(replaceme==1)
    rf.errX<-rf.errTX<-randomForest((errs.te[boot.samp])^2,x=cbind(treat,X)[boot.samp,])
    rf.errX<-randomForest((errs.te[boot.samp])^2,x=cbind(X)[boot.samp,])
    TXpred<-data.frame(treat,X)
    Xpred<-data.frame(sample(treat),X)
    names(Xpred)<-names(TXpred)
    
    if(!splitsamp) replaceme<-replaceme*0+2
    
    errs.diff<-predict(rf.errTX,newdata=TXpred)-predict(rf.errX,newdata=Xpred)
    errs.diff<-abs(errs.diff)
    errs.y<-abs(predict(rf.errTX,newdata=TXpred))
    
    ## Calculate LOO errors ----
    if(!splitsamp) replaceme<-replaceme*0+1
    
    taus.use<-ste.EM$parameters$tausq
    #Dtau.inv<-Dtau<-diag(ste.EM$parameters$tausq)
    Dtau.inv<-Dtau<-diag(taus.use)
    diag(Dtau.inv) <- 1/diag(Dtau)
    A <- arma_invert(XprimeX + Dtau.inv)

    if(!splitsamp) replaceme<-replaceme*0+2
    
    hii<-colSums(t(ste.EM$X)*(A%*%t(ste.EM$X)))
    errs.hii<-errs.te[replaceme==2]#/(1-hii)^.5
    
    ## Gather things ----
    if(!splitsamp) replaceme<-replaceme*0+2
    
    ses.curr<-(errs.diff[replaceme==2])^.5
    
    se.run[replaceme==2,i.split]<-ses.curr
    te.run[replaceme==2,i.split]<-te.curr[replaceme==2]/sd.treat
    y.fits.run[replaceme==2,i.split]<-(y2+fits.curr)[replaceme==2]-errs.hii
      #(y2-resy+fits.curr)[replaceme==2]
    #resy.run[replaceme==2,i.split]<-(resy)[replaceme==2]
    resy.run[replaceme==2,i.split]<-errs.hii#(resy)[replaceme==2]
    
    se.y.run[replaceme==2,i.split]<-(errs.y[replaceme==2])^.5
    
    replaceme.run[,i.split]<-replaceme
    treat.fits.all.run[,i.split]<-fits.curr
    
    betas.temp<-ste.EM$beta/sd.treat
    names(ste.EM$coef)
    coef.run[colnames(ste.EM$coef),i.split]<-ste.EM$coef[1,]/sd.treat
    treat.fits.all.run[replaceme==2,i.split]<-fits.curr[replaceme==2]
    
    if(i.split%%5==1) cat(c("Split", i.split, " of ", splits, ". \n"))
  }
  time.keep[[3]]<-proc.time()
  
  ## End of loop ----
  
  ses.out<-(apply(se.run^2,1,mean,na.rm=T)+apply(te.run,1,var,na.rm=T))^.5
  te.out<-apply(te.run,1,mean,na.rm=T)
  fits.out<-rowMeans(y.fits.run,na.rm=T)
  y.se.run<-(rowMeans(se.y.run^2,na.rm=T)+apply(y.fits.run,1,var,na.rm=T))^.5
  te.sd.out<-apply(te.run,1,sd,na.rm=T)
  ##find output
  alpha.cv<-CI



  ## Invert the confidence interval ----
  ts.all<-abs(y.fits.run-y)/y.se.run
  cvalpha<-apply(ts.all,2,quantile,CI,na.rm=T)
  ts.all<-as.vector(ts.all)
  ts.all<-ts.all[!is.na(ts.all)]
  cvalpha<-quantile(ts.all,CI)
    
  cvalpha2<-(mean(cvalpha)+1)
  if(length(critical.value)>0) cvalpha2<-critical.value
  
  CIs_old<-CIs<-te.out+cvalpha2*cbind(-ses.out,ses.out)
  diff(colMeans(CIs))
  
  
  ses.out.samp<-(apply(se.run^2,1,mean,na.rm=T))^.5
  ses.out.err<-         (apply(te.run,1,var,na.rm=T))^.5
  
  
  iv.obj<-list("replaceme.mat"=replaceme.run,"encourage.mat"=treat.fits.all.run)
  
  time.keep[[4]]<-proc.time()
  
  output<-list("te.ave"= te.out,"se.eff"=ses.out,"CIs"=CIs,
               "coef.run"=coef.run,"fits.run"=fits.out,"y.se.run"=y.se.run,
               "y.partial.run"=rowMeans(resy.run,na.rm=T),"cvalpha"=cvalpha2,
               "treat.fitted"=rowMeans(treat.fits.all.run,na.rm=T),"iv.obj"=iv.obj,
               "old"=CIs_old,"timing"=time.keep, "explainedT"=mean(rsqt.run),
               "tebases" = list("theta"=baseste.full,"tau"=baseste.del.full,"resy"=resy.full,"rest"=rest.full)
  )
  
  class(output) <- c("sparseregTE","sparsereg")
  
  return(output)
}



########################################
##  Functions for TEsplit ----
########################################

##  Check correlation ----
#' @export
check.cort<-function(x,tspline0=tspline,Xt0=Xte,resy0=resy,sample.mat){
  ts1<-tspline0[,x[1]]*Xt0[,x[2]]*Xt0[,x[3]]
    cors.check<-apply(sample.mat,2,FUN=function(ind,a=resy0,b=ts1)
    #pcor(a[ind],b[ind])
      arma_cor(a[ind],b[ind])
    )
    if(sum(is.na(cors.check)) >0) return(0)
    if(abs(mean(sign(cors.check)))<1) return(0)
    abs(colMedians(as.matrix(cors.check)))[1]
}
##  Construct basis  ----

#' @export
construct.cort<-function(x,tspline0=tspline,Xt0=Xte){
  ts1<-tspline0[,x[1]]*Xt0[,x[2]]*Xt0[,x[3]]
  if(arma_var(tspline0[,1])>0){X.pred<-cbind(1,Xt0[,x[2]],Xt0[,x[3]],Xt0[,x[2]]*Xt0[,x[3]])
  ts2<-fastLm(y=ts1,X=X.pred)$res
  } else ts2<-ts1
  return(as.matrix(ts2))
}

##  Add names to treatment interaction bases ----
#' @export
names.cort<-function(x,tspline0=tspline,Xt0=Xte){
  paste(colnames(tspline0)[x[1]],colnames(Xt0)[x[2]],colnames(Xt0)[x[3]],sep=
          "_x_" )
}


##  Convert a covariate to bases ----
#' @export
make.bs<-function(x,name1){
  x<-as.vector(scale(x))
  if(length(unique(x))<=3 ){
    m2<-cbind(1,x)
  } else{
  m2<-cbind(1,bs.me(x))
  
  }
  m2
}

##  Convert an X matrix to bases ----
#' @export
make.Xte<-function(X){
  b.dum<-make.bs(rnorm(nrow(X)),"dum")
  Xte<-matrix(NA,nrow=nrow(X),ncol=ncol(X)*(ncol(b.dum)+2))
  colnames(Xte)<-paste("Xte",1:ncol(Xte),sep="_")
  for(i.X in 1:ncol(X)) {
    x.b<-make.bs(X[,i.X],colnames(X)[i.X])
    first.na<-min(which(is.na(Xte[1,])))
    Xte[,first.na:(first.na+ncol(x.b)-1)]<-x.b
    colnames(Xte)[first.na:(first.na+ncol(x.b)-1)]<-colnames(x.b)
  }
  Xte<-Xte[,colSums(is.na(Xte))==0]
  Xte<-Xte[,check.cor(Xte,thresh=1e-4)$k]
  #Xte<-Xte[,apply(Xte,2,arma_var)>0]
  #Xte<-apply(Xte,2,scale)
  Xte<-cbind(apply(X,2,FUN=function(x) (x-min(x))/(diff(range(x))) ),Xte)
  Xte<-Xte[,apply(Xte,2,arma_var)>0]
  Xte<-cbind(1,Xte)
  colnames(Xte)[1]<-"(Intercept)"
  return(Xte)
  
 
}

##  Cleaning the covariate matrix ----
#' @export
cleanX<-function(X,  interpretable){
  X0<-X
  X<-apply(X,2,rank)
  X<-apply(X,2,FUN=function(x)(x-min(x))/diff(range(x)) )
  #X<-name_X(X)
  if(length(colnames(X))==0 ) colnames(X) <-paste("X",1:ncol(X),sep="_")
  X<-apply(X,2,scale)
  if(!interpretable){
    X<-cbind(X,svd(X)$u)
    colnames(X)[(ncol(X0)+1):ncol(X)]<-paste("pc",1:(ncol(X)-ncol(X0) ),sep="_")
  }
  X<-X[,check.cor(X,thresh=1e-4)$k]
  
}

##  Not sure what this is for ----
#' @export
interpretTE<-function(obj,varname,exact=FALSE){
  if(exact) {
    cols.use<-which(varname==rownames(obj$coef.run))
  } else {
    cols.use<-grep(varname,rownames(obj$coef.run))
  }
  fits<-obj$treatbases[,cols.use]%*%obj$coef.run[cols.use,]
  fits.mean<-rowMeans(fits)
  fits.ci<-t(apply(fits,2,quantile,c(0.05,0.95)))
  
  output<-list("mean"=fits.mean,"cis"=fits.ci)
  return(output)
  
}

##  Partialing random effects out ----
#' @export
residualize_sis<-function(y,X,id,X.lin=NULL){
  if(length(id)==0){
    y.res<-lm(y~cbind(X,X.lin))$res
  }
  if(length(id)>0){
    id<-as.matrix(id)
    which.add<-(1:ncol(id))#sapply(colnames(id),nchar)==0
    if(sum(which.add)>0) colnames(id)[which.add]<-paste("re",1:sum(which.add),sep="_")
    form.re<-paste("(1 |",colnames(id),")", collapse=" + ")
    form.re<-as.formula(paste("y~",form.re,collapse=" "))
    y.res<-y-fitted(lmer(form.re,data=data.frame(X,id)))
  }
  y.res
}

##  Making bases to interact w the treatment ----
#' @export
makeTEbases<-function(y,X, id,X.lin=NULL,gs0=FALSE,SIS.num0=50){
  X<-apply(X,2,scale)
  y.res<-residualize_sis(y,X,id,X.lin)
  Xy <- spline_tensor(y = y.res, X = X, degrees=c(3), SIS.num=100)$X.spline
  basesy1 <- cbind(X,splineinter(y = y.res, X = Xy, SIS.num=SIS.num0)$X.inter)
  basesy1<-basesy1[,check.cor(basesy1,thresh=1e-4)$k]
  
  for(i.dum in 1:5) colnames(basesy1)<-gsub(" ", "", colnames(basesy1))
  for(i.dum in 1:5) colnames(basesy1)<-gsub(":", "_", colnames(basesy1))
  basesy1<-basesy1[,unique(colnames(basesy1))]
  if(gs0) basesy1<-gramschmidt(y.res,basesy1)
  basesy1<-cbind(X,X.lin,basesy1)
}


##  Make matrix of interactions to be used throughout ----
#' @export
make.inter.schedule<-function(Xte,resy,rest,X,n,tspline,sis.tx0=NULL){
  inter.schedule.0<-cbind(
    rep(1:ncol(Xte),each=ncol(Xte)),
    rep(1:ncol(Xte),times=ncol(Xte))
  )
  inter.schedule<-NULL
  for(i.inter in 1:ncol(tspline)){
    inter.schedule<-rbind(inter.schedule,cbind(i.inter,inter.schedule.0))
  }
  
  inter.schedule<-inter.schedule[inter.schedule[,2]<=inter.schedule[,3],]
  
  resy.temp<-resy#as.vector(lm(resy~ rest*apply(X,2,scale)) )$res
  sample.mat<-matrix(NA,nrow=ceiling(length(resy)/2),ncol=5)
  for(i.s in 1:5) sample.mat[,i.s]<-sample(1:n,nrow(sample.mat),FALSE)
  cors.t<-apply(inter.schedule,1,FUN=check.cort,tspline0=tspline,Xt0=Xte,resy0=resy.temp,sample.mat)
  
  if(length(sis.tx0)==0){
    sis.tx<-floor(20*(1+n^(.2)))
    sis.tx<-min(sis.tx,length(cors.t)) 
  } else{
    sis.tx<-sis.tx0
  }
  ##New!
  #sis.tx<-200
  sis.ind<-sort(cors.t,dec=T,ind=T)$ix[1:sis.tx]
  inter.schedule<-inter.schedule[sis.ind,]
  inter.schedule
}


## Make treatment effect basees ----
#' @export
make.baseste<-function(rest,resy,X,inter.schedule,tspline,tspline.del,Xte,X.lin=NULL){
  X2<-cbind(X.lin,X)
  baseste<-cbind(rest,apply(inter.schedule,1,construct.cort,tspline0=tspline,Xt0=Xte))
  baseste.del<-cbind(1,apply(inter.schedule,1,construct.cort,tspline0=tspline.del,Xt0=Xte))
  colnames(baseste.del)<-colnames(baseste)<-c(
    "treat",#paste("treat",colnames(X2),sep=":"),
    apply(inter.schedule,1,names.cort,tspline0=tspline,Xt0=Xte)   
  )
  
  ste<-NULL
  keeps<-check.cor(baseste,thresh=1e-4)$k
  baseste<-baseste[,keeps]
  baseste.del<-baseste.del[,keeps]
  output<-list("baseste.del"=baseste.del,"baseste"=baseste)
  return(output)
}


## Make replaceme vector ----
make.replaceme<-function(y,id){
  n<-length(y)
  replaceme<-sample(rep(1:2,length(y))[1:length(y)])
  if(length(id)>0){
    
    subsamp<-function(z,x){
      sample(rep(1:2,length(y))[1:sum(x==z)])
    }
    replaceme<-rep(0,n)
    for(i.r in unique(id)) replaceme[id==i.r]<-subsamp(id,i.r)
  }
  replaceme
}

##  Make treatment effect bases w/o checking correlations ----
#' @export
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

#' @export
## For making screened out bases ----
make.basesscreen<-function(y,basesy1,id){
  colnames(basesy1)<-paste(colnames(basesy1),"_",1:ncol(basesy1),sep="")
  sy<-sparsereg(y,basesy1,EM=F,verbose=F, tol = 0.001,id=id,gibbs=50,burnin=50,thin=10, sparseregweights = TRUE, iter.initialize = 20)
  basesy1<-basesy1[,colnames(basesy1)%in%rownames(coef(sy))]
  sy$coefficients<-rowMeans(sy$coefficients)
  output<-list("obj"=sy,"bases"=basesy1)
}
#' @export
## Gram Schmidt orthogonalization
#X<-bases.temp; X.cor<-NULL
gramschmidt<-function(y,X,X.cor=NULL){
  if(is.null(X.cor)) X.cor<-X
  n<-length(y)
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


## For creating density off confidence interval ----
#' @export
make.dens<-function(z,p){
  d1<-density(z,na.rm=T)
  ys.adj<-d1$y/sum(d1$y)
  cumsum(ys.adj)
  o2<-d1$x[which(cumsum(ys.adj)>1-p/2)][1]
  o1<-d1$x[rev(which(cumsum(ys.adj)>p/2))][1]
  c(o1,o2)
}

#' @export
##  Iterated sure screening ----
isis.te<-function(rest,X,resy,quants,lms,id=NULL){
  n<-length(rest)
  ## Cubic spline of treatment ----
  tspline<-cbind(rest,bs.me(rest))
  tspline.del<-cbind(rest*0+1,dbs.me(rest))
  nc<-ncol(tspline)-1
  colnames(tspline.del)<-colnames(tspline)<-c("treatlin",paste("treatspline",1:nc,"_"))
  
  Xte<-make.Xte(X)

  ## Create interaction schedule ----
  inter.schedule0<-make.inter.schedule(Xte,resy,rest,X,n,tspline)
  b1<-make.baseste(rest,resy,X,inter.schedule0,tspline,tspline.del,Xte)
  ste0<-sparsereg(resy,b1$baseste,EM=T,verbose=F, tol = 0.001,id=id, sparseregweights = TRUE, iter.initialize = 20)
  resy2<-resy-ste0$fitted.values
  #betas.test<-rev(rev(ste0$coefficients)[1:nrow(inter.schedule0)])
  #betas.keep0<-(abs(betas.test)>(0.0001*sd(resy)))
  
  #cor(cbind(0,b1$baseste.del)%*%ste0$coefficients,m1$mce.true)
  inter.schedule1<-make.inter.schedule(Xte,resy2,rest,X,n,tspline,sis.tx0=10)
  b2<-make.baseste(rest,resy2,X,inter.schedule1,tspline,tspline.del,Xte)
  
  inter.schedule<-unique(rbind(inter.schedule0,inter.schedule1))
  b1<-make.baseste(rest,resy,X,inter.schedule,tspline,tspline.del,Xte)
  
  baseste<-b1$baseste
  baseste.del<-b1$baseste.del
  
  ##Screen out main effects
  
  means.te<-colMeans(baseste)
  sds.te<-apply(baseste,2,sd)
  for(i.b in 1:ncol(baseste)) {
    baseste[,i.b]<-(baseste[,i.b]-means.te[i.b])/sds.te[i.b]
    baseste.del[,i.b]<-(baseste.del[,i.b]-means.te[i.b])/sds.te[i.b]
  }
  baseste.del.0<-baseste.del
  
  
  output<-list("b1"=b1, "inter.schedule"=inter.schedule,"means.te"=means.te,"sds.te"=sds.te,"Xte"=Xte, "full.inter.schedule"=rbind(inter.schedule0,inter.schedule1),"tspline"=tspline,"tspline.del"=tspline.del)
  return(output)
}




##  Bspline functions, adding a curve shooting up on both sides ----
bSpline2<-function(x,...){
  b1<-bSpline(x,...)
  b2<-bSpline(-x,...)
  cbind(b2[,ncol(b2)],b1)
}

dbs2<-function(x,...){
  b1<-dbs(x,...)
  b2<-dbs(-x,...)
  cbind(-b2[,ncol(b2)],b1)
}

##  Make a bspline matrix from a vector ----
#' @export
bs.me<-function(x){
  #sd.x<-sd(x)
  #x<-x/sd.x
  m2<-cbind(x,bSpline2(x,df=3),bSpline2(x,df=4),bSpline2(x,df=5),bSpline2(x,df=6))
 return(m2)
  
}

##  Derivative of bspline from bs.me ----
#' @export
dbs.me<-function(x){
  sd.x<-1
  #sd.x<-sd(x)
  #x<-x/sd.x
   m1<-cbind(1,dbs2(x,df=3),dbs2(x,df=4),dbs2(x,df=5),dbs2(x,df=6))
  return(m1/sd.x)
  
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

