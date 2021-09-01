# library(sparsereg)
library(grf)
library(KRLS2)
# library(randomForest)
# library(splines2)
library(MASS)
# library(splines)
# library(MDEI)
devtools::load_all("~/Dropbox/Github/MDEI")
library(tictoc)

###Functions for the simulations
###
# f1<-c('/Users/ratkovic/Dropbox/Subgroup/MDE_round2/code/sparseregTE_split.R','sparseregTE_split.R')
# sapply(f1,FUN=function(x) ifelse(file.exists(x),source(x),NA))

make.data<-function(n, k, pot.type){
  ##Generate data ----
  var.mat<-diag(k)
  var.mat[var.mat==0]<-.5
  
  X<-mvrnorm(n,rep(0,k),Sig=var.mat)

  sign1<-sign(X[,1])
  g.lin<-(X[,2]^2-1)/4
  
  treat<-g.lin+rnorm(n,sd=1)
  
  if(pot.type==5)  treat<-g.lin*sign1+rnorm(n,sd=1)

  if(pot.type==1) {te.eff<-treat; mce.true=treat*0+1;rmy<-X[,1]+(X[,2]-1)/4 }
  if(pot.type==2) {te.eff<-treat; mce.true= treat*0+1;rmy<-X[,1]+(X[,2]-1)^2/4 }
  if(pot.type==3) {te.eff<-4*sin(treat); mce.true = 4*cos(treat); rmy<-X[,1]+(X[,2]-1)^2/4}
  if(pot.type==4) {te.eff<-4*sin(treat)*X[,1]; mce.true =4*cos(treat)*X[,1];rmy<-(X[,2]-1)^2/4}
  if(pot.type==5) {te.eff<-4*sin(treat)*sign1; mce.true =4*cos(treat)*sign1;rmy<-(X[,2]-1)^2/4}
  

  
  fits.true<-te.eff+rmy
  errs<-rnorm(n,0,1)*sd(fits.true)
  y<- fits.true+errs
  
  obs<-y
  
  output<-list("obs"=obs, "treat"=treat,"X"=X, "mce.true"=as.vector(mce.true), "fits.true"=as.vector(fits.true))
}

runmethods<-function(obj){
  #f1<-c('~/Dropbox/Github/sparsereg/R/sparseregTE_split.R','sparseregTE_split.R')
  #sapply(f1,FUN=function(x) ifelse(file.exists(x),source(x),NA))
  
  y<-obj$obs
  X<-obj$X
  treat<-obj$treat
  mce.true<-obj$mce.true
  
  ## Run sparsereg
  #sapply(f1,FUN=function(x) ifelse(file.exists(x),source(x),NA))
  s1<- MDEI(y, treat, X, splits=20)#sparseregTE(y=y,treat=treat,X=X,splits=40)
  s2<-s1
  
  output<-list(s1,s2)
  return(output)
}

gatherresults<-function(fit,obj=m1){
  bias<-mean(fit$tau.est-mean(obj$mce.true),na.rm=T)
  mae<-mean(abs(fit$tau.est-obj$mce.true))
  cover<-mean(apply((fit$CIs.tau-obj$mce.true),1,prod)<=0,na.rm=T)
  meanci<-mean(abs(fit$CIs.tau[,2]-fit$CIs.tau[,1]),na.rm=T)
  CI.ave<-colMeans(fit$CIs.tau)
  cover.ape<-1*(prod(CI.ave-mean(obj$mce.true))<0)
  #width.ci<-abs(diff(quantile(fit$te.post,c(0.05,.9))))
  #cover.ape<-1*(prod(ci.ave-mean(obj$mce.true))<0)
  width.ci<-diff(CI.ave)
  time.diff<-fit$timing
  output<-c(bias,mae,cover,cover.ape,meanci,width.ci,time.diff)
  names(output)<-c("bias","mae","cover.bci","cover.ape","mean.bci","width.ci","timing")
  return(output)
}

select.run<-k.run<-output.run<-rf.run<-NULL

n <- 250; k <- 5; pot.type <- 1

for(n in c(250,500,750,1000,1500,2500)){
for(k in c(5)){
for(pot.type in c(1,2,3,4,5)){
#for(n in c(1000)){
#for(k in 10){
#for(pot.type in rep(5,1)){
options(warn=-1)

#############################################
#############################################

  ## Make data
    m1<-make.data(n, k, pot.type)
    
  ## Run methods ----
    ## Run us
    tic()
    r1<-runmethods(m1)
    t2<-toc()
    r1[[1]]$timing <- t2[[2]]-t2[[1]]
    
    mean(r1[[1]]$tau.est)
    
    ## Run grf
    pt.rf0<-proc.time()
    c1<-causal_forest(X=m1$X,Y=m1$obs,W=m1$treat,num.trees=4000)
    c2<-predict(c1,estimate.variance=T)
    tau.hat = predict(c1, m1$X, estimate.variance = TRUE)
    sigma.hat = sqrt(tau.hat$variance.estimates)
    calpha<-1.645
    cf.ci<-cbind(tau.hat$predictions-calpha*sigma.hat,tau.hat$pred+calpha*sigma.hat)
    rf.cover<-mean(apply(cf.ci-m1$mce.true,1,prod)<0)
    rf.mae<-mean(abs(tau.hat$predictions-m1$mce.true))
    rf.bias<-mean(average_treatment_effect(c1)[1]-m1$mce.true)
    rf.ci<-mean(abs(calpha*2*sigma.hat))
    ape1<-average_treatment_effect(c1)
    cover.rf<-(mean(m1$mce.true)-ape1[1])/ape1[2] < calpha
    width.apeci.rf<-2*calpha*ape1[2]
    pt.rf1<-proc.time()
    rf.time<-pt.rf1[3]-pt.rf0[3]
    rf.curr<-c(n,k,pot.type,rf.bias,rf.mae,rf.cover,cover.rf,rf.ci,width.apeci.rf,rf.time)
    rf.run<-rbind(rf.run,rf.curr)
    colnames(rf.run)<-c("n","k","pot.type","bias","mae","cover.bci","cover.ape","mean.bci","width.ci","timing")
    
    
    ## Run krls
    pt.k0<-proc.time()
    k1<-with(m1,krls(y=obs,X=cbind(treat,X), epsilon = 0.001))
    sk1<-inference.krls2(k1)
    k1.est<-(sk1$derivatives[,1])
    k1.se<-sk1$var.derivatives[,1]^.5
    k1.ci<-k1.est+cbind(-1.645*k1.se, 1.645*k1.se)
    k1.bias<-sk1$avgderivatives[1]-mean(m1$mce.true)
    k1.cover<-mean(apply(k1.ci-m1$mce.true,1,prod)<0)
    k1.mae<-mean(abs(k1.est-m1$mce.true))
    k1.cover.ape<-1*( (abs(sk1$avgderivatives[1]-mean(m1$mce.true))/sk1$var.avgderivatives[1]^.5)<1.645)
    k1.ci.width<-mean(abs(calpha*2*k1.se))
    k1.ape.ci<-2*calpha*sk1$var.avgderivatives[1]^.5
    pt.k1<-proc.time()
    k.time<-pt.k1[3]-pt.k0[3]
    
    k1.curr<-c(n,k,pot.type,k1.bias,k1.mae,k1.cover,k1.cover.ape,k1.ci.width,k1.ape.ci,k.time)
    
    k.run<-rbind(k.run,k1.curr)
    colnames(k.run)<-c("n","k","pot.type","bias","mae","cover.bci","cover.ape","mean.bci","width.ci","timing")
    
    
    #       bias  mae   cover.bci   cover.ape    mean.bci    width.ci

  ## Gather results
    #o1<-gatherresults(r1[[1]])
    #g1<-t(matrix(unlist(o1),nc=2))
    #colnames(g1)<-names(o1[[1]])
    #g1<-t(as.matrix(cbind(c(1,2),g1)[1,]))
    #colnames(g1)[1]<-"type"
    g1 <- c(1,gatherresults(r1[[1]]))
    names(g1)[1] <- "type"
    
    results.curr<-c(n,k,pot.type,g1)
    names(results.curr)[1:3] <- c("n","k","pot.type")
    
    output.run<-rbind(output.run,results.curr)
    # strings.all<-c("treat","treat:X_1","treat_x_X_1_x_X_2","treat_x_X_1_x_X_2" )
    # string.search<-strings.all[pot.type]
    # s1<-r1[[1]]
    # names.sort<-names(sort(rowMeans(s1$coef.run)))
    # names.select<-c(names.sort[1:2],rev(names.sort)[1:2])
    # select.right<-string.search%in%names.select
    # if(pot.type>=3){
    #   select.both<-intersect(grep("X_1_",names.select),grep("X_2_",names.select))
    #   select.right<-select.right|(length(select.both)>0)
    # }
    # select.run<-rbind(select.run,c(n,k,pot.type,select.right*1))
for(i in 1:20) print("#######################")
print(c(n,k,pot.type))
print("#######################")

}
}
}

ran.num<-round(runif(1)*1e10)
save(output.run,file=paste("results/results",ran.num,sep="_"))
save(rf.run,file=paste("results/rf",ran.num,sep="_"))
save(k.run,file=paste("results/k",ran.num,sep="_"))
# save(select.run,file=paste("results/select",ran.num,sep="_"))

