library(grf)
library(KRLS2)
library(MASS)
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



n <- 5000; k <- 5; pot.type <- 5

for(n in c(250,500,750,1000,1500,2500)){
for(k in c(5)){
for(pot.type in c(5)){
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
    r1<-runmethods(m1)
    cor(r1[[1]]$tau.est,m1$mce.true)
    mean(apply(r1[[1]]$CIs.tau-m1$mce.true,1,prod)<0)
    
    
    g1 <- c(1,gatherresults(r1[[1]]))
    names(g1)[1] <- "type"
    
    mean(g1[[1]]$beta)
    
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

