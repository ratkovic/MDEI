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

results.run <- NULL

n <- 250; k <- 5; pot.type <- 5

for(n in c(250,500,750,1000,1500,2500)){
for(k in c(5)){
for(pot.type in c(5)){

options(warn=-1)

#############################################
#############################################

  ## Make data
    m1<-make.data(n, k, pot.type)
    
  ## Run methods ----
    ## Run us
    md <- list()
      md[[1]] <- MDEI(m1$obs, m1$treat,m1$X,samplesplit = TRUE, conformal=TRUE, splits=20)
      md[[2]] <- MDEI(m1$obs, m1$treat,m1$X,samplesplit = TRUE, conformal=FALSE, splits=20)
      md[[3]] <- MDEI(m1$obs, m1$treat,m1$X,samplesplit = FALSE, conformal=TRUE, splits=20)
      md[[4]] <- MDEI(m1$obs, m1$treat,m1$X,samplesplit = FALSE, conformal=FALSE, splits=20)

cover.func <- function(z){
  diffs<-apply(z$CIs.tau-m1$mce.true,1,prod)
  mean(diffs < 0)
}

width.func <- function(z){
  mean(abs(apply(z$CIs.tau,1,diff)))
  
}

results.curr <- c(n,unlist(lapply(md,cover.func)),unlist(lapply(md,width.func)))
results.run <- rbind(results.run, results.curr)

}
}
}

colnames(results.run) <- c("n", 
                           "both_cover","split_cover","conf_cover","neither_cover",
                           "both_width","split_width","conf_width","neither_width"
                           )

ran.num<-round(runif(1)*1e10)
save(results.run,file=paste("results_split/results",ran.num,sep="_"))

