
devtools::load_all("~/Dropbox/Github/MDEI")
#devtools::install_github('ratkovic/MDEI/tree/development', force=TRUE)
library(tictoc)
# library(MDEI)

options(device="quartz")
plot.tsq <- function(pointest, CI){
  cover.curr <- apply(CI-tau.true,1,prod)<0
  plot(treat,pointest,type="n",ylim=range(c(m1$CIs.tau)), xlab="", ylab="")
  lines(treat,treat*0+tau.true)
  segments(x0=treat,x1=treat,y0=CI[,1],y1=CI[,2],col=ifelse(cover.curr,
                                                            gray(.7),
                                                            "black")
           
  )
  mtext("Estimated Effect",1,line=2.2,cex=1.25)
  mtext("True Effect",2,line=2,cex=1.25)
  
  points(treat,pointest,pch=19,cex=.5)
  print(mean(cover.curr))
  
}

n <- 1000

#set.seed(1); 
X <- matrix(rnorm(n*5), nrow = n)
#set.seed(1); 
treat <- sort(rnorm(n))
# set.seed(1); treat <- X[,1]+rnorm(n)
 
theta.true <- treat#^2/2
theta.true <- theta.true - mean(theta.true)
#theta.true <- X[,1]^2-1
tau.true <- treat*0+1
Ey.x.true <- 0#(X[,1]^2+X[,1])
Ey.x.true <- Ey.x.true - mean(Ey.x.true)

y <- theta.true  + rnorm(n,sd=1)#*sd(theta.true)

tic()
set.seed(1); m1 <- MDEI(y, treat, X, splits=10, alpha=.9)
toc()

  
cor(m1$tau.est, tau.true)
cor(m1$theta.est, theta.true)

plot(treat,tau.true,type="l")
points(treat,m1$tau.est, pch=19, cex=.5)


plot(sort(treat),sort(treat)^2-mean(sort(treat)^2), type="l", ylim=range(m1$theta.est))
points(treat,m1$theta.est-mean(m1$theta.est), pch=19, cex=.5)



mean(apply(m1$CIs.tau-tau.true,1,prod)<0)
table(sign(apply(m1$CIs.tau-tau.true,1,prod)))

mean(m1$tau.est)
cover.mdei <- apply(m1$CIs.tau-tau.true,1,prod)<0

plot.tsq(m1$tau.est,m1$CIs.tau)
