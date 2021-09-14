
#library(MDEI)

devtools::load_all("~/Dropbox/Github/MDEI")
#devtools::install_github('ratkovic/MDEI/tree/development', force=TRUE)
library(tictoc)
# library(MDEI)

library(grf)
library(KRLS2)

plot.tsq <- function(pointest, CI){
  cover.curr <- apply(CI-tau.true,1,prod)<0
  plot(treat,pointest,type="n",ylim=range(m1$CIs.tau), xlab="", ylab="")
  lines(treat,treat*0+tau.true)
  segments(x0=treat,x1=treat,y0=CI[,1],y1=CI[,2],col=ifelse(cover.curr,
                                                            gray(.7),
                                                            "black")
  )
  
  points(treat,pointest,pch=19,cex=.5)
  print(mean(cover.curr))
  
}

n <- 1000

set.seed(1); X <- matrix(rnorm(n*5), nrow = n)
set.seed(1); treat <- sort(rnorm(n))

theta.true <- treat^2
theta.true <- theta.true - mean(theta.true)
tau.true <- 2 * treat
Ey.x.true <- 0#(X[,1]^2+X[,1])
Ey.x.true <- Ey.x.true - mean(Ey.x.true)

y <- theta.true  + rnorm(n,sd=1)*sd(theta.true )

# set.seed(100)
tic()
set.seed(1); m1 <- MDEI(y, treat, X, splits=100, alpha=.9)
toc()



cor(m1$tau.est, tau.true)
cor(m1$theta.est, theta.true)

plot(treat,tau.true,type="l")
points(treat,m1$tau.est, pch=19, cex=.5)


plot(sort(treat),sort(treat)^2-mean(sort(treat)^2), type="l", ylim=range(m1$theta.est))
points(treat,m1$theta.est-mean(m1$theta.est), pch=19, cex=.5)



mean(apply(m1$CIs.tau-tau.true,1,prod)<0)
mean(m1$tau.est)
cover.mdei <- apply(m1$CIs.tau-tau.true,1,prod)<0



set.seed(1); grf1 <- causal_forest(X, y, treat, num.trees=10000)
predict.grf  <- predict(grf1, estimate.variance=TRUE )
grf.CIs <- predict.grf$predictions+1.645*cbind(-predict.grf$variance.estimates^.5, predict.grf$variance.estimates^.5)

set.seed(1); krls1 <- krls(y=y,X=cbind(treat,X), epsilon = 0.001)
  sk1<-inference.krls2(krls1)
  k1.est<-(sk1$derivatives[,1])
  k1.se<-sk1$var.derivatives[,1]^.5
  k1.ci<-k1.est+cbind(-1.645*k1.se, 1.645*k1.se)

options(device="quartz")


pdf("tsqex.pdf",h=4,w=16)
par(mfrow=c(1,4))


plot(treat,m1$tau.est,type="n",ylim=range(m1$CIs.tau), xlab="", ylab="")
lines(treat,2*treat)
lines(treat,treat^2-mean(treat^2), lty=2)
points(treat,y,pch=19,cex=.5)
mtext(side=3,"Data Setup")





plot.tsq(m1$tau.est, m1$CIs.tau)
mtext(side=3,"MDEI")

plot.tsq(grf1$predictions, grf.CIs)
mtext(side=3,"GRF")

plot.tsq(k1.est,k1.ci)
mtext(side=3,"KRLS")

dev.off()

