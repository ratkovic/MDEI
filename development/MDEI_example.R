#library(MDEI)

devtools::load_all("~/Dropbox/Github/MDEI")
#devtools::install_github('ratkovic/MDEI/tree/development', force=TRUE)
library(tictoc)
# library(MDEI)
n <- 1000

X <- matrix(rnorm(n*5), nrow = n)
treat <- rnorm(n)

theta.true <- treat*2
theta.true <- theta.true - mean(theta.true)
tau.true <- 2
Ey.x.true <- 0#(X[,1]^2+X[,1])
Ey.x.true <- Ey.x.true - mean(Ey.x.true)

y <- theta.true  + rnorm(n,sd=1)*sd(theta.true )

# set.seed(100)
tic()
m1 <- MDEI(y, treat, X, splits=100, alpha=.9)
toc()



cor(m1$tau.est, tau.true)
cor(m1$theta.est, theta.true)

plot(treat,tau.true,type="l")
points(treat,m1$tau.est, pch=19, cex=.5)


plot(sort(treat),sort(treat)^2-mean(sort(treat)^2), type="l", ylim=range(m1$theta.est))
points(treat,m1$theta.est-mean(m1$theta.est), pch=19, cex=.5)



mean(apply(m1$CIs.tau-tau.true,1,prod)<0)
mean(m1$tau.est)


options(device="quartz")
plot(treat,m1$tau.est,type="n",ylim=range(m1$CIs.tau))
lines(treat,tau.true)
segments(x0=treat,x1=treat,y0=m1$CIs.tau[,1],y1=m1$CIs.tau[,2],col=gray(.7))
points(treat,m1$tau.est,pch=19,cex=.5)
