library(MDEI)
library(tictoc)
n <- 1000

X <- matrix(rnorm(n*5), nrow = n)
treat <- rnorm(n)+X[,2]

theta.true <- treat^2
theta.true <- theta.true - mean(theta.true)
tau.true <- 2*treat
Ey.x.true <- X[,1]^2+X[,1]
Ey.x.true <- Ey.x.true - mean(Ey.x.true)

y <- theta.true + X[,2] + rnorm(n)

set.seed(1)
tic()
m1 <- MDEI(y, treat, X, splits=10, alpha=.9)
toc()

cor(m1$tau.est, tau.true)
cor(m1$theta.est, theta.true)

plot(treat,tau.true,type="l")
points(treat,m1$tau.est, pch=19, cex=.5)


plot(sort(treat),sort(treat)^2-mean(sort(treat)^2), type="l", ylim=range(m1$theta.est))
points(treat,m1$theta.est-mean(m1$theta.est), pch=19, cex=.5)


mean(apply(m1$CIs.tau-2*treat,1,prod)<0)

