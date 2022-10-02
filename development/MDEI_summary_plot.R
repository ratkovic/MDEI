devtools::load_all("~/Dropbox/Github/MDEI")
plot.tsq <- function(pointest, CI){
  cover.curr <- apply(CI-tau.true,1,prod)<0
  plot(treat,pointest,type="n",ylim=range(m1$CIs.tau), xlab="", ylab="")
  lines(treat,treat*0+tau.true)
  segments(x0=treat,x1=treat,y0=CI[,1],y1=CI[,2],col=ifelse(cover.curr,
                                                            gray(.7),
                                                            "black")
           )
}

n <- 500

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
