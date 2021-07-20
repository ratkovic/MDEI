n<-1000
p<-50

X <- matrix(rnorm(n*p),nrow=n)
X <-X
X<- apply(X,2,scale)
beta.true <- rep(0,p)
beta.true[1:4] <- c(-1,1,-.5,.5)

y <- X%*%beta.true+rnorm(n)

alpha.max <- 10*max(abs(t(X)%*%(y-mean(y))))
#alpha.max <- max(alpha.max,p*10)
alpha.seq <- seq(alpha.max,p,length=20)


bayesLasso(y,cbind(1,X),10*p,1e-8)$coef[1:10]
lm(y~X)$coef

g1<-GCV(y,cbind(1,X),alpha.seq,1e-4)
g1$beta[1:5]


#microbenchmark(GCV(y,cbind(1,X),alpha.seq[1:10],1e-4),bayesLasso(y,cbind(1,X),2*p,1e-8) )

g.fun <- function(a) bayesLasso(y,cbind(1,X),a,1e-8)$GCV

g2<-sapply(alpha.seq,g.fun)

minalpha <- which.min(g2)
bayesLasso(y,cbind(1,X),alpha.seq[minalpha],1e-8)
#bayesLasso(y,cbind(1,X),12,1e-8)$GCV
#microbenchmark(bayesLasso(y,cbind(1,X),alpha.seq[2],1e-8),GCV(y,cbind(1,X),alpha.seq,1e-8),sparsereg(y,X,EM=T,verbose=F,usesparseregweights=T),times=10)


#bayesLasso(y,cbind(1,X),0.01,0.001)$coef[2:5]


#1) For GCV function, can you diagnose why the singular error is showing?  (should be in
# bayeslasso too)

#2) For GCV, start at highest value of alpha and go down until GCV has increased
# 3 times, then return coefficients, Etausqinv etc.
