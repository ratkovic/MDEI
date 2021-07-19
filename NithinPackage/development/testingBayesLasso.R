n<-100
p<-100

X <- matrix(rnorm(n*p),nrow=n)
X <-X %*% cov(X)
X<- apply(X,2,scale)
beta.true <- rep(0,p)
beta.true[1:4] <- c(-1,1,-.5,.5)

X<- apply(X,2,scale)
y <- X%*%beta.true+rnorm(n)

alpha.max <- 2*max(abs(t(X)%*%(y-mean(y))))
alpha.max <- max(alpha.max,p*10)
alpha.seq <- seq(alpha.max,p,length=100)


bayesLasso(y,cbind(1,X),1,0.001)$coef
lm(y~X)$coef

g1<-GCV(y,cbind(1,X),alpha.seq,1e-4)

#microbenchmark(GCV(y,cbind(1,X),alpha.seq[1:10],1e-4),bayesLasso(y,cbind(1,X),alpha.seq[which.min(g1)],1e-8) )

bayesLasso(y,cbind(1,X),alpha.seq[which.min(g1)],0.001)$coef[2:5]
#sparsereg(y,X,EM=T,verbose=F)$coef[1:5]
#microbenchmark(bayesLasso(y,cbind(1,X),alpha.seq[which.min(g1)],1e-8),GCV(y,cbind(1,X),alpha.seq[1:20],1e-8),sparsereg(y,X,EM=T,verbose=F))


#bayesLasso(y,cbind(1,X),0.01,0.001)$coef[2:5]
