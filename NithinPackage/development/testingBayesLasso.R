n<-200
p<-100

X <- matrix(rnorm(n*p),nrow=n)
beta.true <- rep(0,p)
beta.true[1:4] <- c(-1,1,-.5,.5)

X<- apply(X,2,scale)
y <- X%*%beta.true+rnorm(n)

alpha.max <- 2*max(abs(t(X)%*%(y-mean(y))))
alpha.max <- max(alpha.max,p*10)
alpha.seq <- seq(alpha.max,p/10,length=10)


bayesLasso(y,cbind(1,X),1,0.001)$coef
lm(y~X)$coef

g1<-GCV(y,cbind(1,X),alpha.seq,0.001)

bayesLasso(y,X,alpha.seq[which.min(g1)],0.001)$coef
