devtools::load_all('~/Dropbox/Github/warmup/NithinPackage')

n <- 500
p <- 5

X <- matrix(rnorm(n*p), nrow = n)


X2<-NithinPackage:::bsme(X)
