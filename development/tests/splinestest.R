devtools::load_all('~/Dropbox/Github/warmup/NithinPackage')

library(testthat)

n <- 500
p <- 5

X <- matrix(rnorm(n*p), nrow = n)
X2 <- NithinPackage:::bsme(X)

dim(X)
dim(X2)

y<- rnorm(n)

X2 <- NithinPackage::splineBasesAndCorrs(n, p, X, letters[1:p], y, y, "treat", 100)
