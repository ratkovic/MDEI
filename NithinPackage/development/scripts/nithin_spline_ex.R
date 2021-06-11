library(splines2)
x<-(0:100)/100

splineMat <- bSpline(x, df = 5, knots = median(x))

head(splineMat)

plot(x,splineMat[,1],type="n", ylim=range(splineMat) )
for(i in 1:ncol(splineMat)) lines(x,splineMat[,i])

# derivative of b spline
dsplineMat <- dbs(x, df = 5, knots = median(x))


y

y<- x^2+ rnorm(n)


## Steps:
# 1) Transform x to  a spline matrix
# 2) Transform a set of x's to a concatenated spline matrix (add an intercept, a vector of all 1's)
# 3) Find correlation between y and all two-way interactions of spline matrix

