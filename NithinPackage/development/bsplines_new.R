bs.me <- function(x) {
  x <- rank(x)  ## do this externally, assume x comes in as already rank transformed
  x <- x- length(x)/2
  # n <- length(x)
  # x <- (rank(x) - 1) / (n)
  # x <- 2 * x - 1
  # basis.out <- NULL
  # for (i.basis in 1:degree) {
  #   # basis.out<-cbind(basis.out,cos(i.basis*acos(x)),cos(-i.basis*acos(x)))
  # }
  # basis.out <-
  #   cbind(x, bs(x, degree = 3, knots = median(x)), bs(-x, degree = 3, knots = median(x)))
  basis.out <-
    cbind(x, bSpline(x, degree = 3, knots = 0), bSpline(-x, degree = 3, knots = 0),
          bSpline(x, degree = 3)[,-3], bSpline(-x, degree = 3)[,-3] )
  basis.out <- basis.out[, check.cor(basis.out, 0.000)$k]
}


##
#if(FALSE){
  x<- 1:50
  x<- x- median(x)
  basis.out <- bs.me(x)
  plot(x,basis.out[,2],type="l")
  for(i in 2:ncol(basis.out)) lines(x,basis.out[,i])
#}