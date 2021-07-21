


##
##  Make a bspline matrix from a vector ----
#' @export
bs.me<-function(x){
  m2<-cbind(x,bSpline2(x,df=3),bSpline2(x,df=4))#,bSpline2(x,df=5),bSpline2(x,df=6))
  return(m2)
  
}

##  Derivative of bspline from bs.me ----
#' @export
dbs.me<-function(x){
  m1<-cbind(1,dbs2(x,df=3),dbs2(x,df=4))#,dbs2(x,df=5),dbs2(x,df=6))
  return(m1)
  
}

