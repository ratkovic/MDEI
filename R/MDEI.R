## New approach.  By split:
## 3: Partial out, learn bases
## 1: Estimation
## 2: Conformal
## Then switch 'er up.

## Primary fitting function; wrapper function for this one is below
fit.singlesubsample <- function(y0, treat0, X0, replaceme0, Xmat0) {
  ## Partial out X's ----
  y <- y0
  treat <- treat0
  X <- X0
  replaceme <- replaceme0
  Xmat <- Xmat0
  
  treat.partial <- partialOut(treat, X, replaceme)
  # treat.pscore <- treat-treat.partial#lm(treat~treat.partial, weights=1*(replaceme==1))$fit
  y.partial <- partialOut(y, cbind( X), replaceme)
  Ey.x <- y - y.partial

  treatmat.theta <- bs.me(treat.partial, "treatment")
  treatmat.tau <- dbs.me(treat.partial)
  
  keeps <- which(as.vector(checkcor(treatmat.theta, .9) == 1))
  treatmat.theta <- treatmat.theta[, keeps]
  treatmat.tau <- treatmat.tau[, keeps]
  
  colnames.treat <- paste("treat", 1:ncol(treatmat.theta), sep = "")
  
  ## Calculate correlations
  # tic("Creating bases")
  bases.obj <-
    createBases(replaceme,
                Xmat,
                y.partial,
                treatmat.theta,
                treatmat.tau,
                ratio = 50)
  # toc()
  
  ##Partial out treatment?
  for(i.p in 2:ncol(bases.obj$MConstruct)){
    lm.p <- lm(bases.obj$MConstruct[,i.p]~treat.partial, weights=1*(replaceme==2))
    bases.obj$MConstruct[,i.p] <- lm.p$res
    bases.obj$MConstructDerivative[,i.p] <- bases.obj$MConstructDerivative[,i.p]-lm.p$coef[2]
    }
  
  n.a <- sum(replaceme == 2)
  p.a <- ncol(bases.obj$Msubsamp)
  alpha.seq <- seq(max(n.a * log(p.a), 10 * p.a), p.a, length = 10)
  
  X.Construct <- cbind(1, bases.obj$MConstruct)[replaceme==2,]
  X.Construct1 <- cbind(1, bases.obj$MConstruct)[replaceme==1,]
  X.Construct2 <- cbind(1, bases.obj$MConstruct)[replaceme==2,]
  
  XpX.Construct <-crossprod(X.Construct1)
  
  #m1 <- mget(ls())
  #save(m1,file="diagnose.Rda")
  g1 <-
    GCV(y.partial[replaceme == 1],
        X.Construct1,
        alphas = alpha.seq,
        tol = 1e-6*sd(y.partial))

  beta.sp <- as.vector(g1$beta)
  errs.loo <- 
    as.vector(y.partial[replaceme == 2]-X.Construct%*%beta.sp)
  
  # tic("Gathering coefficients")
  beta.sparse <- beta.sp[-1][abs(beta.sp[-1]) > 1e-2*sd(y)]
  cormat.sparse <- as.matrix(bases.obj$cormat[,abs(beta.sp[-1]) > 1e-2*sd(y)]+1)
  cormat.sparse[4, ] <- beta.sparse
  
   coefnames.sparse <- apply(as.matrix(cormat.sparse[1:3,]), 2, FUN=function(z){
    c1 <- colnames(treatmat.theta)[z[1]]
    c2 <- colnames(Xmat)[z[2]]
    c3 <- colnames(Xmat)[z[3]]
    
    paste(c1,c2,c3,sep=" x ")
  })
   
   
   colnames(cormat.sparse) <- coefnames.sparse
    # toc()
  
  
  fits.curr <- cbind(1, bases.obj$MConstruct) %*% beta.sp
  te.curr <- cbind(0, bases.obj$MConstructDerivative) %*% beta.sp 
  ## Variance calculations ----
  
  # Variance of fitted value
  
  res.sq <- (y.partial - fits.curr)
  res.sq <- (res.sq - mean(res.sq[replaceme == 1])) ^ 2
  treat2 <- treat
  var1 <-
    ranger(res.sq ~ .,
           data = data.frame(treat2, X),
           case.weights = 1 * (replaceme == 1))
  
  numvar <- 25
  var.treatperm <- 0
  for (i in 1:numvar) {
    treat2 <- sample(treat)
    var.treatperm <-
      var.treatperm + predict(var1, data = data.frame(treat2, X))[[1]] / numvar
  }
  # var.tau <- abs(var.treatperm[replaceme==2]-var1$predictions[replaceme==2])
  var.tau <-
    pmax(var.treatperm[replaceme == 2] - var1$predictions[replaceme == 2], 0)
  output <-
    list(
      "theta.pred" = fits.curr[replaceme == 2],
      "tau.pred" = te.curr[replaceme == 2],
      "var.theta" = var1$predictions[replaceme == 2],
      "var.tau" = var.tau,
      "y.partial" = y.partial[replaceme == 2],
      "Ey.x" = Ey.x[replaceme == 2],
      "errs.loo" = errs.loo,
      "cormat.sparse" = cormat.sparse
    )
  
  
  return(output)
}


#' MDEI function
#'
#' Implements the Method of Direct Estimation and Inference
#' @param y The outcome variable, a vector.
#' @param treat The treatment variable, a vector.
#' @param X A matrix of covariates.
#' @param splits Number of repeated cross-fitting steps to implement.
#' @param alpha The desired level of the confidence band.
#' @export
#'
#' @examples
#' n<-200
#'
#' X <- matrix(rnorm(n*3), nrow = n)
#' treat <- rnorm(n)
#' y <- treat^2 + X[,2] + rnorm(n)
#'
#' set.seed(1)
#' m1 <- MDEI(y, treat, X, splits=5, alpha=.9)
#'
#' # Accuracy
#' cor(m1$tau.est, treat*2)
#' cor(m1$theta.est, treat^2)
#'
#' # Coverage
#' mean(apply(m1$CIs.tau-2*treat,1,prod)<0)


MDEI <- function(y,
                 treat,
                 X,
                 splits = 10,
                 alpha = .9) {
  n <- length(treat)
  X <- apply(X, 2, rank)
  X <- apply(X ,2, FUN=function(x) pnorm(x/(max(x)+1)))
  if(length(colnames(X))!=ncol(X)) colnames(X) <- paste("X",1:ncol(X), sep="_")
  
  Xmat.spline <-
    matrix(NA, nrow = n, ncol = ncol(X) * ncol(bs.me(rnorm(nrow(X)),"rnorm" )) + 10)
  colnames(Xmat.spline) <- paste("init",1:ncol(Xmat.spline),sep="_")
  for (i.X in 1:ncol(X)) {
    col.start <- which(is.na(Xmat.spline[1, ]))[1]
    bmat <- as.matrix(bs.me(X[, i.X], colnames(X)[i.X] ))
    if(ncol(bmat)==1) Xmat.spline[ ,col.start]  <- bmat
    col.stop <- ncol(bmat) + col.start - 1
    Xmat.spline[, col.start:col.stop] <- bmat
    colnames(Xmat.spline)[col.start:col.stop] <- colnames(bmat)
  }
  
  Xmat.spline <- Xmat.spline[, is.finite(Xmat.spline[1, ])]
  Xmat <- Xmat.spline
  keeps <- which(as.vector(checkcor(Xmat, .9)) == 1)
  Xmat <- cbind(1, Xmat[, keeps])
  colnames(Xmat)[1] <- "Intercept"
  
  ## Containers ----
  
  errs.loo.run <- Ey.x.run <-
    y.partial.run <-
    theta.run <-
    tau.run <-
    thetavar.run <- tauvar.run <- matrix(NA, nrow = n, ncol = splits)
  
  coefmat <- NULL
  ## Now, start split sample here ----
  
  for (i.runs in 1:splits) {
    replaceme <- rep(sample(1:3), n/2)[1:n]
    replaceme <- sample(replaceme)
    
    for(i.inner in 1:3){
    singlefit.2 <- fit.singlesubsample(y, treat, X, replaceme, Xmat)
 
    theta.run[replaceme == 2, i.runs] <- singlefit.2$theta.pred
    tau.run[replaceme == 2, i.runs] <- singlefit.2$tau.pred
    thetavar.run[replaceme == 2, i.runs] <- singlefit.2$var.theta
    tauvar.run[replaceme == 2, i.runs] <- singlefit.2$var.tau
    y.partial.run[replaceme == 2, i.runs] <- singlefit.2$y.partial
    Ey.x.run[replaceme == 2, i.runs] <- singlefit.2$Ey.x
    errs.loo.run[replaceme == 2, i.runs] <- singlefit.2$errs.loo
    
    coefmat <- cbind(coefmat,singlefit.2$cormat.sparse)
    replaceme <- (replaceme%%3)+1
    }
    cat("Finished with cross-fit", i.runs, "\n")
  }
  
  se.theta <-
    (apply(thetavar.run, 1, hl.mean) + apply(theta.run, 1, hl.var)) ^ .5
  # ts.theta <- (y.partial.run - theta.run) / se.theta
   ts.theta <- errs.loo.run / se.theta
  
  critical.value.theta <- quantile(abs(ts.theta), alpha)
  
  CIs.theta <-
    apply(y.partial.run, 1, hl.mean) + critical.value.theta * cbind(-se.theta, se.theta)
  
  se.tau <- (apply(tauvar.run, 1, hl.mean) + apply(tau.run, 1, hl.var)) ^
    .5
  critical.value.tau <- critical.value.theta+1#(critical.value.theta ^ 2 + 1) ^ .5
  CIs.tau <-
    rowMeans(tau.run) + critical.value.tau * cbind(-se.tau, se.tau)
  
  coefs.ave <- sort(tapply(coefmat[4,],colnames(coefmat), sum))/(splits*2)
  prop.count <- sort(tapply(coefmat[4,],colnames(coefmat), length))/(splits*2)
  
  coefs.return <- data.frame(names(coefs.ave),coefs.ave,prop.count)
  coefs.return <- coefs.return[sort(prop.count,decreasing = T, index.return = T)$ix, ]
  rownames(coefs.return) <- NULL
  
  allobjs <- mget(ls())
  output <- list(
    "tau.est" = apply(tau.run, 1, hl.mean),
    "CIs.tau" = CIs.tau,
    "theta.est" = apply(theta.run, 1, hl.mean),
    "CIs.theta" = CIs.theta,
    "critical.values" = list("theta" = critical.value.theta, "tau" =
                               critical.value.tau),
    "Ey.x" = rowMeans(Ey.x.run),
    "coefficients"=coefs.return,
    "internal"= allobjs
  )
  return(output)
}
