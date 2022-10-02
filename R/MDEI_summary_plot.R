#' Plot an MDEI Object
#'
#' Plots an object of class MDEI
#' @param obj An object of class MDEI.
#' @param nvars Label for y-axis of figure.  Default is \code{""}.
#' @param ... Additional values to be passed to \code{summary}.
#' @export
#'
#'@rdname plot
plot.MDEI <-
  function(obj,
           xvar = "treat",
           sigval = 0,
           target = "tau",
           colors = c(gray(.7), gray(0)),
           cex.point = 0.5,
           xlabel = "",
           ylabel = "",
           ...) {
    # Set up point estimates and CI
    pointest <- obj$tau.est
    CI <- obj$CIs.tau
    if (target == "theta") {
      pointest <- obj$theta.est
      CI <- obj$CIs.theta
    }
    
    # set up variable to be plotted
    if (xvar[1] == "treat") {
      xvar <- obj$internal$treat
    } else{
      if (is.character(xvar[1]))
        xvar <- obj$internal$X0[, xvar]
    }
    
    cover.curr <- apply(CI - sigval, 1, prod) < 0
    plot(
      xvar,
      pointest,
      type = "n",
      ylim = range(CI),
      xlab = xlabel,
      ylab = ylabel
    )
    lines(xvar[sort(xvar, index.return = T)$ix], (xvar * 0 + sigval)[sort(xvar, index.return =
                                                                            T)$ix])
    segments(
      x0 = xvar,
      x1 = xvar,
      y0 = CI[, 1],
      y1 = CI[, 2],
      col = ifelse(cover.curr,
                   colors[1],
                   colors[2])
    )
    
    points(xvar, pointest, pch = 19, cex = cex.point)
    print("Proporaiton of the time that the confidence interval contains sigval:")
    print(mean(cover.curr))
    return(cover.curr)
  }

#' Summary of an MDEI Object
#'
#' Summary of an object of class MDEI.  
#' @param obj An object of class MDEI.
#' @param features Number of spline bases to include.
#' @export
#' 
#' @return \describe{
#' \item{coeftable}{A table with three columns: the names of selected spline interations, the average coefficient, and
#' proportion of time it was included in the model. Averages over taken over subsamples in the
#' split sample strategy.Note that the coefficients are interactions between
#' spline interactions that can be accessed through obj$internal$Xmat.spline.}
#'}
#'
#'
#'
#'@rdname summary
summary.MDEI <-
  function(obj,
           features = 10){
    coeftable<- obj$coefficients[1:features,]
    return(coeftable)
  }
