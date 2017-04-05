interp.dataset <- function(y, x, xout, method=c("linear","loess","sspline"), rep.negt=TRUE, span=0.25, df=min(20, nrow(y)*.7), ...) {
# Smooth and interpolate a core (x) on depths (y) to new y-intervals (yout)
   method <- match.arg(method)
   lin.f <- function(y, x1, xout) {
      approx(x1, y, xout)$y
   }
   lo.f <- function(y, x1, xout, span, ...) {
       fit <- loess(y ~ x1, span=span, ...)
       predict(fit, newdata=data.frame(x1=xout))
   }
   ss.f <- function(y, x1, xout, df, ...) {
       fit <- smooth.spline(y ~ x1, df=df, ...)
       predict(fit, x=data.frame(x1=xout))$y[, 1]
   }
   if (is.null(method))
      stop("Interpolation method not recognised")
   
   if (method=="linear") {
      res <- apply(y, 2, lin.f, x1=x, xout=xout, ...)
   } else if (method=="loess") {
      res <- apply(y, 2, lo.f, x1=x, xout=xout, span=span, ...)
   } else {
      res <- apply(y, 2, ss.f, x1=x, xout=xout, df=df, ...)
   } 
   if (rep.negt) {
      res[res<0] <- 0
   }
   res
}
