\name{interp.dataset}
\alias{interp.dataset}

\title{Interpolate a dataset}
\description{
Given a data frame of variables measured along a temporal or spatial gradient, interpolate each variable to new values of the gradient.  Useful for interpolating sediment core data to the depths ot ages of another sequences, or to evenly spaced intervals.
}
\usage{
interp.dataset(y, x, xout, method=c("linear","loess","sspline"), 
     rep.negt=TRUE, span=0.25, df=min(20, nrow(y)*.7), \dots)
}

\arguments{
  \item{y}{ data frame to be interpolated. }
  \item{x}{ numeric vector giving ages, depths (ie. x-values( for data frame to be interpolated. }
  \item{xout}{ numeric vector of values to interpolate to. }
  \item{method}{ interpolation method, should be an unambiguous abbreviation of either linear, loess, sspline
  or aspline. See details.}
  \item{rep.negt}{ logical to indicate whether or not to replace negative values with zero in the interpolated data. }
  \item{span}{ span for loess, default=0.25. }
  \item{df}{ degress of freedome for smoothing spline, default is the lower of 20 or 0.7 * number of samples. }
  \item{\dots}{ additional arguments to \code{loess} and \code{smooth.spline}. }
}

\details{
Function \code{interp.dataset} interpolates the columns of data frame with rows measured at intervals given by \code{x}, to new intervals given by \code{xout}. This function is useful to interpolation one set of sediment core data to the depth or ages of another, or to a regular set of intervals. Interpolation can be done using linear interpolation between data points in the original series (default) using \code{\link[stats]{approx}}, using \code{\link[stats]{loess}} locally weighted regression, or by \code{\link[stats]{smooth.spline}}. The latter two methods will also smooth the data and additional arguments may be passed to these functions to control the amount of smoothing. 

}
\value{
Function \code{interp.datasets} returns a data frame of the input data interpolated to the values given in \code{xout}. Values of \code{xout} outside the range of the original data are replaced by \code{NA}.
} 

\author{ 
Steve Juggins
}

\seealso{ 
\code{\link[stats]{approx}}, \code{\link[stats]{loess}}, and \code{\link[stats]{smooth.spline}} for details of interpolation methods.
}

\examples{
data(RLGH)
spec <- RLGH$spec
depth <- RLGH$depths$Depth

# interpolate new dataset to every 0.5 cm
# using default method (linear)
x.new <- seq(0, 20, by=0.5)
sp.interp <- interp.dataset(y=spec, x=depth, xout=x.new)
\dontrun{
# examine the results and compare to original data 
strat.plot.simple(spec, depth, sp.interp, x.new)
}
}
\keyword{ hplot }
