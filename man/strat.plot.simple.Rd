\name{strat.plot.simple}
\alias{strat.plot.simple}

\title{ Plot a simple stratigraphic diagram }
\description{
Plots a simple stratigraphic diagram using \code{lattice}. The function is most useful for visualising the effect of interpolating and / or smoothing a dataset.
}
\usage{
strat.plot.simple(y1, x1, y2=NULL, x2=NULL, col=c("blue", "red"), 
     sort.vars=c("original","wa", "alphabetical"), 
     ylim=range(x1), y.rev=FALSE, type=c("b", "l"), 
     subset=c(1:ncol(y1)), \dots)
}
\arguments{
  \item{y1, y2, x1, x2}{ data frames to be plotted and corresponding x-values. If y2 and x2 are specified columns of y2 will be plotted overlain on those of y1. }
  \item{col, type}{ colours and line types for the plots. }
  \item{sort.vars}{ how to order variables in the plot). }
  \item{ylim}{ limits for the y-axis. }
  \item{y.rev}{ logical to reverse y-axis (which is actually the x-variable). }
  \item{subset}{ a vector of column indices or a logical vector giving a subset of columns. }
  \item{\dots}{ additional arguments to \code{xyplot}. }
}

\details{
\code{strat.plot.simple} plots all or a subset of variables in a data frame (y) against a single secondary variable (x). The function is appropriate for plotting stratigraphic data, with the x-variable is plotted on the y-axis and variables in y arrayed across the page. The function is primarily intended to visualise the effect of smoothing by plotting the original data with the interpolated or smoothed version superimposed.  

 The number of rows and columns in the plot may be controlled using the \code{layout} argument passed to \code{link{xyplot}}. For example \code{layout=c(10,1)} will produce a plot with 1 row and 10 columns.
 
 By default, the x-axis has the same scale for all variables.  To allow independent x-scales for each variable include the additional argument: \code{scales=list(x="free"))}.

The function uses lattice to produce the plots so if you call it from a loop or using \code{source} you will need to wrap it in a \code{print} or \code{plot} statement.
}

\value{
Function \code{strat.plot.simple} returns an object of class \code{link{trellis}}.
}

\author{ 
Steve Juggins
}

\seealso{ 
\code{\link{strat.plot}} for better looking stratigraphics diagrams and \code{\link{interp.dataset}} for interpolation and / or smoothing. 
}

\examples{
data(RLGH)
spec <- RLGH$spec
depth <- RLGH$depth$Depth

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
