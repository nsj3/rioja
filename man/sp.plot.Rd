\name{sp.plot}
\alias{sp.plot}
\alias{sp.plot.default}
\alias{sp.plot.formula}

\title{ Plots species distributions along an environmental gradient }
\description{
Plots a matrix of species distributions along a single environmental gradient using \code{lattice}. 
}
\usage{
\method{sp.plot}{default}(y, x, n.cut=5, sort.vars=c("original","wa", 
        "alphabetical"), subset=NULL, sp.scale=c("fixed", "free", 
        "exaggerated"), xlab=NULL, ylab=NULL, sp.max=NULL, 
        cex.max=NULL, as.table=TRUE, \dots)

\method{sp.plot}{formula}(formula, data=NULL, n.cut=5, 
        sort.vars=c("original","wa", "alphabetical"), subset=NULL, 
        sp.scale=c("fixed", "free", "exaggerated"), xlab=NULL, 
        ylab=NULL, sp.max=NULL, cex.max=10, as.table=TRUE, \dots)
}

\arguments{
  \item{y}{ data frame or matrix os species abundances. }
  \item{x}{ vector of environmental values, or a matrix or data.frame with one or two columns (see details). }
  \item{formula}{ model formula, where the left hand side gives a data frame or matrix of species' abundances and the right hand side gives one or two environmental variables. }
  \item{data}{ data frame containing the variables on the right-hand side of the formula. }
  \item{n.cut}{ cut.off for species occurrences. Species with fewer than n.cut occurrences will be omitted from plot. }
  \item{sort.vars}{ how to order variables in the plot). }
  \item{subset}{ a logical vector or list of column indices to include in plot. }
  \item{sp.scale}{ how to scale species abundances in the plot (see details). }
  \item{xlab, ylab}{ labels for the x- and y-axes. }
  \item{sp.max}{ maximum value for scaling abundances (see details). }
  \item{cex.max}{ maximum symbol size for plotting abundances (see details). }
  \item{as.table}{ order of plots on page.  Defualt \code{TRUE} strat top-left.  Set to \code{FALSE} to start bommot left. }
  \item{\dots}{ additional arguments to \code{xyplot}. }
}

\details{
Function \code{sp.plot} plots all or a subset of variables in a data frame (\code{y}) (usually species abundances, with a column for each species) against one or two environmental gradients given in \code{x}. If \code{x} contains a single environmental variable, either as a vector or a matrix / data frame with one column, a matrix of scatterplots is produced, with species abundance on the y-axis and the environmental values on the x-axis. If \code{x} is a matrix or data frame containing two columns, a matrix of plots is produced with the first variable on the x-axis and the second on the y-axis, and symbols scaled according to the species' abundances.  

Individual plots may be scaled in one of three ways: \code{"fixed"} uses a common scaling for species abundances for all plots; \code{"exaaggerated"} uses a common scaling, but any species with a maximum abundance of less than a tenth of the overall maximum is exaggerated by 10, and the text "x10" indicated in the plot; \code{free} uses a separate scaling for each species, and displays the maximum value if plotting against two environmental variables.

\code{cex.max} sets the maximum symbol size for the plots. \code{sp.max} sets the maximum species value that will be set to \code{cex.max}.  By default this will be the overall species maximum but this can be overridden to provide a uniform scaling: for example, with percentage data, setting cex.max to 5 and sp.max to 100 will scale points proportionately so 100 percent gets a size of 5.

The order of variables in the plot may be changed using \code{sort.vars} - \code{sort.vars="wa"} sorts variables according to their wa optima along the first environmental variable.

 The number of rows and columns in the plot may be controlled using the \code{layout} argument passed to \code{\link[lattice]{xyplot}}. For example \code{layout=c(10,1)} will produce a plot with 1 row and 10 columns. Additional arguments may be passed to the underlying \code{\link[lattice]{lattice}} routines: for example \code{col}, \code{pch} and \code{fill} to control the symbol size and colour, and \code{between=list(y=0.5, x=0.5)} to add a space between figures.
 
The function uses lattice to produce the plots so if you call it from a loop or using \code{source} you will need to wrap it in a \code{print} or \code{plot} statement.
}

\value{
Function \code{sp.plot} returns an object of class \code{trellis}.
}

\author{ 
Steve Juggins
}

\seealso{ 
\code{\link[lattice]{lattice}}, \code{\link[lattice]{xyplot}}, and \code{\link{inkspot}} for a more compact way to visualise multiple species distributions along a single gradient. 
}

\examples{
\dontrun{
data(IK)
spec <- IK$spec
env <- IK$env

#basic scatter plot vs SumSST
sp.plot(spec, env$SumSST)

# exaggerate less abundant taxa
sp.plot(spec, env$SumSST, sp.scale="exag")

# separate y-scale for each
sp.plot(spec, env$SumSST, sp.scale="free")

# order by SumSST preference
sp.plot(spec, env$SumSST, sp.scale="free", sort.vars="wa")

# add a y-space
sp.plot(spec, env$SumSST, sp.scale="free", sort.vars="wa", 
between=list(y=0.5))

# two environmental variables
sp.plot(spec, env[, c(1, 3)])

# exaggerate less abundant taxa
sp.plot(spec, env[, c(1, 3)], sp.scale="exag")

# separate y-scale for each
sp.plot(spec, env[, c(1, 3)], sp.scale="free")

# same using the formula interface
sp.plot(spec ~ SumSST + Salinity, data=env, sp.scale="free")

# Exaggerated scaling using transparent colours
sp.plot(spec ~ SumSST + Salinity, data=env, sp.scale="exag", 
col="salmon", pch=21, fill=rgb(0, 0, 1, alpha=0.3), sp.max=100)

# Change the layout (3 * 3 matrix) 
# and start plots bottom left (as.table=FALSE
sp.plot(spec ~ SumSST + Salinity, data=env, sp.scale="exag", 
col="salmon", pch=21, fill=rgb(0, 0, 1, alpha=0.3), sp.max=100, 
layout=c(3,3), as.table=FALSE)
}
}
\keyword{ hplot }
