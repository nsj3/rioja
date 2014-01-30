\name{gutils}
\alias{hulls}
\alias{figCnvt}

\title{Graphic utilities. }
\description{
Functions to perform simple graphics or enhance existing plots.
}
\usage{
hulls(x, y, gr, col.gr=NULL)

figCnvt(fig1, fig2)
}
\arguments{
  \item{x, y}{ vectors of x, y coordinates. }
  \item{gr}{ factor to grop observations. }
  \item{col.gr}{ a single colour or a vector of colours of length nG, where nG is the number of groups. }
  \item{fig1, fig2}{ original \code{fig} dimensions (\code{fig1}) and new \code{fig2} dimensions (\code{fig2}). See details.}
}

\details{
 Function \code{hulls} is a wrapper for \code{\link[grDevices]{chull}} to add convex hulls to a scatterplot, optionally specifying a different colour for each hull.
  
 Function \code{figCnvt} projects a set of \code{fig} dimensions \code{fig2} with respect to an original set \code{fig1}.  Useful for laying out plots where the ploting region has already been partitioned using \code{fig}.
}

\value{
  Function \code{figCnvt} returns a vector of 4 values specifying the new new figure dimensions.
}

\author{ 
Steve Juggins
}

\examples{
data(iris)
with(iris, plot(Sepal.Width, Sepal.Length, col=as.integer(Species)))
with(iris, hulls(Sepal.Width, Sepal.Length, gr=(Species)))
}
\keyword{ aplot }
