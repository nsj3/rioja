\name{inkspot}
\alias{inkspot}

\title{Two-way ordered bubble plot of a species by sites data table}

\description{
Plots a two-way ordered bubble plot of a species by sites data table, where species are rows and sites are columns.  The sites can be ordered and the functions will sort species to cluster abundances on the diagonal.
}

\usage{
inkspot(data, gradient=1:nrow(data), use.rank=FALSE, 
      reorder.species = TRUE, x.axis=c("sites", "gradient", 
      "none"), x.axis.top=FALSE, site.names=NULL, spec.names=NULL, 
      pch=21, cex.max=3, col="black", bg="darkgrey", 
      legend.values=c(2, 5, 10, 20, 50), \dots)
}

\arguments{
  \item{data}{ data frame to be plotted. }
  \item{gradient}{ a vector for ordering sites along the x-axis. }
  \item{use.rank}{ logical to indicate that the rank rather than absolute values of gradient should be used to plot site positions. Defaults to \code{FALSE}. }
  \item{reorder.species}{ should species be reordered to reflect pattern in site ordering?  Defaults to \code{TRUE}. }
  \item{x.axis}{ controls labelling of x.axis. \code{sites} uses site names, \code{gradient} plots an axis relecting values of the supplied gradient, and \code{none} omits labels and draws ticks at the site positions. }
  \item{x.axis.top}{ logical to include  an x-axis on the top of the figure labelled with values of the gradient. }
  \item{site.names, spec.names}{ character vectors of site or species names to annotate the axes.  Defaults to row and column names. }
  \item{cex.max}{ maximum size of plotting symbol.  Symbols are scaled so maximum species abundance has a symbol of this size. }
  \item{pch, col, bg}{ plotting symbol and line / fill colour. }
  \item{\dots}{ additional arguments to \code{plot}. }
  \item{legend.values}{ if not null, places a legend in the top-left corner displaying the listed values. }
}

\details{
 Function \code{inkspot} plots a two-way table of species by sites as a bubble plot, with sites ordered along the x-axis, species on the y-axis, and abundance indiacted by scaled symbols ("bubbles"). It is a useful way to visualise species distribution along an envionmental, spatial or temporal gradient. If \code{gradient} is not given sites are plotting in the order they appear in the input data. Otherwise sites are plotting according to the values in \code{gradient}. If site labels overlap (multiple sites at similar values of \code{gradient}), labels can be suppressed \code{x.axis= "none"}, or replaced with the gradient axis \code{x.axis="gradient"}. A gradient axis can also be added to the top x.axis (\code{x.axis.top=TRUE}.  Symbols are scaled so that the maximu abundance has a symbol size of \code{cex.max}. If sites are spaced unvenly along the gradient, or if many labels overlap, sites may be plotted evenly spaced using \code{use.rank=TRUE}. In this case the function will place top axis labels (if requested) at the appropriate positions along the gradient.
}
\value{
   Function \code{inkspot} returns a list with two named elements:
   \item{spec}{index of the species order. } 
   \item{site}{index of the site order. }
}

\author{ 
Steve Juggins
}

\seealso{ 
\code{\link[vegan]{vegemite}} in package \code{vegan} for a tabular alternative.
}

\examples{
data(SWAP)
mx <- apply(SWAP$spec, 2, max)
spec <- SWAP$spec[, mx > 10]
#basic plot of data with legend
inkspot(spec, cex.axis=0.6)

#order sites by pH
pH <- SWAP$pH
inkspot(spec, pH, cex.axis=0.6)

# add a top axis
inkspot(spec, pH, x.axis.top=TRUE, cex.axis=0.6)

# order by pH but plot sites at regular intervals to avoid label overlap
inkspot(spec, pH, use.rank=TRUE, x.axis.top=TRUE, cex.axis=0.6)

# or add long taxon names
oldmar <- par("mar")
par(mar=c(3,12,2,1))
nms <- SWAP$names[mx > 10, 2]
inkspot(spec, pH, spec.names=as.character(nms), use.rank=TRUE, 
x.axis.top=TRUE, cex.axis=0.6)
par(mar=oldmar)
}

\keyword{ hplot }
