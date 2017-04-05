\name{chclust}
\alias{chclust}
\alias{plot.chclust}
\alias{bstick}
\alias{bstick.chclust}
\title{Constrained hierarchical clustering}
\description{
Constrained hierarchical clustering.
}
\usage{
chclust(d, method = "coniss")

\method{plot}{chclust}(x, labels = NULL, hang = 0.1, axes = TRUE,
      xvar=1:(length(x$height)+1), xlim=NULL, ylim=NULL, 
      x.rev = FALSE, y.rev=FALSE, horiz=FALSE, \dots)

bstick(n, \dots)

\method{bstick}{chclust}(n, ng=10, plot=TRUE, \dots)
}
\arguments{
  \item{d}{ a dissimilarity structure as produced, for example, by \code{dist} or \code{vegdist}. }
  \item{method}{ the agglomeration method to be used. This should be (an unambiguous abbreviation of) either "coniss" or "conslink". }
  \item{x, n}{ a constrained cluster object of class \code{chclust} produced by chclust. }
  \item{xvar}{ numeric vector containing x-coordinates for the leaves of the dendrogram (see \emph{details} below). }
  \item{x.rev,y.rev}{ logical flags to reverse the x- or y-axis (and dendrogram labels). Defaults to \code{FALSE}. }
  \item{horiz}{ logical indicating if the dendrogram should be drawn horizontally or not. Note that y-axis still refers to the dendrogram height even after rotating. }
  \item{xlim, ylim}{ optional x- and y-limits of the plot, passed to the underlying plto function. The defaults for these show the full dendrogram. }
  \item{labels, hang, axes}{ further arguments as in \code{hclust}. }
  \item{ng}{ number of groups to display. }
  \item{plot}{ logical to plot a broken stick model. Defaults to \code{TRUE}. }
  \item{\dots}{ further graphical arguments. Use \code{cex} to change the text size of the x-axis labels, and \code{cex.axis} to change size of the y-axis values. }
}

\details{
\code{chclust} performs a constrained hierarchical clustering of a distance matrix, with clusters constrained by sample order.  Returns an object of class \code{chclust} which can be plotted and interrogated.  See Grimm (1987), Gordon & Birks (1972) and Birks & Gordon (1985) for discusssiom of the coniss and conslink algorithms. The resulting dendrogram can be plotted with \code{plot}.  This is an extension of \code{plclust} that allows the dendrogram to be plotted horizontally or vertically (default). \code{plot} also accepts a numeric vector coordinates for x-axis positions of the leaves of the dendrogram.  These could, for example, be the stratigraphic depths of core samples or geographic distances along a line transect.

\code{bstick.chclust} compares the dispersion of a hierarchical classification to that obtained from a broken stick model and displays the results graphically.  See Bennett (1996) for details. \code{bstick} is a generic function and the default method is defined in package \code{vegan}. If this package is loaded the function may be called using \code{bstick}, otherwise use \code{bstick.chclust}.
}

\value{
Function \code{chclust} returns an object of class \code{chclust}, derived from \code{\link[stats]{hclust}}.
}

\author{ 
Steve Juggins
}

\references{
Bennett, K. (1996) Determination of the number of zones in a biostratigraphic sequence. \emph{New Phytologist}, \bold{132}, 155-170.

Birks, H.J.B. & Gordon, A.D. (1985) \emph{Numerical Methods in Quaternary Pollen Analysis} Academic Press, London.

Gordon, A.D. & Birks, H.J.B. (1972) Numerical methods in Quaternary palaeoecology I. Zonation of pollen diagrams. \emph{New Phytologist}, \bold{71}, 961-979.

Grimm, E.C. (1987) CONISS: A FORTRAN 77 program for stratigraphically constrained cluster analysis by the method of incremental sum of squares. \emph{Computers & Geosciences}, \bold{13}, 13-35.
}

\seealso{ 
\code{\link[stats]{hclust}}, \code{\link[stats]{plclust}}, \code{\link[stats]{cutree}}, \code{\link[stats]{dendrogram}},
\code{\link[vegan]{bstick}}.
}

\examples{
data(RLGH)
diss <- dist(sqrt(RLGH$spec/100))
clust <- chclust(diss)
bstick(clust, 10)
# Basic diagram
plot(clust, hang=-1)
# Rotated through 90 degrees
plot(clust, hang=-1, horiz=TRUE)
# Rotated and observations plotted according to sample depth.
plot(clust, xvar=RLGH$depths$Depth, hang=-1, horiz=TRUE, x.rev=TRUE)

# Conslink for comparison
clust <- chclust(diss, method = "conslink")
plot(clust, hang=-1)
}

\keyword{ cluster }
\keyword{ hplot }
