\name{strat.plot}
\alias{strat.plot}
\alias{addZone}
\alias{addClustZone}

\title{ Plot a stratigraphic diagram }
\description{
Plots a diagram of multiple biological, physical or chemical parameters agains depth or time, as used in geology & palaeoecology.
}
\usage{
strat.plot (d, yvar=NULL, scale.percent=FALSE, graph.widths=1, minmax=NULL, 
            scale.minmax=TRUE, xLeft=0.07, xRight=1, yBottom=0.07, 
            yTop=0.8, title="", cex.title=1.8, y.axis=TRUE, x.axis=TRUE,
            min.width=5, ylim=NULL, y.rev=FALSE, y.tks=NULL, y.tks.labels=NULL, 
            ylabel="", cex.ylabel=1, cex.yaxis=0.8, xSpace=0.01, x.pc.inc=10, 
            x.pc.lab=TRUE, x.pc.omit0=TRUE, wa.order="none", plot.line=TRUE, 
            col.line="black", lwd.line=1, col.symb="black", plot.bar=TRUE, 
            lwd.bar=1, col.bar="grey", sep.bar=FALSE, bar.back=FALSE,
            plot.poly=FALSE, col.poly="grey", col.poly.line=NA, lwd.poly=1,
            plot.symb=FALSE, symb.pch=19, symb.cex=1, x.names=NULL, cex.xlabel=1.1, 
            srt.xlabel=90, mgp=NULL, ylabPos=2, cex.axis=.8, clust=NULL, clust.width=0.1,
            orig.fig=NULL, exag=FALSE, exag.mult=5, col.exag="grey90", exag.alpha=0.2, 
            col.bg=NULL, fun1=NULL, fun2=NULL, add=FALSE, omitMissing=TRUE, ...)

addZone (x, upper, lower=NULL, ...)

addClustZone(x, clust, nZone, ...)
}
\arguments{
  \item{d}{ a matrix or data frame of variables to plot. }
  \item{yvar}{ a vector of depths or ages to use for the y-axis (defaults to sample number). }
  \item{scale.percent}{ logical to scale x-axes for (biological) percentage data. }
  \item{graph.widths}{ a vector of relative widths for each curve, used if \code{scale.percent=FALSE}. }
  \item{minmax}{ 2 * nvar matrix of min and max values to scale each curve if \code{scale.percent=FALSE}. }
  \item{scale.minmax}{ logical to show only min and max values on x-axes (to avoid label crowding). }
  \item{xLeft, xRight, yBottom, yTop }{ x, y position of plot on page, in relative units. }
  \item{title}{ main title for plot. }
  \item{x.names}{ character vector of names for each graph, of same length as \code{ncol(d)}. }
  \item{cex.title}{ size of label for title. }
  \item{y.axis}{ logical to control drawing of left-hand y-axis scale. Defaults to TRUE. }
  \item{x.axis}{ logical or logical vector to control drawing of x-axes. Defaults to TRUE. }
  \item{min.width}{ minimum upper value of x-axis when scaled for percent data. }
  \item{ylim}{ numeric vector of 2 values to control limist of y-axis. Defaults to data range. }
  \item{y.rev}{ logical to reverse y-axis. Defaults to FALSE. }
  \item{y.tks}{ numerical vector listing values of y-axis ticks. }
  \item{y.tks.labels}{ character vector listing values of y-axis labels. }  
  \item{ylabel}{ label for y-axis. }
  \item{ylabPos}{ position for y-axis label.}
  \item{cex.ylabel, cex.yaxis}{ text size for y-axis labels and values. }
  \item{xSpace}{ space between graphs, in relative units. }
  \item{x.pc.inc}{ increment for x-axis values when \code{scale.percent} is TRUE. }
  \item{x.pc.lab}{ logical to control drawing of x-axis values when \code{scale.percent} is TRUE. }
  \item{x.pc.omit0}{ logical to omit initial zero x-axis label when \code{scale.percent} is TRUE. }
  \item{wa.order}{ "none", "topleft" or "bottomleft", to sort variables according to the weighted average with y. }
  \item{plot.line, plot.poly, plot.bar, plot.symb}{ logical flags to plot graphs as lines, silhouettes, bars or symbols. }
  \item{col.line, col.poly.line}{ colour of lines and silhouette outlines. Can be a single colour or a vector of colours, one for each graph. }
  \item{col.poly}{ silhouette fill colour. Can be a single colour or a vector of colours, one for each graph. }
  \item{lwd.line, lwd.poly, lwd.bar}{ line widths for line, silhouette or bar graphs. }
  \item{col.bar}{ colour of bars in a bar graph. \code{col.bar} can be a vector to specify colours of individual bars or graphs.}
  \item{col.symb}{ symbol colour.}
  \item{sep.bar}{ If true, colours in \code{col.bar} are applied to individual bars, otherwise individual graphs. }
  \item{bar.back}{ logical to plot bars behind (TRUE) or on top (FALSE: default) of curves. }
  \item{cex.xlabel}{ size of label for variable names. }
  \item{srt.xlabel}{ rotation angle for variable names. }
  \item{symb.pch, symb.cex}{ symbol type / size. }
  \item{exag}{ logical to add exaggerated curves when \code{plot.poly=TRUE}. Can be a single value or a vector to add exaggeration to individual curves.}
  \item{exag.mult}{ multiplier for exaggerated curves. Can be a single value or a vector to control exaggeration to individual curves.}
  \item{col.exag}{ colour for exaggerated curves. Can be a single value, a vector to control colour of individual curves, or \code{"auto"} for transparent version of main curve.}
  \item{exag.alpha}{ alpha channel for transparent exaggerated curves when \code{col.exag="auto"}. }
  \item{mgp}{ value of mgp for x-axes.  See \code{par} for details. }
  \item{cex.axis}{ text size for x-axis labels.  See \code{par} for details. }
  \item{clust}{ an constrained classification object of class \code{chclust} to add to plot. }
  \item{fun1, fun2}{ custom functions to add additional features to curve. Can be a single function applied to all curves or a vector to apply individual functions to individual curves. fun1 draws behind curves, fun2 draws on top of curves. }
  \item{clust.width}{ width of dendrogram to add to right of plot, in relative units. }
  \item{orig.fig}{ fig values to specify area of window in which to place diagram.  See \code{par} for details. Defaults to whole window. }
  \item{add}{ logical to contol drawing of new page. See \code{par} for details. Defaults to FALSE in which a call to strat.plot will start a new diagram. Set to  TRUE to add a diagram to an existing plot. }
  \item{x}{ a stratigraphic diagram object produced by strat.plot. }
  \item{upper, lower}{ upper and (optional) lower limits of a zone to add to an existing stratigraphic diagram. }
  \item{nZone}{ number of zones to draw. }
  \item{omitMissing}{ remove missing values before plotting.  Defaults to TRUE. }
  \item{col.bg}{ background colour for each curve. }
  \item{...}{ further graphical arguments. }
}
\details{
\code{strat.plot} plots a series of variables in a stratigraphic diagram. Diagrams can be plotted as line graphs and / or bar charts. Samples are plotted on the y-axis by sample number by default but may be plotted against sample age or depth by specifying a variable for yvar. Margins of the plotting area can be changed using xLeft, xRight, yBottom and yTop.  A dendrogram produced by \code{chclust} can be added to the right of the diagram.

The function \code{addZone} can be used to add a horizontal line or box to an existing plot, and 
\code{addClustZone} will add a specified number of zones from a dendrogram (see examples).

The function uses fig to split the screen and may be incompatible with \code{par(mfrow)} and 
\code{split.screen}.
}

\value{
   Returns a list containing the following objects:
   \item{box }{ Vector of 4 values giving the coordinates of the left, right, bottom
    and top of the plotting area, in relative units. }
   \item{usr }{ Ranges of the plotting area, in data units. }
   \item{yvar }{ Variable used for the y-axis. }
   \item{ylim }{ Limits of the y-axis. }
}

\author{ 
Steve Juggins
}

\seealso{ 
\code{\link{chclust}}. 
}

\examples{
library(vegan) ## decorana
data(RLGH)
\dontrun{
# create appropriately sized graphics window
windows(width=12, height=7) # quartz() on Mac, X11 on linux
}
# remove less abundant taxa
mx <- apply(RLGH$spec, 2, max)
spec <- RLGH$spec[, mx > 3]
depth <- RLGH$depths$Depth
#basic stratigraphic plot
strat.plot(spec, y.rev=TRUE)
#scale for percentage data
strat.plot(spec, y.rev=TRUE, scale.percent=TRUE)
# plot by sample depth
strat.plot(spec, yvar = depth, y.rev=TRUE, scale.percent=TRUE,
title="Round Loch of Glenhead", ylabel="Depth (cm)")
# add a dendromgram from constrained cluster analysis
diss <- dist(sqrt(RLGH$spec/100)^2)
clust <- chclust(diss, method="coniss")
# broken stick model suggest 3 significant zones
bstick(clust)
x <- strat.plot(spec, yvar = depth, y.rev=TRUE,
scale.percent=TRUE, title="Round Loch of Glenhead", ylabel="Depth (cm)",
clust=clust)
# add zones
addClustZone(x, clust, 3, col="red")
# use fig to contol diagram size and position
x <- strat.plot(spec, xRight = 0.7, yvar = depth, y.rev=TRUE,
scale.percent=TRUE, title="Round Loch of Glenhead", ylabel="Depth (cm)")
# add curves for first two DCA components of diatom data
dca <- decorana(spec, iweigh=1)
sc <- scores(dca, display="sites", choices=1:2)
strat.plot(sc, xLeft = 0.7, yvar = depth, y.rev=TRUE, xRight=0.99,
y.axis=FALSE, clust=clust, clust.width=0.08, add=TRUE)

# Use custom function to add smooth to curve

sm.fun <- function(x, y, i, nm) {
  tmp <- data.frame(x=y, y=x)
  tmp <- na.omit(tmp)
  lo <- lowess(tmp, f=0.3)
  lines(lo$y, lo$x, col="red", lwd=1)
}

x <- strat.plot(spec, yvar = depth, y.rev=TRUE, scale.percent=TRUE, 
title="Round Loch of Glenhead", ylabel="Depth (cm)", fun1=sm.fun)

# Pollen diagram using built-in Abernethy Forest dataset
data(aber)
depth <- aber$ages$Age
spec <- aber$spec

# basic silhouette plot
strat.plot(spec, yvar = depth, y.rev=TRUE, scale.percent=TRUE, ylabel="Depth (cm)", 
plot.poly=TRUE, col.poly="darkgreen", col.poly.line=NA)

# now with horizontal lines at sample positions
strat.plot(spec, yvar = depth, y.rev=TRUE, scale.percent=TRUE, ylabel="Depth (cm)", 
plot.poly=TRUE, col.poly="darkgreen", plot.bar="Full", col.poly.line=NA)

# add exaggerated curves
strat.plot(spec, yvar = depth, y.rev=TRUE, scale.percent=TRUE, ylabel="Depth (cm)", 
plot.poly=TRUE, col.poly="darkgreen", plot.bar="Full", col.poly.line=NA, exag=TRUE)

# use different colours for trees
xx <- 1:ncol(spec)
cc <- ifelse(xx < 8, "darkgreen", "darkred")
strat.plot(spec, yvar = depth, y.rev=TRUE, scale.percent=TRUE, ylabel="Depth (cm)", 
plot.poly=TRUE, col.poly=cc, plot.bar="Full", col.poly.line=NA, exag=TRUE, col.exag="auto")
}

\keyword{ hplot }
