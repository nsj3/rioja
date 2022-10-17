\name{riojaPlot}
\alias{riojaPlot}
\alias{addZone2}
\alias{listStyles}

\title{ Plot a stratigraphic diagram }
\description{
Plots a diagram of multiple biological, physical or chemical parameters agains depth or time, as used in geology & palaeoecology.
}
\usage{
riojaPlot(x, y, selVars=NULL, groups=NULL, ...)

addZone2(rp, upper, lower=NULL, ...)

listStyles()

}
\arguments{
  \item{x}{ a data frame or tibble containing variables to plot. }
  \item{y}{ a data frame or tibble containing depth and/or age variables for the y-axis. }
  \item{selVars}{ a character vector of variable names to include in the plot. }
  \item{groups}{ a data frame or tibble with 2 columns giving the variable name (col 1) and factor containing group membership (col 2), used to colour variables by group and / or to generate a summary cumulative plot of group totals. }
  \item{...}{ additional arguments to control the style of the figure (see details below) or colour of zones. }
  \item{rp}{ object returned by riojaPlot. }
  \item{upper, lower}{ upper and lower values of a zone to add. }
}

\details{
\code{riojaPlot} plots a set of variables in a stratigraphic diagram. Diagrams can be plotted as line graphs, silhouettes or bar charts.  {riojaPlot} extends strat.plot in several ways: variables can be grouped and displayed with different colours, a cumulative summary plot can be added automatically based on the grouping, a cluster dendrogram and resulting zones can be added automatically, a secondary y-axis can be added to plot a depth and age axis, and the figures margins are determined automatically depending on the size and length of labels.  

\itemize{
  \item{\code{yvarName}}{ name of the variable to use for primary y-axis scale. If not specified the first column in the y table will be used. }
  \item{\code{secYvarName}}{ name of the variable to use for an optional secondary y-axis scale. }
  \item{\code{showSecAxis}}{ show secondary y-axis.  Default FALSE. }
  \item{\code{yLabel, secyLabel}} { labels for the primary and secondary y-axes.  Uses column names if not specified. }
  \item{\code{scalePC}}{ scale the diagram for percentage data (default TRUE). }
  \item{\code{scaleMinMax}}{ only label min and max values of x-axes to avoid crowding (default FALSE). }
  \item{\code{yrev}}{ reverse the y-axis (default TRUE). }
  \item{\code{yMin, yMax, secYMin, secYMax}}{ minimum and maximum of primary and secondary y-axes. }
  \item{\code{yInterval, secYInterval}}{ intervals for primary and secondary y-axes. }
  \item{\code{autoOrder}}{ rearrange the curves to show trend from bottom left to top right (default FALSE). }
  \item{\code{showBars, showLines, showPoly, showSymbols}}{ show bars, lines, polygons (silhouettes) and / or symbols (defaults to polygons with thin horizintal lines).}
  \item{\code{lwdBar}}{ thickness of bars (1=line). }
  \item{\code{colBar, outlineCol, fillCol, symbCol, zoneCol}}{ colour of bars, lines, polygon outlines and fills, symbols and zones. }
  \item{\code{barTop}}{, plot horizontal bars on top or below shaded silhouettes (default TRUE). }
  \item{\code{xAxisFontSize, yAxisFontSize, yLabelFontSize, nameFontSize, cumulSize}}, { font sizes (size relative to default of 1). }
  \item{\code{tickLen}}{ tick length for axes, in logical units (default -0.3). }
  \item{\code{nameAngle}}{ angle (in degrees) of names (default=90). }
  \item{\code{showClust}}{ add a CONISS cluster dendrogram. }
  \item{\code{dataTransformation}}{ transform data before cluster analysis ("none", "sqrt", or "scale"). }
  \item{\code{ClustUseSelected}}{ use selected variables (ie. those listed in selVars) or all variables when calculating cluster dendrogram (default FALSE). }
  \item{\code{showExag}}{ show exaggereted curves (default FALSE. }
  \item{\code{exagCol}}{ colour for exaggerations, either a colour or "auto". }
  \item{\code{exagMult}}{ multiplication factor for exaggerations (default 2). }
  \item{\code{showGroups}}{ show curves with different colours (default FALSE). }
  \item{\code{showCumul}}{ add a cumulative summary based on groups on right of diagram (default FALSE). }
  \item{\code{groupCol1-6}}{ colours for a maximum of 6 groups. }
  \item{\code{nameStyleBreakLong}}{ break long names (default TRUE). }
  \item{\code{nameStylenBreak}}{ number of characters to break after (default 20). }
  \item{\code{nameStyleItalicize}}{ italicize names. }
}

\code{addZone2} adds horizontal zone liens to an existing diagram and \code{listStyles} simply returns a list of modifiable styles and their default values.

These function uses base graphics and \code{fig} to split the screen and may be incompatible with \code{par(mfrow)} and 
\code{split.screen}.
}

\value{
Returns (invisibly) a list containing the following objects:

box	Vector of 4 values giving the coordinates of the left, right, bottom and top of the plotting area, in relative units.

usr	
Ranges of the plotting area, in data units.

yvar	
Variable used for the y-axis.

ylim	
Limits of the y-axis.
}

\author{ 
Steve Juggins
}

\seealso{ 
\code{\link{strat.plot}}. 
}

\examples{

data(aber)

spec <- aber$spec
# replace species codes with full taxon names
colnames(spec) <- aber$names$Name 
yvar <- aber$ages
types <- aber$types
# convert pollen types to a factor
types$Group <- factor(types$Group, levels=c("Trees", "Shrubs", "Herbs"))

# remove taxa < 5% max abundance
mx <- apply(spec, 2, max)
selTaxa <- names(mx[mx > 5])

riojaPlot(spec, yvar, selTaxa, types,
   yvarName="Depth (cm)", 
   secYvarName="Age (14C years BP)",
   showSecAxis=TRUE, 
   secyLabel = "Age", 
   showGroups=TRUE, 
   showCumul=TRUE,
   showExag=FALSE, 
   exagCol="grey90",
   outlineCol="grey50",
   cumulFontSize=0.7)

}

\keyword{ hplot }
