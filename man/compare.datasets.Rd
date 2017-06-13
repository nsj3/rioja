\name{compare.datasets}
\alias{compare.datasets}
\alias{plot.compare.datasets}

\title{Compare datasets for matching variables (species)}
\description{
Compare two datasets and summarise species occurrance and abundance of species recorded in dataset one across dataset two. Useful for examining the conformity between sediment core and training set species data.
}
\usage{
compare.datasets(y1, y2, n.cut=c(5, 10, 20, 50), 
      max.cut=c(2, 5, 10, 20, 50))

\method{plot}{compare.datasets}(x, y, subset=1:nrow(x$obs), \dots) 
}
\arguments{
 \item{y1, y2}{ two data frames or matrices, usually of biological species abundance data, to compare. }
 \item{x}{ an object of class \code{compare.datasets} produce by function 
           \code{compare.datasets}. }
 \item{y}{ original dataset (ie. y1 above) used in comparison. }
 \item{n.cut}{ vector of abundances to be used for species occurrence calculations (see details). }
 \item{max.cut}{ vector of occurences to be used for species maximum abundance calculations (see details). }
 \item{subset}{ a vector giving row indices to plot. Used to limit the number of plots with larges datsets. }
 \item{\dots}{ additional arguments to \code{\link[lattice]{xyplot}}. }
} 

\details{
 Function \code{compare.datasets} compares two datasets. It summarise the species profile (number of occurences etc.) and sample profile (number of species in each sample etc.) of dataset 1.  For those species recorded in dataset 1 it also provides summaries of their occurence and abundance in dataset 2.  It is useful diagnostic for checking the conformity between core and training set data, specifically for identifying core taxa absent from the training set, and core samples with portions of their assemblage missing from the training set.
  
 \code{\link{plot.compare.datasets}} provides a simple visualisation of the comparisons.  It produces a matrix of plots, one for each sample in dataset 1, showing the abundance of each taxon in dataset 1 (x-axis) against the N2 value of that taxon in dataset 2 (y-axis, with symbols scaled according to abundance in dataset 2.  The plots shouls aid identification of samples with high abundance of taxa that are rare (low N2) or have low abundance in the training set. Taxa thar are absent from the training set are indicated with a red "+".
}

\value{
 Function \code{compare.datasets} returns a list with two names elements:
   \item{vars}{ data frame listing for each variable in the first dataset: N.occur = number of occurences in dataset 1, N2, Hill's N2 for species in dataset 1, Max = maximum value in dataset 1, N.2 = number of occurences in dataset 2, N2.2 = Hill's N2 for species in dataset 2, Max.2 = maximum value in dataset 2, N.005, number of occurences where the species is greater than 5 etc. }
   \item{objs}{ data frame listing for each observation in the first dataset: N.taxa = number of species greater than zero abundance, N2, Hill's N2 for samples, Max = maximum value, total =  sample total, M.002 = number of taxa with a maximum abundance greater than 2 2 etc., N2.005 = number of taxa in dataset 1 with more than 5 occurences in 2 dataset 2 etc., Sum.N2.005 = sample total including only those taxa with at least 5 occurrences in dataset 2 etc., M2.005 = number of taxa in dataset 1 with maximum abundance greater than 2 in dataset 2 etc., and Sum. M2.005 = sample total including only those taxa with a maximum abundance greater than 2 in dataset 2 etc. }
   
 Function \code{plot.compare.datasets} returns an object of class \code{trellis} which may be plotted.
}

\author{ 
Steve Juggins
}

\examples{
# compare diatom data from core from Round Loch of Glenhead
# with SWAP surface sample dataset
data(RLGH)
data(SWAP)
result <- compare.datasets(RLGH$spec, SWAP$spec)
result

\dontrun{
#visualise the comparison
plot(result, RLGH$spec)
}
}
\keyword{ utilities }
