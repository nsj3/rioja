\name{Merge}
\alias{Merge}
\title{Merges two or more data frames on the basis of common column names. }
\description{
Merges two or more data frames on the basis of common column names.
}
\usage{
Merge(..., join="outer", fill=0, split=FALSE, verbose=TRUE)
}

\arguments{
  \item{...}{ two or more data frames to merge. }
  \item{join}{ type of join to perform.  Should be an unambiguous abbreviation of either "outer", "inner", or "leftouter". An outer join produces a data frame that contains all the unique column names of the input data, ie, the union of all input column names. An inner join produces a data frame containing only column names that are common across the input data, ie. the intersection of the input column names. A left outer join produces a data frame containing all column names of the first data frame only: column names that occur in subsequent data frames are omitted. }
  \item{fill}{ value to use to fill non-matched columns. Defaults to zero which is appropriate for species abundance data. }
  \item{split}{ logical to return a single data frame (TRUE) or a named list containing separate (original) data frames with a common set of merged columns (FALSE). Defaults to TRUE (a single data frame). }
  \item{verbose}{ logical to suppress warning messages. }
}

\details{
  \code{Merge} is a utilty function for combining separate datasets of biological count data that have only a subset of taxa (column names) in common.  The outer join is appropriate for merging prior to a joint ordination or for merging a training set and core data prior to environmental reconstruction using the modern analogue technique (MAT).  A left outer join should be used to prepare data for an ordination of a training set and subsequent projection of a second onto the ordination axes. The function is capitalised to distinguish it from \code{merge} in the base R.
}
\value{
If split is set to FALSE the function returns a single data frame with the number of rows equal to the combined rows of the input data and columns sorted alphabetically according to the join type.  Otherwise returns a named list of the merged data frames.
}

\author{ 
Steve Juggins
}

\seealso{ 
\code{\link[base]{merge}}. 
}

\examples{
data(RLGH)
data(SWAP)
# Merge RLGH core data with SWAP training set
# Extract species data from datasets
SWAPsp <- SWAP$spec
RLGHsp <- RLGH$spec
# full outer join for joint ordination of both datasets
comb <- Merge(SWAPsp, RLGHsp)

\dontrun{
# superimpose core trajectory on ordination plot
library(vegan) # decorana
ord <- decorana(comb, iweigh=1)
par(mfrow=c(1,2))
plot(ord, display="sites")
sc <- scores(ord, display="sites")
sc <- sc[(nrow(SWAPsp)+1):nrow(comb), ]
lines(sc, col="red")
title("Joint DCA ordination of surface and core")

# Do the same but this time project core passively
# Note we cannot use data from the outer join since decorana
# will delete taxa only present in the core - the resulting
# ordination model will then not match the taxa in the core
comb2 <- Merge(SWAPsp, RLGHsp, join="leftouter", split=TRUE)
ord2 <- decorana(comb2$SWAPsp, iweigh=1)
sc2 <- predict(ord2, comb2$RLGHsp, type="sites")
plot(ord2, display="sites")
lines(sc2, col="red")
title("DCA with core added \"passively\"")
}
}
\keyword{ utilities }
