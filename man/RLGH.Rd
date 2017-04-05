\name{RLGH}
\alias{RLGH}

\docType{data}
\title{Diatom stratigraphic data from the Round Loch of Glenhead, Galloway, Southwest Scotland}
\usage{
  data(RLGH)
}
\description{
Diatom stratigraphic data from the Round Loch of Glenhead, Galloway, Southwest Scotland from core K05, first published in Allott et al. (1992) and re-analysed in Juggins et al. (1996) and Battarbee et al. (2005).  Data are relative abundances (percentages) of a subset of 41 diatom taxa in 20 samples, and includes all taxa with a maximum abundance of 1 percent in any core sample. Dataset is a list with the following named elements: \code{spec} diatom relative abundances, \code{depths} associated sediment core depths and 210Pb ages. Column names in \code{RLGH$spec} are short, 6-character alphanumeric codes for each diatom taxon. \code{RLGH$names} contains the full names for each taxon, in the correct order).  Note that some rare and low abundance taxa have been removed so the percentages do not sum to 100.
}
\references{
Battarbee, R.W., Monteith, D.T., Juggins, S. Evans, C.D., Jenkins, A. & Simpson, G.L. (2005) Reconstructing pre-acidification pH for an acidified Scottish loch: A comparison of palaeolimnological and modelling approaches. \emph{Environmental Pollution}, \bold{137}, 135-149.

Allott, T.E.H., Harriman, R., & Battarbee, R.W. (1992) Reversibility of acidification at the Round Loch of Glenhead, Galloway, Scotland. \emph{Environmental Pollution}, \bold{77}, 219-225.

Juggins, S., Flower, R., & Battarbee, R. (1996) Palaeolimnological evidence for recent chemical and biological changes in UK Acid Waters Monitoring Network sites. \emph{Freshwater Biology}, \bold{36}, 203-219.
}
\examples{
data(RLGH)
names(RLGH$spec)
names(RLGH$depths)
}
\keyword{datasets}
