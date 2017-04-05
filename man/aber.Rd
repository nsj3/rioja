\name{aber}
\alias{aber}
\docType{data}
\title{Abernethy Forest pollen data}
\usage{
  data(aber)
}
\description{
Pollen stratigraphic data from Abernethy Forest, Scotland, spanning approximately 5500 - 12100 BP (from Birks & Mathews 1978).  The data are a list with the following named components: \code{spec} Data are percentages of 36 dryland pollen taxa in 49 samples, (\code{ages}) core depths and ages for the 49 stratigraphic levels, and (\code{names}) codes and full names for the 36 taxa.
}
\source{
  Birks, HH & Mathews, RW (1978). Studies in the vegetational history of Scotland V. Late Devensian and early Flandrian macrofossil stratigraphy at Abernethy Forest, Invernessshire. \emph{New Phytologist} \bold{80}, 455-84.
}
\examples{
data(aber)
strat.plot(aber$spec, scale.percent=TRUE, y.rev=TRUE)
}
\keyword{datasets}
