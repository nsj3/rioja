\name{aber}
\alias{aber}
\docType{data}
\title{Abernethy Forest pollen data}
\usage{
  data(aber)
}
\description{
Pollen stratigraphic data from Abernethy Forest, Scotland, spanning approximately 5500 - 12100 BP (from Birks & Mathews 1978). Data are percentages of 36 dryland pollen taxa in 49 samples.
}
\source{
  Birks, HH & Mathews, RW (1978). Studies in the vegetational history of Scotland V. Late Devensian and early Flandrian macrofossil stratigraphy at Abernethy Forest, Invernessshire. \emph{New Phytologist} \bold{80}, 455-84.
}
\examples{
data(aber)
strat.plot(aber, scale.percent=TRUE, y.rev=TRUE)
}
\keyword{datasets}
