\name{Ponds}
\alias{Ponds}

\docType{data}
\title{Southeast England ponds and pools diatom and water chemistry dataset. }
\usage{
  data(Ponds)
}
\description{
Diatom and associated water chemistry data for 30 small ponds & pools from SE England collected by, and described in Bennion (1994). Dataset is a list with the following named elements: (\code{spec}) diatom relative abundances for 48 selected common taxa, (\code{env}) lake names, UK GB grid references, lake depth (m) and mean lake-water chemistry. Units are ueq/l except pH, conductivity (uS/cm), alkalinity (meq/l), total phoshporus and chlorophyll-a (ug/l), and nitrate (mg/l). Column names in \code{spec} are short, 6-character alphanumeric codes for each diatom taxon. \code{Ponds$names} contains the full names for each taxon, in the correct order).  
}
\source{
Bennion, H. (1994) A diatom-phosphorus transfer function for shallow, eutrophic ponds in southeast England. \emph{Hydrobiologia}, \bold{275/276}, 391-410.
}
\examples{
data(Ponds)
names(Ponds$spec)
hist(Ponds$env$TP)
}
\keyword{datasets}
