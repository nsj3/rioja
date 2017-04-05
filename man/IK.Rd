\name{IK}
\alias{IK}

\docType{data}
\title{Imbrie and Kipp foraminifera data}
\usage{
  data(IK)
}
\description{
Core-top foraminifera data from the Atlantic and Indian Oceans and core V12.122 from the Carribean published by Imbrie and Kipp (1971).  The data are a list with the following named components: \code{spec} relative abundances (percentages) of 22 foraminifera taxa in 61 core-top samples, (\code{env}) sea surface temperature and salinity measurements for the core-top samples, and (\code{core}) relative abundances of 28 foraminifer taxa in 110 samples from core V12.122.
}
\references{
Imbrie, J. & Kipp, N.G. (1971). A new micropaleontological method for quantitative paleoclimatology: application to a Late Pleistocene Caribbean core. In \emph{The Late Cenozoic Glacial Ages} (ed K.K. Turekian), pp. 77-181. Yale University Press, New Haven.
}
\examples{
data(IK)
names(IK$spec)
pairs(IK$env)
}
\keyword{datasets}
