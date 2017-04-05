\name{SWAP}
\alias{SWAP}

\docType{data}
\title{SWAP surface sediment diatom data and lake-water pH. }
\usage{
  data(SWAP)
}
\description{
SWAP (Surface Water Acidification Programme) surface sediment diatom data from Birks et al. (1990) and Stevenson et al. (1990). Dataset is a list with the following named elements: (\code{spec}) diatom relative abundances for 277 taxa in 167 surface samples, (\code{pH}) associated lake-water pH. Column names in \code{spec} are short, 6-character alphanumeric codes for each diatom taxon. \code{SWAP$names} contains the full names for each taxon, in the correct order).
}
\source{
  Birks, H.J.B., Line, J.M., Juggins, S., Stevenson, A.C., & ter Braak, C.J.F. (1990) Diatoms and pH reconstruction. \emph{Philosophical Transactions of the Royal Society of London}, \bold{B 327}, 263-278.
  
  Stevenson, A.C., Juggins, S., Birks, H.J.B., Anderson, D.S., Anderson, N.J., Battarbee, R.W., Berge, F., Davis, R.B., Flower, R.J., Haworth, E.Y., Jones, V.J., Kingston, J.C., Kreiser, A.M., Line, J.M., Munro, M.A.R., & Renberg, I. (1991) \emph{The Surface Waters Acidification Project Palaeolimnology Programme: Modern Diatom / Lake-Water Chemistry Data-Set} ENSIS Ltd, London.
}
\examples{
data(SWAP)
names(SWAP$spec)
hist(SWAP$pH)
}
\keyword{datasets}
