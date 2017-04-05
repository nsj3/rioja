\name{read.Tilia}
\alias{read.Tilia}

\title{Read data in Tilia format}
\description{
Reads data in old (*.til) Tilia format used by the Tilia program for plotting stratigraphic diagrams.
}
\usage{
read.Tilia(fName, removeChron=FALSE)
}

\arguments{
  \item{fName}{ filename to read. }
  \item{removeChron}{ Remove chron variables from data, defaults to FALSE (chron vars will be include in sample list). }
}

\details{
 \code{read.Tila} reads data in the Tilia format. Tilia is a program written by Eric Grimm to plot stratigraphic diagrams and is popular with palynologists. Note that any spaces or other illegal characters in the species codes will be converted to periods on import.
}
\value{
Returns a list with three names elements:
\item{data}{ data frame with species as columns and sites as rows.  Column and row names are taken from the Tilia file. }
\item{vars}{ names and types of each variable in the dataset. }
\item{levels}{ names depths, and optionally ages of each level in the core. }
}

\note{
The code uses a combination of C and C++ and needs 1-byte structure alignment to read the Tilia binary file. It may crash some 64 bit versions of R.

See function \code{\link[neotoma]{read.tilia}} in the package \code{neotoma} for a version to read Tilia *.tlx files.
}

\author{ 
Steve Juggins
}

\seealso{ 
\code{\link[rioja]{read.CEP}}.
}

\examples{
\dontrun{
pth <- system.file("example.datasets/WOLSFELD.TIL", package="rioja")
WOLS <- read.Tilia(pth)
sel <- WOLS$vars$Sums == "A" | WOLS$vars$Sums == "B"
spec <- WOLS$data[, sel]
totals <- apply(spec, 1, sum)
spec.pc <- spec / totals * 100
mx <- apply(spec.pc, 2, max)
spec.sub <- spec.pc[, mx > 5]
age <- WOLS$levels$Chron2
strat.plot(spec.sub, scale.percent=TRUE, yvar=age, y.rev=TRUE,
      wa.order="bottomleft")
}
}
\keyword{ IO }
\keyword{ file }
