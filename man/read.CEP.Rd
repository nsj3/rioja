\name{read.CEP}
\alias{read.CEP}
\alias{write.CEP}

\title{Reads and writes data in Cornell Ecology Program (CEP) format}
\description{
Reads and writes data in Cornell Ecology Program (CEP) format used by CANOCO and other programs.
}
\usage{
read.CEP(fName, mValue=-99.9, impZero=0)

write.CEP(x, fName, fTitle=fName, type = "condensed",
    nSig=5, nCouplets=4, mValue=-99.9)
}

\arguments{
  \item{fName}{ filename to read. }
  \item{mValue}{ missing value code in the CEP file, converted to NA on import. }
  \item{impZero}{ value to assign to "missing" taxa in a condensed CEP file. These are usually implicitly zero and the default is to read them as such. }
  \item{x}{ data frame to save in CEP format. }
  \item{fTitle}{ title for header in CEP file. }
  \item{type}{ format of CEP file. Should be an unambiguous abbreviation of either \code{condensed} or \code{full}. \code{condensed} is usually used for species data and \code{full} for environmental data}
  \item{nSig}{ number of significant digits to include in output. }
  \item{nCouplets}{ number of couplets per line to include in condensed file. }
}

\details{
\code{read.CEP} and \code{write.CEP} are functions to read and write data in the Cornell Ecology Program (CEP) format. CEP \emph{condensed} format was originally designed by Mark Hill as an efficient way of storing species abundance data for use with his DECORANA and TWINSPAN programs (Hill 1979a, b), and subsequently adopted by CANOCO (ter Braak & Smilauer, 1998) where \emph{full} and \emph{free} format variations were added. \code{read.CEP} and \code{write.CEP} can read/write data in condensed and full formats. The CEP in \code{read.CEP} is capitalised to distinguish it from a similar function (\code{read.cep)} in the \code{vegan} package, which can also import data in CEP \emph{free} 
  
\code{write.CEP} takes the site and species names from the \code{rownames} and \code{colnames} of the input data frame.  Names in CEP format are limited to a maximum of 8 characters - if they excede this limit they will be abbreviated using the function \code{make.cepnames} in the package \code{vegan} and a warning issued.  The function invisibly returns a list of original and saved names.
}
\value{
\code{read.CEP} returns a data frame with species as columns and sites as rows. Columns with zero totals (ie. species with no occurrences) are deleted. Column and row names are taken from the species and site names in the CEP file: invalid characters in the names (including spaces) are replaced by "." and names are made unique using \code{\link[base]{make.names}}.

\code{write.CEP} returns a list with two named elements, \code{site.names} and \code{sp.names}, each character matrix of two columns giving the \emph{original} and \emph{saved} site or species names.
}

\note{
The code uses a combination of C and C++, with several function converted to C from FORTRAN using f2c. It works on 32 bit systems but currently generates an error on some 64-bit systems.  Use \code{\link[vegan]{read.cep}} in the package \code{vegan} as an alternative.
}

\author{ 
Steve Juggins
}

\references{
Hill, M.O. (1979a). DECORANA - A FORTRAN program for detrended correspondence analysis and reciprocal averaging. Cornell University, Ithaca, New York. Program manual..

Hill, M.O. (1979b). TWINSPAN - A FORTRAN program for arranging multivariate data in an ordered two-way table by classification of the individuals and attributes. Cornell University, Ithaca, New York. Program manual. 90pp.

ter Braak, C.J.F. & Smilauer, P. (1998) CANOCO Reference Manual and User's Guide to CANOCO for Windows: Software for Canonical Community Ordination (version 4) Microcomputer Power, Ithaca, NY, USA.
}

\seealso{ 
\code{\link[vegan]{read.cep}}.
}

\examples{
\dontrun{
pth <- system.file("example.datasets/RLGHLongCore.cep", package="rioja")
rlgh <- read.CEP(pth)
depth <- as.numeric(gsub("r_", "", rownames(rlgh)))
mx <- apply(rlgh, 2, max)
rlgh.sub <- rlgh[, mx > 5]
strat.plot(rlgh.sub, scale.percent=TRUE, yvar=depth, y.rev=TRUE,
wa.order="bottomleft")
}
}
\keyword{ IO }
\keyword{ file }
