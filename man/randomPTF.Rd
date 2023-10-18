\name{randomPTF}
\alias{randomPTF}
\alias{plot.randomPTF}
\alias{print.randomPTF}

\title{Random transfer functions to calculate variable importance}
\description{Function for calculating the important of each taxon (predictor) in palaeoecological transfer functions}
\usage{
randomPTF(spec, env, fun, ncol = 1, nVar, nTF = 500, verbose = TRUE, 
          do.parallel = FALSE, \dots)

\method{plot}{randomPTF}(x, use.pointLabel=TRUE, \dots)

\method{print}{randomPTF}(x, \dots)
}

\arguments{
  \item{spec}{ a data frame or matrix of biological abundance data. }
  \item{env}{ a vector of environmental values to be modelled. }
  \item{fun}{ a transfer function method.  Additional arguments can be passed with \dots) }
  \item{ncol}{ some transfer functions return more than one column of results, for example with different WAPLS components. col selects which column to use. See the relevant transfer function method help file. }
  \item{nVar}{ number of variables (ie. species) to use in each randomisation (defaults to nsp/3). }
  \item{nTF}{ number of random transfer functions to create (default=500). }
  \item{verbose}{ logical show feedback during cross-validation. }
  \item{do.parallel}{ logical to run in parallel on multi-core machines.  If true a suitable parallel back-end should be installed (see examples). }
  \item{\dots}{ additional parameters to the transfer function call. }
  \item{x}{ an object of class randomPTF. }
  \item{use.pointLabel}{ argument is deprecated. }
}
\details{
Function \code{randomPTF} calculates taxon importance values using a method analogous to that used in random forests and described in Juggins et al. (2015).

The parallel version can give c. 3 times speed-up on a quad-core machine.
}

\value{
Function \code{randomPTF} returns an object of class \code{randomPTF} with the following named elements:
\item{VI}{ taxon importance values, ordered form high to low. }
\item{spec}{ original species data frame. }
\item{env}{ original vector of environmental values. }
}

\author{
Steve Juggins
}

\references{
Juggins S, Simpson GL, Telford RJ. Taxon selection using statistical learning techniques to improve transfer function prediction. The Holocene 2015; 25: 130-136.}

\examples{
\dontrun{
data(SWAP)
result <- randomPTF(SWAP$spec, SWAP$pH, fun=WA)
plot(result, cex=0.6)
print(result)
# parallel version
if (.Platform$OS.type=='windows') {
  library(doParallel)
  registerDoParallel(cores=4)
} else {
  library(doMC)
  registerDoMC(cores=4)
}
system.time(result <- randomPTF(SWAP$spec, SWAP$pH, fun=WA, do.parallel=TRUE, nTF=5000))}
}

\keyword{ multivariate }
\keyword{ models }
