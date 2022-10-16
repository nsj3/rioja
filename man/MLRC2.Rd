\name{MLRC2}
\alias{MLRC2}
\alias{MLRC2.fit}
\alias{predict.MLRC2}
\alias{performance.MLRC2}
\alias{print.MLRC2}
\alias{summary.MLRC2}
\alias{coef.MLRC2}
\alias{fitted.MLRC2}
\alias{residuals.MLRC2}

\title{Palaeoenvironmental reconstruction using Maximum Likelihood Response Surfaces}
\description{
Functions for reconstructing (predicting) environmental values from biological assemblages using Maximum Likelihood response Surfaces.
}
\usage{
MLRC2(y, x, n.out=100, expand.grad=0.1, use.gam=FALSE, check.data=TRUE, 
       lean=FALSE, n.cut=5, verbose=TRUE, \dots)

MLRC2.fit(y, x, n.out=100, expand.grad=0.1, use.gam=FALSE, check.data=TRUE, 
       lean=FALSE, n.cut=5, verbose=TRUE, \dots)

\method{predict}{MLRC2} (object, newdata=NULL, sse=FALSE, nboot=100,
      match.data=TRUE, verbose=TRUE, \dots)

\method{performance}{MLRC2}(object, \dots)

\method{print}{MLRC2}(x, \dots)

\method{summary}{MLRC2}(object, full=FALSE, \dots)

\method{residuals}{MLRC2}(object, cv=FALSE, \dots)

\method{coef}{MLRC2}(object, \dots)

\method{fitted}{MLRC2}(object, \dots)
}
\arguments{
  \item{y}{ a data frame or matrix of biological abundance data. }
  \item{x, object}{ a vector of environmental values to be modelled or an object of class \code{wa}. }
  \item{n.cut}{ cutoff value for number of occurrences.  Species with fewer than n.cut occurrences will be excluded from the analysis. }
  \item{n.out}{ to do }
  \item{expand.grad}{ to do }
  \item{use.gam}{ logical to use \code{gam} to fit responses rather than internal code.  Defaults to \code{FALSE}.}
  \item{newdata}{ new biological data to be predicted. }
  \item{check.data}{ logical to perform simple checks on the input data. }
  \item{match.data}{ logical indicate the function will match two species datasets by their column names.  You should only set this to \code{FALSE} if you are sure the column names match exactly. }
  \item{lean}{ logical to exclude some output from the resulting models (used when cross-validating to speed calculations). }
  \item{full}{ logical to show head and tail of output in summaries. }
  \item{verbose}{ logical to show feedback during cross-validation. }
  \item{nboot}{ number of bootstrap samples. }  
  \item{sse}{ logical indicating that sample specific errors should be calculated. }    
  \item{cv}{ logical to indicate model or cross-validation residuals. }
  \item{\dots}{ additional arguments. }
}

\details{
Function \code{MLRC2} Maximim likelihood reconstruction using 2D response curves.

}

\value{
Function \code{MLRC2} returns an object of class \code{MLRC2} with the following named elements:
}


\author{ 
Steve Juggins
}

\references{
Birks, H.J.B., Line, J.M., Juggins, S., Stevenson, A.C., & ter Braak, C.J.F. (1990) Diatoms and pH reconstruction. \emph{Philosophical Transactions of the Royal Society of London}, \bold{B, 327}, 263-278.

Juggins, S. (1992) Diatoms in the Thames Estuary, England: Ecology, Palaeoecology, and Salinity Transfer Function. \emph{Bibliotheca Diatomologica}, \bold{Band 25}, 216pp.

Oksanen, J., Laara, E., Huttunen, P., & Merilainen, J. (1990) Maximum likelihood prediction of lake acidity based on sedimented diatoms. \emph{Journal of Vegetation Science}, \bold{1}, 49-56.

ter Braak, C.J.F. & van Dam, H. (1989) Inferring pH from diatoms: a comparison of old and new calibration methods. \emph{Hydrobiologia}, \bold{178}, 209-223.
}

\seealso{ 
\code{\link{WA}}, \code{\link{MAT}}, \code{\link{performance}}, and \code{\link{compare.datasets}} for diagnostics.
}


\keyword{ models }
\keyword{ multivariate }
