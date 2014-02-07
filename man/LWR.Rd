\name{LWR}
\alias{LWR}
\alias{LWR.fit}
\alias{predict.LWR}
\alias{crossval.LWR}
\alias{performance.LWR}
\alias{print.LWR}
\alias{summary.LWR}
\alias{plot.LWR}
\alias{coef.LWR}
\alias{fitted.LWR}
\alias{residuals.LWR}

\title{Weighted averaging (LWR) regression and calibration}
\description{
Functions for reconstructing (predicting) environmental values from biological assemblages using weighted averaging (LWR) regression and calibration.
}
\usage{
LWR(y, x, FUN=WA, dist.method="sq.chord", k=30, lean=TRUE, 
       fit.model=TRUE, check.data=TRUE, verbose=TRUE, \dots)

\method{predict}{LWR} (object, newdata=NULL, k = object$k, sse=FALSE, 
      nboot=100, match.data=TRUE, verbose=TRUE, lean=TRUE, \dots)

\method{crossval}{LWR}(object, k=object$k, cv.method="lgo", verbose=TRUE, 
      ngroups=10, nboot=100, h.cutoff=0, h.dist=NULL, \dots)

\method{performance}{LWR}(object, \dots)

\method{print}{LWR}(x, \dots)

\method{summary}{LWR}(object, full=FALSE, \dots)

\method{residuals}{LWR}(object, cv=FALSE, \dots)

\method{fitted}{LWR}(object, \dots)
}
\arguments{
  \item{y}{ a data frame or matrix of biological abundance data. }
  \item{x, object}{ a vector of environmental values to be modelled or an object of class \code{LWR}. }
  \item{dist.method}{ distance measure used to derfine closest analogues. }
  \item{k}{ number of close analogues to use in calibration function. }
  \item{FUN}{ calibration function (e.g. \code{WA}, \code{WAPLS} etc). }
  \item{newdata}{ new biological data to be predicted. }
  \item{fit.model}{ TRUE fits model to training set.  FALSE omist this step and builds a LWR object than can be used for prediction. }
  \item{check.data}{ logical to perform simple checks on the input data. }
  \item{full}{ logical to show head and tail of output in summaries. }
  \item{match.data}{ logical indicate the function will match two species datasets by their column names.  You should only set this to \code{FALSE} if you are sure the column names match exactly. }
  \item{lean}{ logical to exclude some output from the resulting models (used when cross-validating to speed calculations). }
  \item{cv.method}{ cross-validation method, either "lgo" or "bootstrap". }
  \item{verbose}{ logical to show feedback during cross-validaton. }
  \item{nboot}{ number of bootstrap samples. }  
  \item{ngroups}{ number of groups in leave-group-out cross-validation. }
  \item{h.cutoff}{ cutoff for h-block cross-validation.  Only training samples greater than \code{h.cutoff} from each test sample will be used. }
  \item{h.dist}{ distance matrix for use in h-block cross-validation.  Usually a matrix of geographical distances between samples. }
  \item{sse}{ logical indicating that sample specific errors should be calculated. }    
  \item{cv}{ logical to indicate model or cross-validation residuals. }
  \item{\dots}{ additional arguments. }
}

\details{
Function \code{LWR} performs ... To do.
}

\value{
Function \code{LWR} returns an object of class \code{LWR} with the following named elements:
}

\author{ 
Steve Juggins
}


\seealso{ 
\code{\link{WAPLS}}, \code{\link{MAT}}, and \code{\link{compare.datasets}} for diagnostics.
}


\keyword{ models }
\keyword{ multivariate }
