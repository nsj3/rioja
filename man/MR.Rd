\name{MR}
\alias{MR}
\alias{MR.fit}
\alias{predict.MR}
\alias{crossval.MR}
\alias{performance.MR}
\alias{print.MR}
\alias{summary.MR}
\alias{plot.MR}
\alias{coef.MR}
\alias{fitted.MR}
\alias{residuals.MR}

\title{Multiple regression}
\description{
Functions for reconstructing (predicting) environmental values from biological assemblages using multiple regression.
}
\usage{
MR(y, x, check.data=TRUE, lean=FALSE, \dots)

MR.fit(y, x, lean=FALSE)

\method{predict}{MR} (object, newdata=NULL, sse=FALSE, nboot=100,
      match.data=TRUE, verbose=TRUE, \dots)

\method{crossval}{MR}(object, cv.method="loo", verbose=TRUE, ngroups=10,
      nboot=100, h.cutoff=0, h.dist=NULL, \dots)

\method{performance}{MR}(object, \dots)

\method{print}{MR}(x, \dots)

\method{summary}{MR}(object, full=FALSE, \dots)

\method{plot}{MR}(x, resid=FALSE, xval=FALSE, xlab="", 
      ylab="", ylim=NULL, xlim=NULL, add.ref=TRUE,
      add.smooth=FALSE, \dots)

\method{residuals}{MR}(object, cv=FALSE, \dots)

\method{coef}{MR}(object, \dots)

\method{fitted}{MR}(object, \dots)

}
\arguments{
  \item{y}{ a data frame or matrix of biological abundance data. }
  \item{x, object}{ a vector of environmental values to be modelled or an object of class \code{wa}. }
  \item{newdata}{ new biological data to be predicted. }
  \item{check.data}{ logical to perform simple checks on the input data. }
  \item{match.data}{ logical indicate the function will match two species datasets by their column names.  You should only set this to \code{FALSE} if you are sure the column names match exactly. }
  \item{lean}{ logical to exclude some output from the resulting models (used when cross-validating to speed calculations). }
  \item{full}{ logical to show head and tail of output in summaries. }
  \item{resid}{ logical to plot residuals instead of fitted values. }
  \item{xval}{ logical to plot cross-validation estimates. }
  \item{xlab, ylab, xlim, ylim}{ additional graphical arguments to \code{plot.wa}. }
  \item{add.ref}{ add 1:1 line on plot. }
  \item{add.smooth}{ add loess smooth to plot. }
  \item{cv.method}{ cross-validation method, either "loo", "lgo", "bootstrap" or "h-block". }
  \item{verbose}{ logical to show feedback during cross-validation. }
  \item{nboot}{ number of bootstrap samples. }  
  \item{ngroups}{ number of groups in leave-group-out cross-validation, or a vector contain leave-out group menbership. }  
  \item{h.cutoff}{ cutoff for h-block cross-validation. Only training samples greater than \code{h.cutoff} from each test sample will be used. }
  \item{h.dist}{ distance matrix for use in h-block cross-validation. Usually a matrix of geographical distances between samples. }
  \item{sse}{ logical indicating that sample specific errors should be calculated. }    
  \item{cv}{ logical to indicate model or cross-validation residuals. }
  \item{\dots}{ additional arguments. }
}

\details{
Function \code{MR} performs multiple regrssion.  It is a wrapper to \code{lm}. 

Function \code{predict} predicts values of the environemntal variable for \code{newdata} or returns the fitted (predicted) values from the original modern dataset if \code{newdata} is \code{NULL}. Variables are matched between training and newdata by column name (if \code{match.data} is \code{TRUE}).  Use \code{\link{compare.datasets}} to assess conformity of two species datasets and identify possible no-analogue samples. 

\code{MR} has methods \code{fitted} and \code{rediduals} that return the fitted values (estimates) and residuals for the training set, \code{performance}, which returns summary performance statistics (see below), \code{coef} which returns the species coefficients, and \code{print} and \code{summary} to summarise the output. \code{MR} also has a \code{plot} method that produces scatter plots of predicted vs observed measurements for the training set. }

\value{
Function \code{MR} returns an object of class \code{MR} with the following named elements:
\item{coefficients}{ species coefficients (the updated "optima"). }
\item{fitted.values}{ fitted values for the training set. }
\item{call}{ original function call. }
\item{x}{ environmental variable used in the model. }

Function \code{crossval} also returns an object of class \code{MR} and adds the following named elements:
\item{predicted}{ predicted values of each training set sample under cross-validation. }
\item{residuals.cv}{ prediction residuals. }

If function \code{predict} is called with \code{newdata=NULL} it returns the fitted values of the original model, otherwise it returns a list with the following named elements:
\item{fit}{ predicted values for \code{newdata}. }
if sample specific errors were generated the list will also include:
\item{fit.boot}{ mean of the bootstrap estimates of newdata. }
\item{v1}{ squared standard error of the bootstrap estimates for each new sample. }
\item{v2}{ mean squared error for the training set samples, across all bootstrap samples. }
\item{SEP}{ standard error of prediction, calculated as the square root of v1 + v2. }

Function \code{performance} returns a matrix of performance statistics for the MR model.  See \code{\link{performance}}, for a description of the summary.
}

\author{ 
Steve Juggins
}

\seealso{ 
\code{\link{WA}}, \code{\link{MAT}}, \code{\link{performance}}, and \code{\link{compare.datasets}} for diagnostics.
}

\examples{
data(IK)
spec <- IK$spec
SumSST <- IK$env$SumSST
core <- IK$core

fit <- IKFA(spec, SumSST)
fit
# cross-validate model
fit.cv <- crossval(fit, cv.method="lgo")
# How many components to use?
screeplot(fit.cv)

#predict the core
pred <- predict(fit, core, npls=2)

#plot predictions - depths are in rownames
depth <- as.numeric(rownames(core))
plot(depth, pred$fit[, 2], type="b")

# fit using only factors 1, 2, 4, & 5
# and using polynomial terms
# as Imbrie & Kipp (1971)
fit2 <- IKFA(spec, SumSST, ccoef=c(1, 2, 4, 5), IsPoly=TRUE)
fit2.cv <- crossval(fit2, cv.method="lgo")
screeplot(fit2.cv)

\dontrun{
# predictions with sample specific errors
# takes approximately 1 minute to run
pred <- predict(fit, core, sse=TRUE, nboot=1000)
pred
}
}
\keyword{ models }
\keyword{ multivariate }
