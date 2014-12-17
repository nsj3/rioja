\name{IKFA}
\alias{IKFA}
\alias{IKFA.fit}
\alias{predict.IKFA}
\alias{crossval.IKFA}
\alias{performance.IKFA}
\alias{print.IKFA}
\alias{summary.IKFA}
\alias{plot.IKFA}
\alias{coef.IKFA}
\alias{fitted.IKFA}
\alias{residuals.IKFA}
\alias{screeplot.IKFA}
\alias{rand.t.test.IKFA}
\alias{communality}

\title{Imbrie & Kipp Factor Analysis}
\description{
Functions for reconstructing (predicting) environmental values from biological assemblages using Imbrie & Kipp Factor Analysis (IKFA), as used in palaeoceanography.
}
\usage{
IKFA(y, x, nFact = 5, IsPoly = FALSE, IsRot = TRUE, 
      ccoef = 1:nFact, check.data=TRUE, lean=FALSE, \dots)

IKFA.fit(y, x, nFact = 5, IsPoly = FALSE, IsRot = TRUE, 
      ccoef = 1:nFact, lean=FALSE)

\method{predict}{IKFA} (object, newdata=NULL, sse=FALSE, nboot=100,
      match.data=TRUE, verbose=TRUE, \dots)

communality(object, y)

\method{crossval}{IKFA}(object, cv.method="loo", verbose=TRUE, ngroups=10,
      nboot=100, h.cutoff=0, h.dist=NULL, \dots)

\method{performance}{IKFA}(object, \dots)

\method{rand.t.test}{IKFA}(object, n.perm=999, \dots)

\method{screeplot}{IKFA}(x, rand.test=TRUE, \dots)

\method{print}{IKFA}(x, \dots)

\method{summary}{IKFA}(object, full=FALSE, \dots)

\method{plot}{IKFA}(x, resid=FALSE, xval=FALSE, nFact=max(x$ccoef), 
      xlab="", ylab="", ylim=NULL, xlim=NULL, add.ref=TRUE,
      add.smooth=FALSE, \dots)

\method{residuals}{IKFA}(object, cv=FALSE, \dots)

\method{coef}{IKFA}(object, \dots)

\method{fitted}{IKFA}(object, \dots)

}
\arguments{
  \item{y}{ a data frame or matrix of biological abundance data. }
  \item{x, object}{ a vector of environmental values to be modelled or an object of class \code{wa}. }
  \item{newdata}{ new biological data to be predicted. }
  \item{nFact}{ number of factor to extract.}
  \item{IsRot}{ logical to rotate factors. }
  \item{ccoef}{ vector of factor numbers to include in the predictions. }
  \item{IsPoly}{ logical to include quadratic of the factors as predictors in the regression. }
  \item{check.data}{ logical to perform simple checks on the input data. }
  \item{match.data}{ logical indicate the function will match two species datasets by their column names.  You should only set this to \code{FALSE} if you are sure the column names match exactly. }
  \item{lean}{ logical to exclude some output from the resulting models (used when cross-validating to speed calculations). }
  \item{full}{ logical to show head and tail of output in summaries. }
  \item{resid}{ logical to plot residuals instead of fitted values. }
  \item{xval}{ logical to plot cross-validation estimates. }
  \item{xlab, ylab, xlim, ylim}{ additional graphical arguments to \code{plot.wa}. }
  \item{add.ref}{ add 1:1 line on plot. }
  \item{add.smooth}{ add loess smooth to plot. }
  \item{cv.method}{ cross-validation method, either "loo", "lgo" or "bootstrap". }
  \item{verbose}{ logical to show feedback during cross-validaton. }
  \item{nboot}{ number of bootstrap samples. }  
  \item{ngroups}{ number of groups in leave-group-out cross-validation, or a vector contain leave-out group menbership. }  
  \item{h.cutoff}{ cutoff for h-block cross-validation. Only training samples greater than \code{h.cutoff} from each test sample will be used. }
  \item{h.dist}{ distance matrix for use in h-block cross-validation. Usually a matrix of geographical distances between samples. }
  \item{sse}{ logical indicating that sample specific errors should be calculated. }    
  \item{rand.test}{ logical to perform a randomisation t-test to test significance of cross validated factors. } 
  \item{n.perm}{ number of permutations for randomisation t-test. } 
  \item{cv}{ logical to indicate model or cross-validation residuals. }
  \item{\dots}{ additional arguments. }
}

\details{
Function \code{IKFA} performs Imbrie and Kipp Factor Analysis, a form of Principal Components Regrssion (Imbrie & Kipp 1971). 

Function \code{predict} predicts values of the environemntal variable for \code{newdata} or returns the fitted (predicted) values from the original modern dataset if \code{newdata} is \code{NULL}. Variables are matched between training and newdata by column name (if \code{match.data} is \code{TRUE}).  Use \code{\link{compare.datasets}} to assess conformity of two species datasets and identify possible no-analogue samples.

\code{IKFA} has methods \code{fitted} and \code{rediduals} that return the fitted values (estimates) and residuals for the training set, \code{performance}, which returns summary performance statistics (see below), \code{coef} which returns the species coefficients, and \code{print} and \code{summary} to summarise the output. \code{IKFA} also has a \code{plot} method that produces scatter plots of predicted vs observed measurements for the training set.

Function \code{\link{rand.t.test}} performs a randomisation t-test to test the significance of the cross-validated components after van der Voet (1994).

Function \code{screeplot} displays the RMSE of prediction for the training set as a function of the number of factors and is useful for estimating the optimal number for use in prediction. By default \code{screeplot} will also carry out a randomisation t-test and add a line to scree plot indicating percentage change in RMSE with each component annotate with the p-value from the randomisation test.
}

\value{
Function \code{IKFA} returns an object of class \code{IKFA} with the following named elements:
\item{coefficients}{ species coefficients (the updated "optima"). }
\item{fitted.values}{ fitted values for the training set. }
\item{call}{ original function call. }
\item{x}{ environmental variable used in the model. }
\item{standx, meanT sdx}{ additional information returned for a PLSif model. }

Function \code{crossval} also returns an object of class \code{IKFA} and adds the following named elements:
\item{predicted}{ predicted values of each training set sample under cross-validation. }
\item{residuals.cv}{ prediction residuals. }

If function \code{predict} is called with \code{newdata=NULL} it returns the fitted values of the original model, otherwise it returns a list with the following named elements:
\item{fit}{ predicted values for \code{newdata}. }
if sample specific errors were generated the list will also include:
\item{fit.boot}{ mean of the bootstrap estimates of newdata. }
\item{v1}{ squared standard error of the bootstrap estimates for each new sample. }
\item{v2}{ mean squared error for the training set samples, across all bootstrap samples. }
\item{SEP}{ standard error of prediction, calculated as the square root of v1 + v2. }

Function \code{performance} returns a matrix of performance statistics for the IKFA model.  See \code{\link{performance}}, for a description of the summary.

Function \code{\link{rand.t.test}} returns a matrix of performance statistics together with columns indicating the p-value and percentage change in RMSE with each higher component (see van der Veot (1994) for details).
}

\author{ 
Steve Juggins
}

\references{
Imbrie, J. & Kipp, N.G. (1971). A new micropaleontological method for quantitative paleoclimatology: application to a Late Pleistocene Caribbean core. In \emph{The Late Cenozoic Glacial Ages} (ed K.K. Turekian), pp. 77-181. Yale University Press, New Haven.

van der Voet, H. (1994) Comparing the predictive accuracy of models uing a simple randomization test. \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{25}, 313-323.
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
