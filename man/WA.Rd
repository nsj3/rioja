\name{WA}
\alias{WA}
\alias{WA.fit}
\alias{predict.WA}
\alias{crossval.WA}
\alias{performance.WA}
\alias{print.WA}
\alias{summary.WA}
\alias{plot.WA}
\alias{coef.WA}
\alias{fitted.WA}
\alias{residuals.WA}
\alias{rand.t.test.WA}

\title{Weighted averaging (WA) regression and calibration}
\description{
Functions for reconstructing (predicting) environmental values from biological assemblages using weighted averaging (WA) regression and calibration.
}
\usage{
WA(y, x, mono=FALSE, tolDW = FALSE, use.N2=TRUE, tol.cut=.01, 
      check.data=TRUE, lean=FALSE)

WA.fit(y, x, mono=FALSE, tolDW=FALSE, use.N2=TRUE, tol.cut=.01,
       lean=FALSE)

\method{predict}{WA} (object, newdata=NULL, sse=FALSE, nboot=100,
      match.data=TRUE, verbose=TRUE, \dots)

\method{crossval}{WA}(object, cv.method="loo", verbose=TRUE, ngroups=10, 
      nboot=100, h.cutoff=0, h.dist=NULL, \dots)

\method{performance}{WA}(object, \dots)

\method{rand.t.test}{WA}(object, n.perm=999, \dots)

\method{print}{WA}(x, \dots)

\method{summary}{WA}(object, full=FALSE, \dots)

\method{plot}{WA}(x, resid=FALSE, xval=FALSE, tolDW=FALSE, deshrink="inverse",
      xlab="", ylab="", ylim=NULL, xlim=NULL, add.ref=TRUE,
      add.smooth=FALSE, \dots)

\method{residuals}{WA}(object, cv=FALSE, \dots)

\method{coef}{WA}(object, \dots)

\method{fitted}{WA}(object, \dots)
}
\arguments{
  \item{y}{ a data frame or matrix of biological abundance data. }
  \item{x, object}{ a vector of environmental values to be modelled or an object of class \code{WA}. }
  \item{newdata}{ new biological data to be predicted. }
  \item{mono}{ logical to perform monotonic curvilinear deshrinking. }
  \item{tolDW}{ logical to include regressions and predictions using tolerance downweighting. }
  \item{use.N2}{ logical to adjust tolerance by species N2 values. }
  \item{tol.cut}{ tolerances less than \code{tol.cut} are replaced by the mean tolerance. }
  \item{check.data}{ logical to perform simple checks on the input data. }
  \item{lean}{ logical to exclude some output from the resulting models (used when cross-validating to speed calculations). }
  \item{full}{ logical to show head and tail of output in summaries. }
  \item{match.data}{ logical indicate the function will match two species datasets by their column names.  You should only set this to \code{FALSE} if you are sure the column names match exactly. }
  \item{resid}{ logical to plot residuals instead of fitted values. }
  \item{xval}{ logical to plot cross-validation estimates. }
  \item{xlab, ylab, xlim, ylim}{ additional graphical arguments to \code{plot.WA}. }
  \item{deshrink}{ deshrinking type to show in plot. }
  \item{add.ref}{ add 1:1 line on plot. }
  \item{add.smooth}{ add loess smooth to plot. }
  \item{cv.method}{ cross-validation method, either "loo", "lgo", "bootstrap" or "h-block". }
  \item{verbose}{ logical to show feedback during cross-validaton. }
  \item{nboot}{ number of bootstrap samples. }  
  \item{ngroups}{ number of groups in leave-group-out cross-validation. } 
  \item{h.cutoff}{ cutoff for h-block cross-validation.  Only training samples greater than \code{h.cutoff} from each test sample will be used. }
  \item{h.dist}{ distance matrix for use in h-block cross-validation.  Usually a matrix of geographical distances between samples. }
  \item{sse}{ logical indicating that sample specific errors should be calculated. }    
  \item{n.perm}{ number of permutations for randomisation t-test. } 
  \item{cv}{ logical to indicate model or cross-validation residuals. }
  \item{\dots}{ additional arguments. }
}

\details{
Function \code{WA} performs weighted average (WA) regression and calibration. Weighted averaging has a long history in ecology and forms the basis of many biotic indices.  It WAs popularised in palaeolimnology by ter Brakk and van Dam (1989) and Birks et al. (1990) follwoing ter Braak & Barendregt (1986) and ter Braak and Looman (1986) who demonstrated it's theroetical properties in providing a robust and simple alternative to species response modelling using Gaussian logistic regression. Function \code{WA} predicts environmental values from sub-fossil biological assemblages, given a training dataset of modern species and envionmental data. It calculates estimates using inverse and classical deshrinking, and, optionally, with taxa downweighted by their tolerances.  Prediction errors and model complexity (simple or tolerance downweighted WA) can be estimated by cross-validation using \code{crossval} which implements leave-one out, leave-group-out, or bootstrapping. With leave-group out one may also supply a vector of group memberships for more carefully designed cross-validation experiments.

Function \code{predict} predicts values of the environemntal variable for \code{newdata} or returns the fitted (predicted) values from the original modern dataset if \code{newdata} is \code{NULL}. Variables are matched between training and newdata by column name (if \code{match.data} is \code{TRUE}).  Use \code{\link{compare.datasets}} to assess conformity of two species datasets and identify possible no-analogue samples.

Function \code{rand.t.test} performs a randomisation t-test to test the significance of the difference in cross-validation RMSE between tolerance-downweighted and simple WA, after van der Voet (1994).

\code{WA} has methods \code{fitted} and \code{rediduals} that return the fitted values (estimates) and residuals for the training set, \code{performance}, which returns summary performance statistics (see below), \code{coef} which returns the species coefficients (optima and tolerances), and \code{print} and \code{summary} to summarise the output. \code{WA} also has a \code{plot} method that produces scatter plots of predicted vs observed measurements for the training set.
}

\value{
Function \code{WA} returns an object of class \code{WA} with the following named elements:
\item{coefficients}{ species coefficients ("optima" and, optionally, "tolerances"). }
\item{deshrink.coefficients}{ deshrinking coefficients. }
\item{tolDW}{ logical to indicate tolerance downweighted results in model. }
\item{fitted.values}{ fitted values for the training set. }
\item{call}{ original function call. }
\item{x}{ environmental variable used in the model. }

If function \code{predict} is called with \code{newdata=NULL} it returns the fitted values of the original model, otherwise it returns a list with the following named elements:
\item{fit}{ predicted values for \code{newdata}. }

If sample specific errors were requested the list will also include:
\item{fit.boot}{ mean of the bootstrap estimates of newdata. }
\item{v1}{ squared standard error of the bootstrap estimates for each new sample. }
\item{v2}{ mean squared error for the training set samples, across all bootstram samples. }
\item{SEP}{ standard error of prediction, calculated as the square root of v1 + v2. }

Function \code{crossval} also returns an object of class \code{WA} and adds the following named elements:
   \item{predicted}{ predicted values of each training set sample under cross-validation. }
   \item{residuals.cv}{ prediction residuals. }
  
Function \code{performance} returns a matrix of performance statistics for the WA model.  See \code{\link{performance}}, for a description of the summary.
}

\author{ 
Steve Juggins
}

\references{
Birks, H.J.B., Line, J.M., Juggins, S., Stevenson, A.C., & ter Braak, C.J.F. (1990) Diatoms and pH reconstruction. \emph{Philosophical Transactions of the Royal Society of London}, \bold{B, 327}, 263-278.

ter Braak, C.J.F. & Barendregt, L.G. (1986) Weighted averaging of species indicator values: its efficiency in environmental calibration. Mathematical Biosciences, 78, 57-72.

ter Braak, C.J.F. & Looman, C.W.N. (1986) Weighted averaging, logistic regression and the Gaussian response model. \emph{Vegetatio}, \bold{65}, 3-11.

ter Braak, C.J.F. & van Dam, H. (1989) Inferring pH from diatoms: a comparison of old and new calibration methods. \emph{Hydrobiologia}, \bold{178}, 209-223.

van der Voet, H. (1994) Comparing the predictive accuracy of models uing a simple randomization test. \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{25}, 313-323.
}

\seealso{ 
\code{\link{WAPLS}}, \code{\link{MAT}}, and \code{\link{compare.datasets}} for diagnostics.
}

\examples{
# pH reconstruction of core K05 from the Round Loch of Glenhead,
# Galloway, SW Scotland. This lake has become acidified over the 
# last c. 150 years

data(SWAP)
data(RLGH)
spec <- SWAP$spec
pH <- SWAP$pH
core <- RLGH$spec
age <- RLGH$depths$Age

fit <- WA(spec, pH, tolDW=TRUE)
# plot predicted vs. observed
plot(fit)
plot(fit, resid=TRUE)

# RLGH reconstruction
pred <- predict(fit, core)

#plot the reconstructio
plot(age, pred$fit[, 1], type="b")

# cross-validation model using bootstrapping
\dontrun{
fit.xv <- crossval(fit, cv.method="boot", nboot=1000)
par(mfrow=c(1,2))
plot(fit)
plot(fit, resid=TRUE)
plot(fit.xv, xval=TRUE)
plot(fit.xv, xval=TRUE, resid=TRUE)

# RLGH reconstruction with sample specific errors
pred <- predict(fit, core, sse=TRUE, nboot=1000)
}
}

\keyword{ models }
\keyword{ multivariate }
