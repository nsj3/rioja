\name{WAPLS}
\alias{WAPLS}
\alias{WAPLS.fit}
\alias{predict.WAPLS}
\alias{crossval.WAPLS}
\alias{performance.WAPLS}
\alias{print.WAPLS}
\alias{summary.WAPLS}
\alias{plot.WAPLS}
\alias{coef.WAPLS}
\alias{fitted.WAPLS}
\alias{residuals.WAPLS}
\alias{screeplot.WAPLS}
\alias{rand.t.test.WAPLS}

\title{Weighted averaging partial least squares (WAPLS) regression and calibration}
\description{
Functions for reconstructing (predicting) environmental values from biological assemblages using weighted averaging partial least squares (WAPLS) regression and calibration.
}
\usage{
WAPLS(y, x, npls=5, iswapls=TRUE, standx=FALSE, lean=FALSE,
      check.data=TRUE, \dots)

WAPLS.fit(y, x, npls=5, iswapls=TRUE, standx=FALSE, lean=FALSE)

\method{predict}{WAPLS} (object, newdata=NULL, sse=FALSE, nboot=100,
      match.data=TRUE, verbose=TRUE, \dots)

\method{crossval}{WAPLS}(object, cv.method="loo", verbose=TRUE, ngroups=10, 
      nboot=100, h.cutoff=0, h.dist=NULL, \dots)

\method{performance}{WAPLS}(object, \dots)

\method{rand.t.test}{WAPLS}(object, n.perm=999, \dots)

\method{screeplot}{WAPLS}(x, rand.test=TRUE, \dots)

\method{print}{WAPLS}(x, \dots)

\method{summary}{WAPLS}(object, full=FALSE, \dots)

\method{plot}{WAPLS}(x, resid=FALSE, xval=FALSE, npls=1, 
      xlab="", ylab="", ylim=NULL, xlim=NULL, add.ref=TRUE,
      add.smooth=FALSE, \dots)

\method{residuals}{WAPLS}(object, cv=FALSE, \dots)

\method{coef}{WAPLS}(object, \dots)

\method{fitted}{WAPLS}(object, \dots)
}
\arguments{
  \item{y}{ a data frame or matrix of biological abundance data. }
  \item{x, object}{ a vector of environmental values to be modelled or an object of class \code{wa}. }
  \item{newdata}{ new biological data to be predicted. }
  \item{iswapls}{ logical logical to perform WAPLS or PLS. Defaults to TRUE = WAPLS.}
  \item{standx}{ logical to standardise x-data in PLS, defaults to FALSE. }
  \item{npls}{ number of pls components to extract. }
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
  \item{verbose}{ logical show feedback during cross-validation. }
  \item{nboot}{ number of bootstrap samples. }  
  \item{ngroups}{ number of groups in leave-group-out cross-validation, or a vector contain leave-out group menbership. } 
  \item{h.cutoff}{ cutoff for h-block cross-validation.  Only training samples greater than \code{h.cutoff} from each test sample will be used. }
  \item{h.dist}{ distance matrix for use in h-block cross-validation.  Usually a matrix of geographical distances between samples. }
  \item{sse}{ logical indicating that sample specific errors should be calculated. }  
  \item{rand.test}{ logical to perform a randomisation t-test to test significance of cross validated components. } 
  \item{n.perm}{ number of permutations for randomisation t-test. }  
  \item{cv}{ logical to indicate model or cross-validation residuals. }
  \item{\dots}{ additional arguments. }
}

\details{
Function \code{WAPLS} performs partial least squares (PLS) or weighted averaging partial least squares (WAPLS) regression. WAPLS was first described in ter Braak and Juggins (1993) and ter Braak et al. (1993) and has since become popular in palaeolimnology for reconstructing (predicting) environmental values from sub-fossil biological assemblages, given a training dataset of modern species and envionmental data. Prediction errors and model complexity (number of components) can be estimated by cross-validation using \code{crossval} which implements leave-one out, leave-group-out, or bootstrapping. With leave-group out one may also supply a vector of group memberships for more carefully designed cross-validation experiments.

Function \code{predict} predicts values of the environemntal variable for \code{newdata} or returns the fitted (predicted) values from the original modern dataset if \code{newdata} is \code{NULL}. Variables are matched between training and newdata by column name (if \code{match.data} is \code{TRUE}).  Use \code{\link{compare.datasets}} to assess conformity of two species datasets and identify possible no-analogue samples.

\code{WAPLS} has methods \code{fitted} and \code{rediduals} that return the fitted values (estimates) and residuals for the training set, \code{performance}, which returns summary performance statistics (see below), \code{coef} which returns the species coefficients, and \code{print} and \code{summary} to summarise the output. \code{WAPLS} also has a \code{plot} method that produces scatter plots of predicted vs observed measurements for the training set.

Function \code{rand.t.test} performs a randomisation t-test to test the significance of the cross-validated components after van der Voet (1994).

Function \code{screeplot} displays the RMSE of prediction for the training set as a function of the number of components and is useful for estimating the optimal number for use in prediction. By default \code{screeplot} will also carry out a randomisation t-test and add a line to scree plot indicating percentage change in RMSE with each component annotate with the p-value from the randomisation test.
}

\value{
Function \code{WAPLS} returns an object of class \code{WAPLS} with the following named elements:
\item{coefficients}{ species coefficients (the updated "optima"). }
\item{meanY}{ weighted mean of the environmental variable. }
\item{iswapls}{ logical indicating whether analysis was WAPLS (TRUE) or PLS (FALSE). }
\item{T}{ sample scores. }
\item{P}{ variable (species) scores. }
\item{npls}{ number of pls components extracted. }
\item{fitted.values}{ fitted values for the training set. }
\item{call}{ original function call. }
\item{x}{ environmental variable used in the model. }
\item{standx, meanT sdx}{ additional information returned for a PLS model. }

Function \code{crossval} also returns an object of class \code{WAPLS} and adds the following named elements:
\item{predicted}{ predicted values of each training set sample under cross-validation. }
\item{residuals.cv}{ prediction residuals. }

If function \code{predict} is called with \code{newdata=NULL} it returns the fitted values of the original model, otherwise it returns a list with the following named elements:
\item{fit}{ predicted values for \code{newdata}. }
if sample specific errors were generated the list will also include:
\item{fit.boot}{ mean of the bootstrap estimates of newdata. }
\item{v1}{ squared standard error of the bootstrap estimates for each new sample. }
\item{v2}{ mean squared error for the training set samples, across all bootstrap samples. }
\item{SEP}{ standard error of prediction, calculated as the square root of v1 + v2. }

Function \code{performance} returns a matrix of performance statistics for the WAPLS model.  See \code{\link{performance}}, for a description of the summary.

Function \code{\link{rand.t.test}} returns a matrix of performance statistics together with columns indicating the p-value and percentage change in RMSE with each higher component (see van der Veot (1994) for details).
}

\author{ 
Steve Juggins
}

\references{
ter Braak, C.J.F. & Juggins, S. (1993) Weighted averaging partial least squares regression (WA-PLS): an improved method for reconstructing environmental variables from species assemblages. \emph{Hydrobiologia}, \bold{269/270}, 485-502.

ter Braak, C.J.F., Juggins, S., Birks, H.J.B., & Voet, H., van der (1993). Weighted averaging partial least squares regression (WA-PLS): definition and comparison with other methods for species-environment calibration. In \emph{Multivariate Environmental Statistics} (eds G.P. Patil & C.R. Rao), pp. 525-560. Elsevier Science Publishers.

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

fit <- WAPLS(spec, SumSST)
fit
# cross-validate model
fit.cv <- crossval(fit, cv.method="loo")
# How many components to use?
rand.t.test(fit.cv)
screeplot(fit.cv)

#predict the core
pred <- predict(fit, core, npls=2)

#plot predictions - depths are in rownames
depth <- as.numeric(rownames(core))
plot(depth, pred$fit[, 2], type="b", ylab="Predicted SumSST", las=1)

# predictions with sample specific errors
\dontrun{
pred <- predict(fit, core, npls=2, sse=TRUE, nboot=1000)
pred
}
}
\keyword{ models }
\keyword{ multivariate }
