\name{MLRC}
\alias{MLRC}
\alias{MLRC.fit}
\alias{predict.MLRC}
\alias{crossval.MLRC}
\alias{performance.MLRC}
\alias{print.MLRC}
\alias{summary.MLRC}
\alias{plot.MLRC}
\alias{coef.MLRC}
\alias{fitted.MLRC}
\alias{residuals.MLRC}

\title{Palaeoenvironmental reconstruction using Maximum Likelihood Response Surfaces}
\description{
Functions for reconstructing (predicting) environmental values from biological assemblages using Maximum Likelihood response Surfaces.
}
\usage{
MLRC(y, x, check.data=TRUE, lean=FALSE, n.cut=5, verbose=TRUE, \dots)

MLRC.fit(y, x, n.cut=2, use.glm=FALSE, max.iter=50, lean=FALSE, verbose=TRUE, \dots)

\method{predict}{MLRC} (object, newdata=NULL, sse=FALSE, nboot=100,
      match.data=TRUE, verbose=TRUE, \dots)

\method{crossval}{MLRC}(object, cv.method="loo", verbose=TRUE, ngroups=10,
      nboot=100, h.cutoff=0, h.dist=NULL, \dots)

\method{performance}{MLRC}(object, \dots)

\method{print}{MLRC}(x, \dots)

\method{summary}{MLRC}(object, full=FALSE, \dots)

\method{plot}{MLRC}(x, resid=FALSE, xval=FALSE, xlab="", ylab="", 
      ylim=NULL, xlim=NULL, add.ref=TRUE, add.smooth=FALSE, \dots)

\method{residuals}{MLRC}(object, cv=FALSE, \dots)

\method{coef}{MLRC}(object, \dots)

\method{fitted}{MLRC}(object, \dots)
}
\arguments{
  \item{y}{ a data frame or matrix of biological abundance data. }
  \item{x, object}{ a vector of environmental values to be modelled or an object of class \code{wa}. }
  \item{n.cut}{ cutoff value for number of occurrences.  Species with fewer than n.cut occurrences will be excluded from the analysis.}
  \item{use.glm}{ logical to use \code{glm} to fit responses rather than internal code.  Defaults to \code{FALSE}.}
  \item{newdata}{ new biological data to be predicted. }
  \item{max.iter}{ maximum iterations of the logit regression algorithm.}
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
  \item{h.cutoff}{ cutoff for h-block cross-validation.  Only training samples greater than \code{h.cutoff} from each test sample will be used. }
  \item{h.dist}{ distance matrix for use in h-block cross-validation.  Usually a matrix of geographical distances between samples. }
  \item{sse}{ logical indicating that sample specific errors should be calculated. }    
  \item{cv}{ logical to indicate model or cross-validation residuals. }
  \item{\dots}{ additional arguments. }
}

\details{
Function \code{MLRC} Maximim likelihood reconstruction using response curves.

Function \code{predict} predicts values of the environemntal variable for \code{newdata} or returns the fitted (predicted) values from the original modern dataset if \code{newdata} is \code{NULL}. Variables are matched between training and newdata by column name (if \code{match.data} is \code{TRUE}).  Use \code{\link{compare.datasets}} to assess conformity of two species datasets and identify possible no-analogue samples.

\code{MLRC} has methods \code{fitted} and \code{rediduals} that return the fitted values (estimates) and residuals for the training set, \code{performance}, which returns summary performance statistics (see below), \code{coef} which returns the species coefficients, and \code{print} and \code{summary} to summarise the output. \code{MLRC} also has a \code{plot} method that produces scatter plots of predicted vs observed measurements for the training set.
}

\value{
Function \code{MLRC} returns an object of class \code{MLRC} with the following named elements:


Function \code{crossval} also returns an object of class \code{MLRC} and adds the following named elements:
\item{predicted}{ predicted values of each training set sample under cross-validation. }
\item{residuals.cv}{ prediction residuals. }

If function \code{predict} is called with \code{newdata=NULL} it returns the fitted values of the original model, otherwise it returns a list with the following named elements:
\item{fit}{ predicted values for \code{newdata}. }
if sample specific errors were generated the list will also include:
\item{fit.boot}{ mean of the bootstrap estimates of newdata. }
\item{v1}{ squared standard error of the bootstrap estimates for each new sample. }
\item{v2}{ mean squared error for the training set samples, across all bootstrap samples. }
\item{SEP}{ standard error of prediction, calculated as the square root of v1 + v2. }

Function \code{performance} returns a matrix of performance statistics for the MLRC model.  See \code{\link{performance}}, for a description of the summary.
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

\examples{
data(IK)
spec <- IK$spec / 100
SumSST <- IK$env$SumSST
core <- IK$core / 100

fit <- MLRC(spec, SumSST)
fit

#predict the core
pred <- predict(fit, core)

#plot predictions - depths are in rownames
depth <- as.numeric(rownames(core))
plot(depth, pred$fit[, 1], type="b")

\dontrun{
# this is slow!
# cross-validate model
fit.cv <- crossval(fit, cv.method="loo", verbose=5)

# predictions with sample specific errors
pred <- predict(fit, core, sse=TRUE, nboot=1000, verbose=5)
}
}
\keyword{ models }
\keyword{ multivariate }
