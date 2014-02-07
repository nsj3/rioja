\name{MAT}
\alias{MAT}
\alias{predict.MAT}
\alias{paldist}
\alias{paldist2}
\alias{performance.MAT}
\alias{crossval.MAT}
\alias{print.MAT}
\alias{summary.MAT}
\alias{plot.MAT}
\alias{fitted.MAT}
\alias{residuals.MAT}
\alias{screeplot.MAT}

\title{Palaeoenvironmental reconstruction using the Modern Analogue Technique (MAT)}
\description{
Functions for reconstructing (predicting) environmental values from biological assemblages using the Modern Analogue Technique (MAT), also know as k nearest neighbours (k-NN).
}
\usage{
MAT(y, x, dist.method="sq.chord", k=5, lean=TRUE)

\method{predict}{MAT}(object, newdata=NULL, k=object$k, sse=FALSE, 
        nboot=100, match.data=TRUE, verbose=TRUE, lean=TRUE, 
        \dots)

\method{performance}{MAT}(object, \dots)

\method{crossval}{MAT}(object, k=object$k, cv.method="lgo", 
        verbose=TRUE, ngroups=10, nboot=100, h.cutoff=0, h.dist=NULL, \dots)

\method{print}{MAT}(x, \dots)

\method{summary}{MAT}(object, full=FALSE, \dots)

\method{plot}{MAT}(x, resid=FALSE, xval=FALSE, k=5, wMean=FALSE, xlab="", 
      ylab="", ylim=NULL, xlim=NULL, add.ref=TRUE,
      add.smooth=FALSE, \dots)

\method{residuals}{MAT}(object, cv=FALSE, \dots)

\method{fitted}{MAT}(object, \dots)

\method{screeplot}{MAT}(x, \dots)

paldist(y, dist.method="sq.chord")

paldist2(y1, y2, dist.method="sq.chord")
}
\arguments{
  \item{y, y1, y2}{ data frame containing biological data. }
  \item{newdata}{ data frame containing biological data to predict from. }
  \item{x}{ a vector of environmental values to be modelled, matched to y. }
  \item{dist.method}{ dissimilarity coefficient. See details for options.}
  \item{match.data}{ logical indicate the function will match two species datasets by their column names.  You should only set this to \code{FALSE} if you are sure the column names match exactly. }
  \item{k}{ number of analogues to use. }
  \item{lean}{ logical to remove items form the output. }
  \item{object}{ an object of class \code{MAT}. }
  \item{resid}{ logical to plot residuals instead of fitted values. }
  \item{xval}{ logical to plot cross-validation estimates. }
  \item{wMean}{ logical to plot weighted-mean estimates. }
  \item{xlab, ylab, xlim, ylim}{ additional graphical arguments to \code{plot.wa}. }
  \item{add.ref}{ add 1:1 line on plot. }
  \item{add.smooth}{ add loess smooth to plot. }
  \item{cv.method}{ cross-validation method, either "lgo" or "bootstrap". }
  \item{verbose}{ logical to show feedback during cross-validaton. }
  \item{nboot}{ number of bootstrap samples. }  
  \item{ngroups}{ number of groups in leave-group-out cross-validation, or a vector contain leave-out group menbership. }  
  \item{h.cutoff}{ cutoff for h-block cross-validation.  Only training samples greater than \code{h.cutoff} from each test sample will be used. }
  \item{h.dist}{ distance matrix for use in h-block cross-validation.  Usually a matrix of geographical distances between samples. }
  \item{sse}{ logical indicating that sample specific errors should be calculated. }  
  \item{full}{ logical to indicate a full or abbreviated summary. }
  \item{cv}{ logical to indicate model or cross-validation residuals. }
  \item{\dots}{ additional arguments. }
}

\details{
\code{MAT} performs an environmental reconstruction using the modern analogue technique. Function \code{MAT} takes a training dataset of biological data (species abundances) \code{y} and a single associated environmental variable \code{x}, and generates a model of closest analogues, or matches, for the modern data data using one of a number of dissimilarity coefficients. Options for the latter are: "euclidean", "sq.euclidean", "chord", "sq.chord", "chord.t", "sq.chord.t", "chi.squared", "sq.chi.squared", "bray". "chord.t" are true chord distances, "chord" refers to the the variant of chord distance using in palaeoecology (e.g. Overpeck et al. 1985), which is actually Hellinger's distance (Legendre & Gallagher 2001). There are various help functions to plot and extract information from the results of a \code{MAT} transfer function. The function \code{predict} takes \code{MAT} object and uses it to predict environmental values for a new set of species data, or returns the fitted (predicted) values from the original modern dataset if \code{newdata} is \code{NULL}. Variables are matched between training and newdata by column name (if \code{match.data} is \code{TRUE}).  Use \code{\link{compare.datasets}} to assess conformity of two species datasets and identify possible no-analogue samples.

\code{MAT} has methods \code{fitted} and \code{rediduals} that return the fitted values (estimates) and residuals for the training set, \code{performance}, which returns summary performance statistics (see below), and \code{print} and \code{summary} to summarise the output. \code{MAT} also has a \code{plot} method that produces scatter plots of predicted vs observed measurements for the training set.

Function \code{screeplot} displays the RMSE of prediction for the training set as a function of the number of analogues (k) and is useful for estimating the optimal value of k for use in prediction. 

\code{paldist} and \code{paldist1} are helper functions though they may be called directly. \code{paldist} takes a single data frame or matrix returns a distance matrix of the row-wise dissimilarities.  \code{paldist2} takes two data frames of matrices and returns a matrix of all row-wise dissimilarities between the two datasets.
}

\value{
Function \code{MAT} returns an object of class \code{MAT} which contains the following items:
\item{call}{ original function call to \code{MAT}. }
\item{fitted.vales}{ fitted (predicted) values for the training set, as the mean and weighted mean (weighed by dissimilarity) of the k closest analogues. }
\item{diagnostics}{ standard deviation of the k analogues and dissimilarity of the closest analogue. }
\item{dist.n}{ dissimilarities of the k closest analogues. }
\item{x.n}{ environmental values of the k closest analogues. }
\item{match.name}{ column names of the k closest analogues. }
\item{x}{ environmental variable used in the model. }
\item{dist.method}{ dissimilarity coefficient. }
\item{k}{ number of closest analogues to use. }
\item{y}{ original species data. }
\item{cv.summary}{ summary of the cross-validation (not yet implemented). }
\item{dist}{ dissimilarity matrix (returned if \code{lean=FALSE}). }

If function \code{predict} is called with \code{newdata=NULL} it returns a matrix of fitted values from the original training set analysis.  If \code{newdata} is not \code{NULL} it returns list with the following named elements:
\item{fit}{ predictions for \code{newdata}. }
\item{diagnostics}{ standard deviations of the k closest analogues and distance of closest analogue. }
\item{dist.n}{ dissimilarities of the k closest analogues. }
\item{x.n}{ environmental values of the k closest analogues. }
\item{match.name}{ column names of the k closest analogues. }
\item{dist}{ dissimilarity matrix (returned if \code{lean=FALSE}). }
If sample specific errors were requested the list will also include:
\item{fit.boot}{ mean of the bootstrap estimates of newdata. }
\item{v1}{ squared standard error of the bootstrap estimates for each new sample. }
\item{v2}{ mean squared error for the training set samples, across all bootstram samples. }
\item{SEP}{ standard error of prediction, calculated as the square root of v1 + v2. }

Functions \code{paldist} and \code{paldist2} return dissimilarity matrices.  \code{performance} returns a matrix of performance statistics for the MAT model, with columns for RMSE, R2, mean and max bias for each number of analogues up to k. See \code{\link{performance}} for a description of the output.
}

\author{ 
Steve Juggins
}

\references{
Legendre, P. & Gallagher, E. (2001) Ecologically meaningful transformations for ordination of species. \emph{Oecologia}, \bold{129}, 271-280.

Overpeck, J.T., Webb, T., III, & Prentice, I.C. (1985) Quantitative interpretation of fossil pollen spectra: dissimilarity coefficients and the method of modern analogs. \emph{Quaternary Research}, \bold{23}, 87-108.
}

\seealso{ 
\code{\link{WAPLS}}, \code{\link{WA}}, \code{\link{performance}}, and \code{\link{compare.datasets}} for diagnostics.
}

\examples{
# pH reconstruction of the RLGH, Scotland, using SWAP training set 
# shows recent acidification history
data(SWAP)
data(RLGH)
fit <- MAT(SWAP$spec, SWAP$pH, k=20)  # generate results for k 1-20
#examine performance
performance(fit)
print(fit)
# How many analogues?
screeplot(fit)
# do the reconstruction
pred.mat <- predict(fit, RLGH$spec, k=10)
# plot the reconstruction
plot(RLGH$depths$Age, pred.mat$fit[, 1], type="b", ylab="pH", xlab="Age")

#compare to a weighted average model
fit <- WA(SWAP$spec, SWAP$pH)
pred.wa <- predict(fit, RLGH$spec)
points(RLGH$depths$Age, pred.wa$fit[, 1], col="red", type="b")
legend("topleft", c("MAT", "WA"), lty=1, col=c("black", "red"))
}
\keyword{ models }
\keyword{ multivariate }
