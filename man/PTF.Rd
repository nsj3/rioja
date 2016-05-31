\name{PTF}
\alias{performance}
\alias{performance.default}
\alias{crossval}
\alias{crossval.default}
\alias{rand.t.test}
\alias{rand.t.test.default}

\title{Palaeoecological transfer functions}
\description{
Functions for diagnosing and interpreting palaeoecological transfer functions.
}
\usage{
\method{performance}{default}(object, \dots)

\method{crossval}{default}(object, \dots)
}
\arguments{
  \item{object}{ a transfer function model from \code{wa}, \code{wapls} etc. }
  \item{\dots}{ additional arguments. }
}

\details{
Package \code{\link{rioja}} implements a number of numerical methods for inferring the value of an environmental variable from a set of sepecies abundances, given a modern training set of species data and  associated environmental values.  In palaeoecology these are known as "transfer functions" or "inference models" and are used to hindcast or "reconstruct" past environmental conditions from sub-fossil species assemblages preserved in sediment cores. The techniques included are weighted averaging (\code{\link{WA}}), partial least squares (PLS) and weighted average partial least squared (\code{\link{WAPLS}}), Imbrie and Kipp Factor Analysis (\code{\link{IKFA}}) a form of principal components regression, Maximum Likelihood Response Curves (\code{\link{MLRC}}), and the Modern Analogue Technique (\code{\link{MAT}}, a form of k-NN non-parametric regression (see Juggins & Birks (2010) for a review). 

The techniques are implemented in a consistent way and include functions for fitting a model to a training set of species and environmental data, with the function name named after the technique: that is, \code{WA} fits a weighted averaging model. Any model can be cross-validated using the \code{crossval} function, which allows internal cross-validation using leave-one-out, leave-n-out, bootstrapping or h-block cross-validation. There are a number of generic functions that can be used to summarise and diagnose the models: (\code{print}, \code{summary}, \code{performance} and \code{plot}. Some techniques have additional diagnostic functions such as \code{screeplot} and \code{rand.t.test} to help estimate the approproate number of components (WAPLS), factors (IKFA) or number of analogues (IKFA).

Predictions for new species data can be made using \code{predict}, with an option to calculate sample-specific errors using bootstrapping, after the method described in Birks et al. (1990).
}
\value{
Function \code{performance} returns a list with a named matrix \code{object} which contains the following columns:
\item{RMSE}{ root mean squared error, defined as the square root of the average sqaured error between the observed and predicted values for the training set.}
\item{R2}{ squared correlation betweenn observed and predicted values.}
\item{Avg.Bias}{ mean bias (mean of the residuals between measured and predicted values). } 
\item{Max.Bias}{ maximum bias, calculated by dividing the environmental gradient into a number of equal spaced segments (10 by default) and calculating the average bias for each segment.  The maximum bias is maximum of these 10 values and quantifies the tendendy for the model to over- or under-estimate at particular part of the gradient (ter Braak & Juggins 1993). }

If the transfer function object has been cross-validated, (ie. is the output of \code{crossval}, the list returned by \code{performance} also contains a matrix named crossval, which contains the above statistics calculated for the cross-validation predictions.

Function \code{crossval} returns an object of the original class and adds the following named elements:
\item{predicted}{ predicted values of each training set sample under cross-validation. }
\item{residuals.cv}{ prediction residuals. }

Function \code{rand.t.test} is a generic function that performs a randomisation t-test to test the significance of a cross-validated model, after van der Voet (1994).  Methods exist for \code{\link{WA}}, \code{\link{WAPLS}} and \code{\link{IKFA}}.
}

\author{ 
Steve Juggins
}

\references{
Birks, H.J.B., Line, J.M., Juggins, S., Stevenson, A.C., & ter Braak, C.J.F. (1990) Diatoms and pH reconstruction. \emph{Philosophical Transactions of the Royal Society of London}, \bold{B, 327}, 263-278.

Juggins, S., & Birks, HJB. (2010) Environmental Reconstructions. In Birks et al. (eds) \emph{Tracking Environmental Change using Lake Sediments: Data Handling and Statistical Techniques}., Kluwer Academic Publishers.

van der Voet, H. (1994) Comparing the predictive accuracy of models uing a simple randomization test. \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{25}, 313-323.
}

\keyword{ utilities }
