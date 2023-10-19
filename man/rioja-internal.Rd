\name{rioja-internal}
\alias{.max.bias}
\alias{.r2}
\alias{.crossval}
\alias{.predict}
\alias{.check.data}
\alias{.print.crossval}
\alias{.performance}
\alias{.print.performance}
\alias{.set.rand.seed}
\alias{.get.rand}
\alias{predict.internal.WA}
\alias{predict.internal.WAPLS}
\alias{predict.internal.IKFA}
\alias{predict.internal.MLRC}
\alias{predict.internal.MR}
\alias{predict.internal.LWR}
\alias{.rmse}
\alias{.mono.fit}
\alias{.mono.predict}

\title{Internal rioja functions. }
\description{
Internal rioja functions.
}
\usage{
.max.bias(resi, x, ng = 10)
.rmse(x)
.r2(x, obs)
.mono.fit(x, y)
.mono.predict(mod, newdata)
}
\arguments{
  \item{resi}{ vector of residuals. }
  \item{x, y}{ vector of observed or predicted values. }
  \item{obs}{ vector of observed values. }
  \item{mod}{ model fitted by .mono.fit. }
  \item{newdata}{ data frame or matrix of new data to predict. }
}
\details{
  Internal \code{rioja} functions, not to be called by the user.
}

\author{ 
Steve Juggins
}

\keyword{ internal }
