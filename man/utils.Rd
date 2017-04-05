\name{utils}
\alias{make.dummy}
\alias{dummy2factor}
\alias{dot}
\alias{dot.default}
\alias{dot.data.frame}
\alias{Hill.N2}
\alias{site.summ}
\alias{sp.summ}

\title{Utility functions. }
\description{
Utility functions to perform simple computations, transformations, formatting etc.
}
\usage{
make.dummy(fact)

dummy2factor(x)

Hill.N2(df, margin=2)

site.summ(y, max.cut=c(2, 5, 10, 20))

sp.summ(y, n.cut=c(5, 10, 20))

}
\arguments{
  \item{fact}{ a factor to convert to a matrix of dummy variables. }
  \item{x}{ a matrix or data frame of dummy variables to convert to a factor. }
  \item{df}{ a data frame of species abundance data. }
  \item{margin}{ margin to calculate over: 1 = by rows, 2 = by columns. }
  \item{y}{ data frame or matrix of species by sites data. }
  \item{n.cut}{ cut levels of abundance for species summary (see below). }
  \item{max.cut}{ cut levels of occurence for species sumamry. }
}

\details{
  Function \code{make.dummy} converts a factor into a matrix of dummy (1/0) variables. \code{dummy2factor} converts a matrix or data frame of dummy variables into a factor.
  
  Function \code{Hill.N2} returns Hill's N2 values for species or samples for a given species by sites dataset (Hill 1973).
}

\value{
  \code{make.dummy} returns a matrix of dummay variables. \code{dummy2factor} returns a factor.

  \code{Hill.N2} returns a numeric vector of N2 values.
  
  \code{sp.summ} returns a matrix with columns for the number of occurences, Hill's N2 and maximum abundance of each species, and the number of occurences at abundance greater than the cut levels given in \code{n.cut}.

  \code{sam.summ} returns a matrix with columns for the number of taxa, Hill's N2, maximum value and site total of each site (sample), and the number of taxa in each site with abundance greater than the cut levels given in \code{max.cut}.
}

\references{
Hill, M.O. (1973) Diversity and evenness: a unifying notation and its consequences. \emph{Ecology}, \bold{54}, 427-432.
}

\author{ 
Steve Juggins
}

\keyword{ utilities }
