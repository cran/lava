% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{Model}
\alias{Model}
\alias{Model<-}
\title{Extract model}
\usage{
Model(x, ...)

Model(x, ...) <- value
}
\arguments{
\item{x}{Fitted model}

\item{\dots}{Additional arguments to be passed to the low level functions}

\item{value}{New model object (e.g. \code{lvm} or \code{multigroup})}
}
\value{
Returns a model object (e.g. \code{lvm} or \code{multigroup})
}
\description{
Extract or replace model object
}
\examples{

m <- lvm(y~x)
e <- estimate(m, sim(m,100))
Model(e)

}
\seealso{
\code{\link{Graph}}
}
\author{
Klaus K. Holst
}
\keyword{models}
