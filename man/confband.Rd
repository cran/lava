\name{confband}
\alias{confband}
\title{Add Confidence limits bar to plot}
\usage{
  confband(x, lower, upper, center = NULL, delta = 0.07,
    centermark = 0.03, pch, blank = TRUE, vert = TRUE, ...)
}
\arguments{
  \item{x}{Position (x-coordinate if vert=TRUE,
  y-coordinate otherwise)}

  \item{lower}{Lower limit (if NULL no limits is added, and
  only the center is drawn (if not NULL))}

  \item{upper}{Upper limit}

  \item{center}{Center point}

  \item{delta}{Length of limit bars}

  \item{centermark}{Length of center bar}

  \item{pch}{Center symbol (if missing a line is drawn)}

  \item{blank}{If TRUE a white ball is plotted before the
  center is added to the plot}

  \item{vert}{If TRUE a vertical bar is plotted. Otherwise
  a horizontal bar is used}

  \item{\dots}{Additional low level arguments (e.g. col,
  lwd, lty,...)}
}
\description{
  Add Confidence limits bar to plot
}
\examples{
plot(0,0,type="n",xlab="",ylab="")
confband(0.5,-0.5,0.5,0,col="darkblue")
confband(0.8,-0.5,0.5,0,col="darkred",vert=FALSE,pch=1,cex=1.5)
}
\author{
  Klaus K. Holst
}
\seealso{
  \code{confband}
}
\keyword{iplot}
