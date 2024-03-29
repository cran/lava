% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotConf.R
\name{plotConf}
\alias{plotConf}
\title{Plot regression lines}
\usage{
plotConf(
  model,
  var1 = NULL,
  var2 = NULL,
  data = NULL,
  ci.lty = 0,
  ci = TRUE,
  level = 0.95,
  pch = 16,
  lty = 1,
  lwd = 2,
  npoints = 100,
  xlim,
  col = NULL,
  colpt,
  alpha = 0.5,
  cex = 1,
  delta = 0.07,
  centermark = 0.03,
  jitter = 0.2,
  cidiff = FALSE,
  mean = TRUE,
  legend = ifelse(is.null(var1), FALSE, "topright"),
  trans = function(x) {
     x
 },
  partres = inherits(model, "lm"),
  partse = FALSE,
  labels,
  vcov,
  predictfun,
  plot = TRUE,
  new = TRUE,
  ...
)
}
\arguments{
\item{model}{Model object (e.g. \code{lm})}

\item{var1}{predictor (Continuous or factor)}

\item{var2}{Factor that interacts with \code{var1}}

\item{data}{data.frame to use for prediction (model.frame is used as default)}

\item{ci.lty}{Line type for confidence limits}

\item{ci}{Boolean indicating wether to draw pointwise 95\% confidence limits}

\item{level}{Level of confidence limits (default 95\%)}

\item{pch}{Point type for partial residuals}

\item{lty}{Line type for estimated regression lines}

\item{lwd}{Line width for regression lines}

\item{npoints}{Number of points used to plot curves}

\item{xlim}{Range of x axis}

\item{col}{Color (for each level in \code{var2})}

\item{colpt}{Color of partial residual points}

\item{alpha}{Alpha level}

\item{cex}{Point size}

\item{delta}{For categorical \code{var1}}

\item{centermark}{For categorical \code{var1}}

\item{jitter}{For categorical \code{var1}}

\item{cidiff}{For categorical \code{var1}}

\item{mean}{For categorical \code{var1}}

\item{legend}{Boolean (add legend)}

\item{trans}{Transform estimates (e.g. exponential)}

\item{partres}{Boolean indicating whether to plot partial residuals}

\item{partse}{.}

\item{labels}{Optional labels of \code{var2}}

\item{vcov}{Optional variance estimates}

\item{predictfun}{Optional predict-function used to calculate confidence limits and predictions}

\item{plot}{If FALSE return only predictions and confidence bands}

\item{new}{If FALSE add to current plot}

\item{\dots}{additional arguments to lower level functions}
}
\value{
list with following members:
\item{x}{Variable on the x-axis (\code{var1})}
\item{y}{Variable on the y-axis (partial residuals)}
\item{predict}{Matrix with confidence limits and predicted values}
}
\description{
Plot regression line (with interactions) and partial residuals.
}
\examples{
n <- 100
x0 <- rnorm(n)
x1 <- seq(-3,3, length.out=n)
x2 <- factor(rep(c(1,2),each=n/2), labels=c("A","B"))
y <- 5 + 2*x0 + 0.5*x1 + -1*(x2=="B")*x1 + 0.5*(x2=="B") + rnorm(n, sd=0.25)
dd <- data.frame(y=y, x1=x1, x2=x2)
lm0 <- lm(y ~ x0 + x1*x2, dd)
plotConf(lm0, var1="x1", var2="x2")
abline(a=5,b=0.5,col="red")
abline(a=5.5,b=-0.5,col="red")
### points(5+0.5*x1 -1*(x2=="B")*x1 + 0.5*(x2=="B") ~ x1, cex=2)

data(iris)
l <- lm(Sepal.Length ~ Sepal.Width*Species,iris)
plotConf(l,var2="Species")
plotConf(l,var1="Sepal.Width",var2="Species")

\dontrun{
## lme4 model
dd$Id <- rbinom(n, size = 3, prob = 0.3)
lmer0 <- lme4::lmer(y ~ x0 + x1*x2 + (1|Id), dd)
plotConf(lmer0, var1="x1", var2="x2")
}
}
\seealso{
\code{termplot}
}
\author{
Klaus K. Holst
}
\keyword{hplot}
\keyword{regression}
