% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.default.R
\name{sim.default}
\alias{sim.default}
\alias{as.sim}
\title{Monte Carlo simulation}
\usage{
\method{sim}{default}(
  x = NULL,
  R = 100,
  f = NULL,
  colnames = NULL,
  seed = NULL,
  args = list(),
  iter = FALSE,
  mc.cores,
  ...
)
}
\arguments{
\item{x}{function or 'sim' object}

\item{R}{Number of replications or data.frame with parameters}

\item{f}{Optional function (i.e., if x is a matrix)}

\item{colnames}{Optional column names}

\item{seed}{(optional) Seed (needed with cl=TRUE)}

\item{args}{(optional) list of named arguments passed to (mc)mapply}

\item{iter}{If TRUE the iteration number is passed as first argument to
(mc)mapply}

\item{mc.cores}{Optional number of cores. Will use parallel::mcmapply instead of future}

\item{...}{Additional arguments to future.apply::future_mapply}
}
\description{
Applies a function repeatedly for a specified number of replications or over
a list/data.frame with plot and summary methods for summarizing the Monte
Carlo experiment. Can be parallelized via the future package (use the
future::plan function).
}
\details{
To parallelize the calculation use the future::plan function (e.g.,
  future::plan(multisession()) to distribute the calculations over the R
  replications on all available cores). The output is controlled via the
  progressr package (e.g., progressr::handlers(global=TRUE) to enable
  progress information).
}
\examples{
m <- lvm(y~x+e)
distribution(m,~y) <- 0
distribution(m,~x) <- uniform.lvm(a=-1.1,b=1.1)
transform(m,e~x) <- function(x) (1*x^4)*rnorm(length(x),sd=1)

onerun <- function(iter=NULL,...,n=2e3,b0=1,idx=2) {
    d <- sim(m,n,p=c("y~x"=b0))
    l <- lm(y~x,d)
    res <- c(coef(summary(l))[idx,1:2],
             confint(l)[idx,],
             estimate(l,only.coef=TRUE)[idx,2:4])
    names(res) <- c("Estimate","Model.se","Model.lo","Model.hi",
                    "Sandwich.se","Sandwich.lo","Sandwich.hi")
    res
}
val <- sim(onerun,R=10,b0=1)
val

val <- sim(val,R=40,b0=1) ## append results
summary(val,estimate=c(1,1),confint=c(3,4,6,7),true=c(1,1))

summary(val,estimate=c(1,1),se=c(2,5),names=c("Model","Sandwich"))
summary(val,estimate=c(1,1),se=c(2,5),true=c(1,1),names=c("Model","Sandwich"),confint=TRUE)

if (interactive()) {
    plot(val,estimate=1,c(2,5),true=1,names=c("Model","Sandwich"),polygon=FALSE)
    plot(val,estimate=c(1,1),se=c(2,5),main=NULL,
         true=c(1,1),names=c("Model","Sandwich"),
         line.lwd=1,col=c("gray20","gray60"),
         rug=FALSE)
    plot(val,estimate=c(1,1),se=c(2,5),true=c(1,1),
         names=c("Model","Sandwich"))
}

f <- function(a=1, b=1) {
  rep(a*b, 5)
}
R <- Expand(a=1:3, b=1:3)
sim(f, R)
sim(function(a,b) f(a,b), 3, args=c(a=5,b=5))
sim(function(iter=1,a=5,b=5) iter*f(a,b), iter=TRUE, R=5)
}
\seealso{
summary.sim plot.sim print.sim
}
