##'Appending \code{Surv} objects
##'
##'\code{rbind} method for \code{Surv} objects
##'
##' 
##'@param ... \code{Surv} objects
##'@return \code{Surv} object
##'@author Klaus K. Holst
##'@keywords utilities
##'@examples
##'
##'y <- yl <- yr <- rnorm(10)
##'yl[1:5] <- NA; yr[6:10] <- NA
##'S1 <- survival::Surv(yl,yr,type="interval2")
##'S2 <- survival::Surv(y,y>0,type="right")
##'S3 <- survival::Surv(y,y<0,type="left")
##'
##'rbind(S1,S1)
##'rbind(S2,S2)
##'rbind(S3,S3)
##'
##' @export
rbind.Surv <- function(...) 
{
  dots <- list(...)
  type <- attributes(dots[[1]])$type
  ncol <- dim(dots[[1]])[2]
  nrow <- unlist(lapply(dots,nrow))
  cnrow <- c(0,cumsum(nrow))
  M <- matrix(ncol=ncol,nrow=sum(nrow))
  for (i in 1:length(dots)) {
    M[(cnrow[i]+1):cnrow[i+1],] <- dots[[i]]
  }
  x <- c(); for (i in 1:ncol(M)) x <- c(x,list(M[,i]))
  x <- c(x,list(type=type))
  do.call(survival::Surv, x)
} 
