##' @S3method simulate lvm
simulate.lvm <- function(x,n=100,p=NULL,sigma=1,rho=.5,
                         X,unlink=FALSE,...) {
  require("mvtnorm")
  if (!missing(X)) {
    n <- nrow(X)
  }

  nn <- setdiff(vars(x),parameter(x))
  mu <- unlist(lapply(x$mean, function(l) ifelse(is.na(l)|is.character(l),0,l)))
  xf <- intersect(unique(parlabels(x)),exogenous(x))
  xfix <- c(randomslope(x),xf); if (length(xfix)>0) normal <- FALSE
  ##  M <- modelVar(x,p,data=NULL)
  A <- x$M; P <- x$cov ##Sigma <- M$P
  mu <- unlist(x$mean)
  mu[is.na(mu)] <- 0
  E <- rmvnorm(n,rep(0,ncol(P)),as.matrix(P)) ## Error term
  
  ## Simulate exogenous variables (covariates)
  res <- matrix(0,ncol=length(nn),nrow=n); colnames(res) <- nn
  res <- as.data.frame(res)
  xx <- unique(c(exogenous(x, latent=TRUE, index=FALSE),xfix))
  X.idx <- match(xx,vars(x))
  res[,X.idx] <- t(mu[X.idx]+t(E[,X.idx]))
  if (missing(X)) {
    if (!is.null(xx) && length(xx)>0)
      for (i in 1:length(xx)) {
        mu.x <- mu[X.idx[i]]
        dist.x <- distribution(x,xx[i])[[1]]
        if (is.function(dist.x)) {
          res[,X.idx[i]] <- dist.x(n=n,mu=mu.x,var=P[X.idx[i],X.idx[i]])
        } else {
          if (is.null(dist.x) || is.na(dist.x)) {
            ##        res[,X.idx[i]] <- rnorm(n,mu.x,sd=Sigma[X.idx[i],X.idx[i]]^0.5)
            ##res[,X.idx[i]] <- mu.x+E[,X.idx[i]]
          } else {
            res[,X.idx[i]] <- dist.x ## Deterministic
          }
        }
      }
  } else {
    res[,X.idx] <- X[,xx]
  }
  simuled <- xx
  resunlink <- NULL
  if (unlink) {
    resunlink <- res
  }
  
  colnames(E) <- vars(x)
  ##  E <- heavytail.sim.hook(x,E)  
  while (length(simuled)<length(nn)) {
    leftovers <- setdiff(nn,simuled)
    
    for (i in leftovers) {
      browser()
      pos <- match(i,vars(x))
      relations <- colnames(A)[A[,pos]!=0]
      if (all(relations%in%simuled)) { ## Only depending on already simulated variable.namess
        mu.i <- mu[pos]        
        for (From in relations) {
          f <- functional(x,i,From)[[1]]
          if (!is.function(f))
            f <- function(x) x
          reglab <- regfix(x)$labels[From,pos]
          mu.i <- mu.i + A[From,pos]*f(res[,From])
        }
        dist.i <- distribution(x,i)[[1]]
        if (!is.function(dist.i)) {
          res[,pos] <- mu.i + E[,pos]
          if (unlink)
            resunlink[,pos] <- res[,pos]
        }
        ##          res[,pos] <- rnorm(n,mu.i,sd=Sigma[pos,pos]^0.5)
        else {
          res[,pos] <- dist.i(n=n,mu=mu.i,var=P[pos,pos])
          if (unlink)
            resunlink[,pos] <- mu.i
        }          
        simuled <- c(simuled,i)
      }
    }
  }
  res <- res[,nn,drop=FALSE]

  myhooks <- gethook("sim.hooks")
  for (f in myhooks) {
    res <- do.call(f, list(x=x,data=res))
  }         
  if (unlink) res <- resunlink
  return(data.frame(res))
}


