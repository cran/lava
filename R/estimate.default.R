##' @export
estimate <- function(x,...) UseMethod("estimate")

##' @export
estimate.list <- function(x,...) {
    if (inherits(x[[1]],"lvm")) return(estimate.lvmlist(x,...))
    lapply(x,function(x) estimate(x,...))
}

##' Estimation of functional of parameters
##'
##' Estimation of functional of parameters.
##' Wald tests, robust standard errors, cluster robust standard errors,
##' LRT (when \code{f} is not a function)...
##' @param x model object (\code{glm}, \code{lvmfit}, ...)
##' @param f transformation of model parameters and (optionally) data, or contrast matrix (or vector)
##' @param ... additional arguments to lower level functions
##' @param data \code{data.frame}
##' @param id (optional) id-variable corresponding to iid decomposition of model parameters.
##' @param iddata (optional) id-variable for 'data'
##' @param stack if TRUE (default)  the i.i.d. decomposition is automatically stacked according to 'id'
##' @param average if TRUE averages are calculated
##' @param subset (optional) subset of data.frame on which to condition (logical expression or variable name)
##' @param score.deriv (optional) derivative of mean score function
##' @param level level of confidence limits
##' @param iid if TRUE (default) the iid decompositions are also returned (extract with \code{iid} method)
##' @param type type of small-sample correction
##' @param keep (optional) index of parameters to keep from final result
##' @param use (optional) index of parameters to use in calculations
##' @param regex If TRUE use regular expression (perl compatible) for keep,use arguments
##' @param contrast (optional) Contrast matrix for final Wald test
##' @param null (optional) null hypothesis to test
##' @param vcov (optional) covariance matrix of parameter estimates (e.g. Wald-test)
##' @param coef (optional) parameter coefficient
##' @param robust if TRUE robust standard errors are calculated. If
##' FALSE p-values for linear models are calculated from t-distribution
##' @param df degrees of freedom (default obtained from 'df.residual')
##' @param print (optional) print function
##' @param labels (optional) names of coefficients
##' @param label.width (optional) max width of labels
##' @param only.coef if TRUE only the coefficient matrix is return
##' @param back.transform (optional) transform of parameters and confidence intervals
##' @param folds (optional) aggregate influence functions (divide and conquer)
##' @param cluster (obsolete) alias for 'id'.
##' @param R Number of simulations (simulated p-values)
##' @param null.sim Mean under the null for simulations
##' @details
##'
##' iid decomposition
##' \deqn{\sqrt{n}(\widehat{\theta}-\theta) = \sum_{i=1}^n\epsilon_i + o_p(1)}
##' can be extracted with the \code{iid} method.
##'
##' @export
##' @examples
##'
##' ## Simulation from logistic regression model
##' m <- lvm(y~x+z);
##' distribution(m,y~x) <- binomial.lvm("logit")
##' d <- sim(m,1000)
##' g <- glm(y~z+x,data=d,family=binomial())
##' g0 <- glm(y~1,data=d,family=binomial())
##'
##' ## LRT
##' estimate(g,g0)
##'
##' ## Plain estimates (robust standard errors)
##' estimate(g)
##'
##' ## Testing contrasts
##' estimate(g,null=0)
##' estimate(g,rbind(c(1,1,0),c(1,0,2)))
##' estimate(g,rbind(c(1,1,0),c(1,0,2)),null=c(1,2))
##' estimate(g,2:3) ## same as cbind(0,1,-1)
##' estimate(g,as.list(2:3)) ## same as rbind(c(0,1,0),c(0,0,1))
##' ## Alternative syntax
##' estimate(g,"z","z"-"x",2*"z"-3*"x")
##' estimate(g,z,z-x,2*z-3*x)
##' estimate(g,"?")  ## Wilcards
##' estimate(g,"*Int*","z")
##' estimate(g,"1","2"-"3",null=c(0,1))
##' estimate(g,2,3)
##'
##' ## Usual (non-robust) confidence intervals
##' estimate(g,robust=FALSE)
##'
##' ## Transformations
##' estimate(g,function(p) p[1]+p[2])
##'
##' ## Multiple parameters
##' e <- estimate(g,function(p) c(p[1]+p[2],p[1]*p[2]))
##' e
##' vcov(e)
##'
##' ## Label new parameters
##' estimate(g,function(p) list("a1"=p[1]+p[2],"b1"=p[1]*p[2]))
##' ##'
##' ## Multiple group
##' m <- lvm(y~x)
##' m <- baptize(m)
##' d2 <- d1 <- sim(m,50)
##' e <- estimate(list(m,m),list(d1,d2))
##' estimate(e) ## Wrong
##' estimate(e,id=rep(seq(nrow(d1)),2))
##' estimate(lm(y~x,d1))
##'
##' ## Marginalize
##' f <- function(p,data)
##'   list(p0=lava:::expit(p["(Intercept)"] + p["z"]*data[,"z"]),
##'        p1=lava:::expit(p["(Intercept)"] + p["x"] + p["z"]*data[,"z"]))
##' e <- estimate(g, f, average=TRUE)
##' e
##' estimate(e,diff)
##' estimate(e,cbind(1,1))
##'
##' ## Clusters and subset (conditional marginal effects)
##' d$id <- rep(seq(nrow(d)/4),each=4)
##' estimate(g,function(p,data)
##'          list(p0=lava:::expit(p[1] + p["z"]*data[,"z"])),
##'          subset=d$z>0, id=d$id, average=TRUE)
##'
##' ## More examples with clusters:
##' m <- lvm(c(y1,y2,y3)~u+x)
##' d <- sim(m,10)
##' l1 <- glm(y1~x,data=d)
##' l2 <- glm(y2~x,data=d)
##' l3 <- glm(y3~x,data=d)
##'
##' ## Some random id-numbers
##' id1 <- c(1,1,4,1,3,1,2,3,4,5)
##' id2 <- c(1,2,3,4,5,6,7,8,1,1)
##' id3 <- seq(10)
##'
##' ## Un-stacked and stacked i.i.d. decomposition
##' iid(estimate(l1,id=id1,stack=FALSE))
##' iid(estimate(l1,id=id1))
##'
##' ## Combined i.i.d. decomposition
##' e1 <- estimate(l1,id=id1)
##' e2 <- estimate(l2,id=id2)
##' e3 <- estimate(l3,id=id3)
##' (a2 <- merge(e1,e2,e3))
##'
##' ## If all models were estimated on the same data we could use the
##' ## syntax:
##' ## Reduce(merge,estimate(list(l1,l2,l3)))
##'
##' ## Same:
##' iid(a1 <- merge(l1,l2,l3,id=list(id1,id2,id3)))
##'
##' iid(merge(l1,l2,l3,id=TRUE)) # one-to-one (same clusters)
##' iid(merge(l1,l2,l3,id=FALSE)) # independence
##'
##'
##' ## Monte Carlo approach, simple trend test example
##'
##' m <- categorical(lvm(),~x,K=5)
##' regression(m,additive=TRUE) <- y~x
##' d <- simulate(m,100,seed=1,'y~x'=0.1)
##' l <- lm(y~-1+factor(x),data=d)
##'
##' f <- function(x) coef(lm(x~seq_along(x)))[2]
##' null <- rep(mean(coef(l)),length(coef(l))) ## just need to make sure we simulate under H0: slope=0
##' estimate(l,f,R=1e2,null.sim=null)
##'
##' estimate(l,f)
##' @aliases estimate estimate.default estimate.estimate merge.estimate
##' @method estimate default
##' @export
estimate.default <- function(x=NULL,f=NULL,...,data,id,
                      iddata,stack=TRUE,average=FALSE,subset,
                      score.deriv,level=0.95,iid=robust,
                      type=c("robust","df","mbn"),
                      keep,use,
                      regex=FALSE,
                      contrast,null,vcov,coef,
                      robust=TRUE,df=NULL,
                      print=NULL,labels,label.width,
                      only.coef=FALSE,back.transform=NULL,
                      folds=0,
                      cluster,
                      R=0,
                      null.sim) {
    cl <- match.call(expand.dots=TRUE)
    if (!missing(use)) {
        p0 <- c("f","contrast","only.coef","subset","average","keep","labels")
        cl0 <- cl
        cl0[c("use",p0)] <- NULL
        cl0$keep <- use
        cl$x <- eval(cl0,parent.frame())
        cl[c("vcov","use")] <- NULL
        return(eval(cl,parent.frame()))
    }
    expr <- suppressWarnings(inherits(try(f,silent=TRUE),"try-error"))
    if (!missing(coef)) {
        pp <- coef
    } else {
        pp <- suppressWarnings(try(stats::coef(x),"try-error"))
        if (inherits(x,"survreg") && length(pp)<NROW(x$var)) {
            pp <- c(pp,scale=x$scale)
        }
    }
    if (!missing(cluster)) id <- cluster
    if (expr || is.character(f) || (is.numeric(f) && !is.matrix(f))) { ## || is.call(f)) {
        dots <- lapply(substitute(placeholder(...))[-1],function(x) x)
        args <- c(list(coef=names(pp),x=substitute(f),regex=regex),dots)
        f <- do.call(parsedesign,args)
    }
    if (!is.null(f) && !is.function(f)) {
        if (!(is.matrix(f) | is.vector(f))) return(compare(x,f,...))
        contrast <- f; f <- NULL
    }

    if (lava.options()$cluster.index) {
        if (!requireNamespace("mets",quietly=TRUE)) stop("'mets' package required")
    }
    if (missing(data)) data <- tryCatch(model.frame(x),error=function(...) NULL)
    ##if (is.matrix(x) || is.vector(x)) contrast <- x
    alpha <- 1-level
    alpha.str <- paste(c(alpha/2,1-alpha/2)*100,"",sep="%")
    nn <- NULL
    if ((( (is.logical(iid) && iid) || length(iid)>0) && robust) && (missing(vcov) || is.null(vcov) || (is.logical(vcov) && vcov[1]==FALSE && !is.na(vcov[1])))) { ## If user supplied vcov, then don't estimate IC
        if (missing(score.deriv)) {
            if (!is.logical(iid)) {
                iidtheta <- iid
                iid <- TRUE
            } else {
                suppressWarnings(iidtheta <- iid(x,folds=folds))
            }
        } else {
            suppressWarnings(iidtheta <- iid(x,score.deriv=score.deriv,folds=folds))
        }
    } else {
        if (missing(vcov) || (is.logical(vcov) && !is.na(vcov)[1])) suppressWarnings(vcov <- stats::vcov(x))
        iidtheta <- NULL
    }

    if (!missing(subset)) {
        e <- substitute(subset)
        expr <- suppressWarnings(inherits(try(subset,silent=TRUE),"try-error"))
        if (expr) subset <- eval(e,envir=data)
        ##subset <- eval(e, data, parent.frame())
        if (is.character(subset)) subset <- data[,subset]
        if (is.numeric(subset)) subset <- subset>0
    }
    idstack <- NULL
    ## Preserve id from 'estimate' object
    if (missing(id) && inherits(x,"estimate") && !is.null(x$id)) id <- x$id
    if (!missing(id) && iid) {
        if (is.null(iidtheta)) stop("'iid' method needed")
        nprev <- nrow(iidtheta)
        if (inherits(id,"formula")) {
            id <- interaction(get_all_vars(id,data))
        }
        ## e <- substitute(id)
        ## expr <- suppressWarnings(inherits(try(id,silent=TRUE),"try-error"))
        ## if (expr) id <- eval(e,envir=data)
        ##if (!is.null(data)) id <- eval(e, data)
        if (is.logical(id) && length(id)==1) {
            id <- if(is.null(iidtheta)) seq(nrow(data)) else seq(nprev)
            stack <- FALSE
        }
        if (is.character(id) && length(id)==1) id <- data[,id,drop=TRUE]
        if (!is.null(iidtheta)) {
            if (length(id)!=nprev) {
                if (!is.null(x$na.action) && (length(id)==length(x$na.action)+nprev)) {
                    warning("Applying na.action")
                    id <- id[-x$na.action]
                } else stop("Dimensions of i.i.d decomposition and 'id' does not agree")
            }
        } else {
            if (length(id)!=nrow(data)) {
                if (!is.null(x$na.action) && (length(id)==length(x$na.action)+nrow(data))) {
                    warning("Applying na.action")
                    id <- id[-x$na.action]
                } else stop("Dimensions of i.i.d decomposition and 'id' does not agree")
            }
        }
        if (stack) {
            N <- nrow(iidtheta)
            clidx <- NULL
            atr <- attributes(iidtheta)
            atr$dimnames <- NULL
            atr$dim <- NULL
            if (!lava.options()$cluster.index) {
                iidtheta <- matrix(unlist(by(iidtheta,id,colSums)),byrow=TRUE,ncol=ncol(iidtheta))
                attributes(iidtheta)[names(atr)] <- atr
                idstack <- sort(unique(id))
            } else {
                clidx <- mets::cluster.index(id,mat=iidtheta,return.all=TRUE)
                iidtheta <- clidx$X
                attributes(iidtheta)[names(atr)] <- atr
                idstack <- id[as.vector(clidx$firstclustid)+1]
            }
            if (is.null(attributes(iidtheta)$N)) {
                attributes(iidtheta)$N <- N
            }
        } else idstack <- id
    } else {
        if (!is.null(data)) idstack <- rownames(data)
    }

    if (!is.null(iidtheta) && (length(idstack)==nrow(iidtheta))) rownames(iidtheta) <- idstack
    if (!robust) {
        if (inherits(x,"lm") && family(x)$family=="gaussian" && is.null(df)) df <- x$df.residual
        if (missing(vcov)) suppressWarnings(vcov <- stats::vcov(x))
    }
    if (!is.null(iidtheta) && robust && (missing(vcov) || is.null(vcov))) {
        ## if (is.null(f))
        V <- crossprod(iidtheta)
        ### Small-sample corrections for clustered data
        K <- NROW(iidtheta)
        N <- attributes(iidtheta)$N
        if (is.null(N)) N <- K
        p <- NCOL(iidtheta)
        adj0 <- K/(K-p) ## Mancl & DeRouen, 2001
        adj1 <- K/(K-1) ## Mancl & DeRouen, 2001
        adj2 <- (N-1)/(N-p)*(K/(K-1)) ## Morel,Bokossa & Neerchal, 2003
        if (tolower(type[1])=="mbn" && !is.null(attributes(iidtheta)$bread)) {
            V0 <- V
            iI0 <- attributes(iidtheta)$bread
            I0 <- Inverse(iI0)
            ##I1 <- crossprod(iidtheta%*%I0)
            delta <- min(0.5,p/(K-p))
            phi <- max(1,tr(I0%*%V0)*adj2/p)
            V <- adj2*V0 + delta*phi*iI0
        }
        if (tolower(type[1])=="df") {
            V <- adj0*V
        }
        if (tolower(type[1])=="df1") {
            V <- adj1*V
        }
        if (tolower(type[1])=="df2") {
            V <- adj2*V
        }
    } else {
        if (!missing(vcov)) {
            if (length(vcov)==1 && is.na(vcov)) vcov <- matrix(NA,length(pp),length(pp))
            V <- vcov
        } else {
            suppressWarnings(V <- stats::vcov(x))
        }
    }

    ## Simulate p-value
    if (R>0) {
        if (is.null(f)) stop("Supply function 'f'")
        if (missing(null.sim)) null.sim <- rep(0,length(pp))
        est <- f(pp)
        if (is.list(est)) {
            nn <- names(est)
            est <- unlist(est)
            names(est) <- nn
        }
        if (missing(labels)) {
            labels <- colnames(rbind(est))
        }
        res <- simnull(R,f,mu=null.sim,sigma=V,labels=labels)
        return(structure(res, class=c("estimate.sim","sim"),
                         coef=pp,
                         vcov=V,
                         f=f,
                         estimate=est))
    }


    if (!is.null(f)) {
        form <- names(formals(f))
        dots <- ("..."%in%names(form))
        form0 <- setdiff(form,"...")
        parname <- "p"

        if (!is.null(form)) parname <- form[1] # unless .Primitive
        if (length(form0)==1 && !(form0%in%c("object","data"))) {
            ##names(formals(f))[1] <- "p"
            parname <- form0
        }
        if (!is.null(iidtheta)) {
            arglist <- c(list(object=x,data=data,p=vec(pp)),list(...))
            names(arglist)[3] <- parname
        } else {
            arglist <- c(list(object=x,p=vec(pp)),list(...))
            names(arglist)[2] <- parname
        }
        if (!dots) {
            arglist <- arglist[intersect(form0,names(arglist))]
        }
        newf <- NULL
        if (length(form)==0) {
            arglist <- list(vec(pp))
            ##newf <- function(p,...) do.call("f",list(p,...))
            newf <- function(...) do.call("f",list(...))
            val <- do.call("f",arglist)
        } else {
            val <- do.call("f",arglist)
            if (is.list(val)) {
                nn <- names(val)
                val <- do.call("cbind",val)
                ##newf <- function(p,...) do.call("cbind",f(p,...))
                newf <- function(...) do.call("cbind",f(...))
            }
        }
        k <- NCOL(val)
        N <- NROW(val)
        D <- attributes(val)$grad
        if (is.null(D)) {
            D <- numDeriv::jacobian(function(p,...) {
                if (length(form)==0) arglist[[1]] <- p
                else arglist[[parname]] <- p
                if (is.null(newf))
                    return(do.call("f",arglist))
                return(do.call("newf",arglist)) }, pp)
        }
        if (is.null(iidtheta)) {
            pp <- structure(as.vector(val),names=names(val))
            V <- D%*%V%*%t(D)
        } else {
            if (!average || (N<NROW(data))) { ## || NROW(data)==0)) { ## transformation not depending on data
                pp <- structure(as.vector(val),names=names(val))
                iidtheta <- iidtheta%*%t(D)
                ##V <- crossprod(iidtheta)
                V <- D%*%V%*%t(D)
            } else {
                if (k>1) { ## More than one parameter (and depends on data)
                    if (!missing(subset)) { ## Conditional estimate
                        val <- apply(val,2,function(x) x*subset)
                    }
                    D0 <- matrix(nrow=k,ncol=length(pp))
                    for (i in seq_len(k)) {
                        D1 <- D[seq(N)+(i-1)*N,,drop=FALSE]
                        if (!missing(subset)) ## Conditional estimate
                            D1 <- apply(D1,2,function(x) x*subset)
                        D0[i,] <- colMeans(D1)
                    }
                    D <- D0
                    iid2 <- iidtheta%*%t(D)
                } else { ## Single parameter
                    if (!missing(subset)) { ## Conditional estimate
                        val <- val*subset
                        D <- apply(rbind(D),2,function(x) x*subset)
                    }
                    D <- colMeans(rbind(D))
                    iid2 <- iidtheta%*%D
                }
                pp <- vec(colMeans(cbind(val)))
                iid1 <- (cbind(val)-rbind(pp)%x%cbind(rep(1,N)))/N
                if (!missing(id)) {
                    if (!lava.options()$cluster.index)
                        iid1 <- matrix(unlist(by(iid1,id,colSums)),byrow=TRUE,ncol=ncol(iid1))
                    else {
                        iid1 <- mets::cluster.index(id,mat=iid1,return.all=FALSE)
                    }
                }
                if (!missing(subset)) { ## Conditional estimate
                    phat <- mean(subset)
                    iid3 <- cbind(-1/phat^2 * (subset-phat)/N) ## check
                    if (!missing(id)) {
                        if (!lava.options()$cluster.index) {
                            iid3 <- matrix(unlist(by(iid3,id,colSums)),byrow=TRUE,ncol=ncol(iid3))
                        } else {
                            iid3 <- mets::cluster.index(id,mat=iid3,return.all=FALSE)
                        }
                    }
                    iidtheta <- (iid1+iid2)/phat + rbind(pp)%x%iid3
                    pp <- pp/phat
                    V <- crossprod(iidtheta)
                } else {
                    if (nrow(iid1)!=nrow(iid2)) {
                        message("Assuming independence between model iid decomposition and new data frame")
                        V <- crossprod(iid1) + crossprod(iid2)
                    } else {
                        iidtheta <- iid1+iid2
                        V <- crossprod(iidtheta)
                    }
                }
            }
        }
    }

    if (is.null(V)) {
        res <- cbind(pp,NA,NA,NA,NA)
    } else {
        if (length(pp)==1) res <- rbind(c(pp,diag(V)^0.5)) else res <- cbind(pp,diag(V)^0.5)
        beta0 <- res[,1]

        if (!missing(null) && missing(contrast)) beta0 <- beta0-null
        if (!is.null(df)) {
            za <- qt(1-alpha/2,df=df)
            pval <- 2*pt(abs(res[,1]/res[,2]),df=df,lower.tail=FALSE)
        } else {
            za <- qnorm(1-alpha/2)
            pval <- 2*pnorm(abs(res[,1]/res[,2]),lower.tail=FALSE)
        }
        res <- cbind(res,res[,1]-za*res[,2],res[,1]+za*res[,2],pval)
    }
    colnames(res) <- c("Estimate","Std.Err",alpha.str,"P-value")

    if (nrow(res)>0)
        if (!is.null(nn)) {
            rownames(res) <- nn
        } else {
            nn <- attributes(res)$varnames
            if (!is.null(nn)) rownames(res) <- nn
            if (is.null(rownames(res))) rownames(res) <- paste0("p",seq(nrow(res)))
        }


    if (NROW(res)==0L) {
        coefs <- NULL
    } else {
        coefs <- res[,1,drop=TRUE]; names(coefs) <- rownames(res)
    }
    res <- structure(list(coef=coefs,coefmat=res,vcov=V, iid=NULL, print=print, id=idstack),class="estimate")
    if (iid) ## && is.null(back.transform))
        res$iid <- iidtheta
    if (length(coefs)==0L) return(res)

    if (!missing(contrast) | !missing(null)) {
        p <- length(res$coef)
        if (missing(contrast)) contrast <- diag(nrow=p)
        if (missing(null)) null <- 0
        if (is.vector(contrast) || is.list(contrast)) {
            contrast <- contr(contrast, names(res$coef), ...)
            ## if (length(contrast)==p) contrast <- rbind(contrast)
            ## else {
            ##     cont <- contrast
            ##     contrast <- diag(nrow=p)[cont,,drop=FALSE]
            ## }
        }
        cc <- compare(res,contrast=contrast,null=null,vcov=V,level=level,df=df)
        res <- structure(c(res, list(compare=cc)),class="estimate")
        if (!is.null(df)) {
            pval <- with(cc,pt(abs(estimate[,1]-null)/estimate[,2],df=df,lower.tail=FALSE)*2)
        } else {
            pval <- with(cc,pnorm(abs(estimate[,1]-null)/estimate[,2],lower.tail=FALSE)*2)
        }
        res$coefmat <- with(cc, cbind(estimate,pval))
        colnames(res$coefmat)[5] <- "P-value"
        rownames(res$coefmat) <- cc$cnames
        if (!is.null(res$iid)) {
            res$iid <- res$iid%*%t(contrast)
            colnames(res$iid) <- cc$cnames
        }
        res$compare$estimate <- NULL
        res$coef <- res$compare$coef
        res$vcov <- res$compare$vcov
        names(res$coef) <- gsub("(^\\[)|(\\]$)","",rownames(res$coefmat))
    }

    if (!is.null(back.transform)) {
        res$coefmat[,c(1,3,4)] <- do.call(back.transform,list(res$coefmat[,c(1,3,4)]))
        res$coefmat[,2] <- NA
    }

    if (!missing(keep) && !is.null(keep)) {
        if (is.character(keep)) {
            if (regex) {
                nn <- rownames(res$coefmat)
                keep <- unlist(lapply(keep, function(x) grep(x,nn, perl=TRUE)))
            } else {
                keep <- match(keep,rownames(res$coefmat))
            }
        }
        res$coef <- res$coef[keep]
        res$coefmat <- res$coefmat[keep,,drop=FALSE]
        if (!is.null(res$iid)) res$iid <- res$iid[,keep,drop=FALSE]
        res$vcov <- res$vcov[keep,keep,drop=FALSE]
    }
    if (!missing(labels)) {
        names(res$coef) <- labels
        if (!is.null(res$iid)) colnames(res$iid) <- labels
        colnames(res$vcov) <- rownames(res$vcov) <- labels
        rownames(res$coefmat) <- labels
    }
    if (!missing(label.width)) {
        rownames(res$coefmat) <- make.unique(unlist(lapply(rownames(res$coefmat),
                                                           function(x) toString(x,width=label.width))))
    }
    if (only.coef) return(res$coefmat)
    res$call <- cl
    res$back.transform <- back.transform
    res$n <- nrow(data)
    res$ncluster <- nrow(res$iid)
    return(res)
}

simnull <- function(R,f,mu,sigma,labels=NULL) {
    X <- rmvn0(R,mu=mu,sigma=sigma)
    est <- f(mu)
    res <- apply(X,1,f)
    if (is.list(est)) {
        nn <- names(est)
        est <- unlist(est)
        names(est) <- nn
        res <- matrix(unlist(res),byrow=TRUE,ncol=length(est))
    } else {
        res <- t(rbind(res))
    }
    if (is.null(labels)) {
        labels <- colnames(rbind(est))
        if (is.null(labels)) labels <- paste0("p",seq_along(est))
    }
    colnames(res) <- labels
    return(res)
}

##' @export
estimate.estimate.sim <- function(x,f,R=0,labels,...) {
    atr <- attributes(x)
    if (R>0) {
        if (missing(f)) {
            val <- simnull(R,f=atr[["f"]],mu=atr[["coef"]],sigma=atr[["vcov"]])
            res <- rbind(x,val)
            for (a in setdiff(names(atr),c("dim","dimnames")))
                attr(res,a) <- atr[[a]]
        } else {
            res <- simnull(R,f=f,mu=atr[["coef"]],sigma=atr[["vcov"]])
            for (a in setdiff(names(atr),c("dim","dimnames","f")))
                attr(res,a) <- atr[[a]]
            attr(f,"f") <- f
            est <- unlist(f(atr[["coef"]]))
            if (missing(labels)) labels <- colnames(rbind(est))
            attr(res,"estimate") <- est
        }
        if (!missing(labels)) colnames(res) <- labels
        return(res)
    }
    if (missing(f)) {
        if (!missing(labels)) colnames(res) <- labels
        return(x)
    }

    est <- f(atr[["coef"]])
    res <- apply(x,1,f)
    if (is.list(est)) {
        res <- matrix(unlist(res),byrow=TRUE,ncol=length(est))
    } else {
        res <- t(rbind(res))
    }
    if (missing(labels)) {
        labels <- colnames(rbind(est))
        if (is.null(labels)) labels <- paste0("p",seq_along(est))
    }
    colnames(res) <- labels
    for (a in setdiff(names(atr),c("dim","dimnames","f","estimate")))
        attr(res,a) <- atr[[a]]
    attr(f,"f") <- f
    attr(res,"estimate") <- unlist(est)
    return(res)
}


##' @export
print.estimate.sim <- function(x,level=.95,...) {
    quantiles <- c((1-level)/2,1-(1-level)/2)
    est <- attr(x,"estimate")
    mysummary <- function(x,INDEX,...) {
        x <- as.vector(x)
        res <- c(mean(x,na.rm=TRUE),
                 sd(x,na.rm=TRUE),
                 quantile(x,quantiles,na.rm=TRUE),
                 est[INDEX],
                 mean(abs(x)>abs(est[INDEX]),na.rm=TRUE))

        names(res) <- c("Mean","SD",paste0(quantiles*100,"%"),
                        "Estimate","P-value")
        res
    }
    env <- new.env()
    assign("est",attr(x,"estimate"),env)
    environment(mysummary) <- env
    print(summary(x,fun=mysummary,...))
}

estimate.glm <- function(x,...) {
    estimate.default(x,...)
}

##' @export
print.estimate <- function(x, type=0L, digits=4L, width=25L,
                    std.error=TRUE, p.value=TRUE,
                    sep="_______",sep.which, sep.labels=NULL,indent=" ",
                    unique.names=TRUE,
                    na.print="", ...) {

    if (!is.null(x$print)) {
        x$print(x,digits=digits,width=width,...)
        return(invisible(x))
    }
    if (type>0 && !is.null(x$call)) {
        cat("Call: "); print(x$call)
        printline(50)
    }
    if (type>0) {
        if (!is.null(x[["n"]]) && !is.null(x[["k"]])) {
            cat("n = ",x[["n"]],", clusters = ",x[["k"]],"\n\n",sep="")
        } else {
            if (!is.null(x[["n"]])) {
                cat("n = ",x[["n"]],"\n\n",sep="")
            }
            if (!is.null(x[["k"]])) {
                cat("n = ",x[["k"]],"\n\n",sep="")
            }
        }
    }

    cc <- x$coefmat
    if (!is.null(rownames(cc)) && unique.names)
        rownames(cc) <- make.unique(unlist(lapply(rownames(cc),
                                                  function(x) toString(x,width=width))))
    if (!std.error) cc <- cc[,-2,drop=FALSE]
    if (!p.value) cc[,-ncol(cc),drop=FALSE]

    sep.pos <- c()
    if (missing(sep.which) && !is.null(x$model.index)) {
        sep.which <- unlist(lapply(x$model.index,function(x) tail(x,1)))[-length(x$model.index)]
    }
    if (missing(sep.which)) sep.which <- NULL

    
    if (!is.null(sep.which)) {
        sep0 <- 0%in%sep.which
        if (sep0) sep.which <- setdiff(sep.which,0)
        cc0 <- c()
        sep.which <- c(0,sep.which,nrow(cc))
        N <- length(sep.which)-1
        for (i in seq(N)) {
            if ((sep.which[i]+1)<=nrow(cc))
                cc0 <- rbind(cc0, cc[seq(sep.which[i]+1,sep.which[i+1]),,drop=FALSE])
            if (i<N) {
                cc0 <- rbind(cc0, NA)
                sep.pos <- c(sep.pos,nrow(cc0))
            }
        }
        if (sep0) {
            sep.pos <- c(1,sep.pos+1)
            cc0 <- rbind(NA, cc0)
        }
        cc <- cc0
    }
    if (!is.null(sep.labels)) {
        sep.labels <- rep(sep.labels, length.out=length(sep.pos))
        rownames(cc)[sep.pos] <- sep.labels
        rownames(cc)[-sep.pos] <- paste0(indent, rownames(cc)[-sep.pos])
    } else {
        if (length(sep.pos)>0) rownames(cc)[sep.pos] <- rep(paste0(rep("_",max(nchar(rownames(cc)))),collapse=""),length(sep.pos))
    }
    print(cc,digits=digits,na.print=na.print,...)

    if (!is.null(x$compare)) {
        cat("\n",x$compare$method[3],"\n")
        cat(paste(" ",x$compare$method[-(1:3)],collapse="\n"),"\n")
        if (length(x$compare$method)>4) {
            out <- character()
            out <- with(x$compare, c(out, paste(names(statistic), "=", format(round(statistic, 4)))))
            out <- with(x$compare, c(out, paste(names(parameter), "=", format(round(parameter,3)))))
            fp  <- with(x$compare, format.pval(p.value, digits = digits))
            out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<") fp else paste("=", fp)))
            cat(" ",strwrap(paste(out, collapse = ", ")), sep = "\n")
        }
    }
}

##' @export
vcov.estimate <- function(object,list=FALSE,...) {
    res <- object$vcov
    nn <- names(coef(object,...))
    if (list && !is.null(object$model.index)) {
        return(lapply(object$model.index, function(x) object$vcov[x,x]))
    }
    dimnames(res) <- list(nn,nn)
    res
}

##' @export
coef.estimate <- function(object,mat=FALSE,list=FALSE,messages=lava.options()$messages,...) {
    if (mat) return(object$coefmat)
    if (messages>0 && !is.null(object$back.transform)) message("Note: estimates on original scale (before 'back.transform')")
    if (list && !is.null(object$model.index)) {
        return(lapply(object$model.index, function(x) object$coef[x]))
    }
    object$coef
}

##' @export
summary.estimate <- function(object,...) {
    p <- coef(object,messages=0)
    test <- estimate(coef=p,vcov=vcov(object,messages=0),
                     contrast=as.list(seq_along(p)),...)
    object$compare <- test$compare
    object <- object[c("coef","coefmat","vcov","call","ncluster","model.index","compare")]
    class(object) <- "summary.estimate"
    object
}

##' @export
coef.summary.estimate <- function(object,...) {
    object$coefmat
}

##' @export
print.summary.estimate <- function(x,...) {
    print.estimate(x,type=2L,...)
}

##' @export
iid.estimate <- function(x,...) {
    if (is.null(x$iid)) return(NULL)
    dimn <- dimnames(x$iid)
    if (!is.null(dimn)) {
        dimn[[2]] <- names(coef(x))
    } else {
        dimn <- list(NULL,names(coef(x)))
    }
    structure(x$iid,dimnames=dimn)
}

##' @export
model.frame.estimate <- function(formula,...) {
    NULL
}
