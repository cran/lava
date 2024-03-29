##' Add an observed event time outcome to a latent variable model.
##'
##' For example, if the model 'm' includes latent event time variables
##' are called 'T1' and 'T2' and 'C' is the end of follow-up (right censored),
##' then one can specify
##'
##' \code{eventTime(object=m,formula=ObsTime~min(T1=a,T2=b,C=0,"ObsEvent"))}
##'
##' when data are simulated from the model
##' one gets 2 new columns:
##'
##' - "ObsTime": the smallest of T1, T2 and C
##' - "ObsEvent": 'a' if T1 is smallest, 'b' if T2 is smallest and '0' if C is smallest
##'
##' Note that "ObsEvent" and "ObsTime" are names specified by the user.
##'
##' @author Thomas A. Gerds, Klaus K. Holst
##' @keywords survival models regression
##' @examples
##'
##' # Right censored survival data without covariates
##' m0 <- lvm()
##' distribution(m0,"eventtime") <- coxWeibull.lvm(scale=1/100,shape=2)
##' distribution(m0,"censtime") <- coxExponential.lvm(rate=1/10)
##' m0 <- eventTime(m0,time~min(eventtime=1,censtime=0),"status")
##' sim(m0,10)
##'
##' # Alternative specification of the right censored survival outcome
##' ## eventTime(m,"Status") <- ~min(eventtime=1,censtime=0)
##'
##' # Cox regression:
##' # lava implements two different parametrizations of the same
##' # Weibull regression model. The first specifies
##' # the effects of covariates as proportional hazard ratios
##' # and works as follows:
##' m <- lvm()
##' distribution(m,"eventtime") <- coxWeibull.lvm(scale=1/100,shape=2)
##' distribution(m,"censtime") <- coxWeibull.lvm(scale=1/100,shape=2)
##' m <- eventTime(m,time~min(eventtime=1,censtime=0),"status")
##' distribution(m,"sex") <- binomial.lvm(p=0.4)
##' distribution(m,"sbp") <- normal.lvm(mean=120,sd=20)
##' regression(m,from="sex",to="eventtime") <- 0.4
##' regression(m,from="sbp",to="eventtime") <- -0.01
##' sim(m,6)
##' # The parameters can be recovered using a Cox regression
##' # routine or a Weibull regression model. E.g.,
##' \dontrun{
##'     set.seed(18)
##'     d <- sim(m,1000)
##'     library(survival)
##'     coxph(Surv(time,status)~sex+sbp,data=d)
##'
##'     sr <- survreg(Surv(time,status)~sex+sbp,data=d)
##'     library(SurvRegCensCov)
##'     ConvertWeibull(sr)
##'
##' }
##'
##' # The second parametrization is an accelerated failure time
##' # regression model and uses the function weibull.lvm instead
##' # of coxWeibull.lvm to specify the event time distributions.
##' # Here is an example:
##'
##' ma <- lvm()
##' distribution(ma,"eventtime") <- weibull.lvm(scale=3,shape=1/0.7)
##' distribution(ma,"censtime") <- weibull.lvm(scale=2,shape=1/0.7)
##' ma <- eventTime(ma,time~min(eventtime=1,censtime=0),"status")
##' distribution(ma,"sex") <- binomial.lvm(p=0.4)
##' distribution(ma,"sbp") <- normal.lvm(mean=120,sd=20)
##' regression(ma,from="sex",to="eventtime") <- 0.7
##' regression(ma,from="sbp",to="eventtime") <- -0.008
##' set.seed(17)
##' sim(ma,6)
##' # The regression coefficients of the AFT model
##' # can be tranformed into log(hazard ratios):
##' #  coef.coxWeibull = - coef.weibull / shape.weibull
##' \dontrun{
##'     set.seed(17)
##'     da <- sim(ma,1000)
##'     library(survival)
##'     fa <- coxph(Surv(time,status)~sex+sbp,data=da)
##'     coef(fa)
##'     c(0.7,-0.008)/0.7
##' }
##'
##'
##' # The following are equivalent parametrizations
##' # which produce exactly the same random numbers:
##'
##' model.aft <- lvm()
##' distribution(model.aft,"eventtime") <- weibull.lvm(intercept=-log(1/100)/2,sigma=1/2)
##' distribution(model.aft,"censtime") <- weibull.lvm(intercept=-log(1/100)/2,sigma=1/2)
##' sim(model.aft,6,seed=17)
##'
##' model.aft <- lvm()
##' distribution(model.aft,"eventtime") <- weibull.lvm(scale=100^(1/2), shape=2)
##' distribution(model.aft,"censtime") <- weibull.lvm(scale=100^(1/2), shape=2)
##' sim(model.aft,6,seed=17)
##'
##' model.cox <- lvm()
##' distribution(model.cox,"eventtime") <- coxWeibull.lvm(scale=1/100,shape=2)
##' distribution(model.cox,"censtime") <- coxWeibull.lvm(scale=1/100,shape=2)
##' sim(model.cox,6,seed=17)
##'
##' # The minimum of multiple latent times one of them still
##' # being a censoring time, yield
##' # right censored competing risks data
##'
##' mc <- lvm()
##' distribution(mc,~X2) <- binomial.lvm()
##' regression(mc) <- T1~f(X1,-.5)+f(X2,0.3)
##' regression(mc) <- T2~f(X2,0.6)
##' distribution(mc,~T1) <- coxWeibull.lvm(scale=1/100)
##' distribution(mc,~T2) <- coxWeibull.lvm(scale=1/100)
##' distribution(mc,~C) <- coxWeibull.lvm(scale=1/100)
##' mc <- eventTime(mc,time~min(T1=1,T2=2,C=0),"event")
##' sim(mc,6)
##'
##'
##' @export
##' @aliases eventTime<-
##' @param object Model object
##' @param formula Formula (see details)
##' @param eventName Event names
##' @param \dots Additional arguments to lower levels functions
eventTime <- function(object,formula,eventName="status",...) {
    if (missing(formula)) return(object$attributes$eventHistory)
    if (inherits(eventName,"formula")) eventName <- all.vars(eventName)
    ff <- as.character(formula)
    timeName <- all.vars(update.formula(formula,"~1"))
    if (length(timeName)==0){
        timeName <- "observedTime"
        rhs <- ff[[2]]
    }else{
        rhs <- ff[[3]]
    }
    ## rhs <- tolower(rhs)
    latentTimes <- strsplit(rhs,"[(,)]")[[1]]
    if (latentTimes[1]!="min")
        stop(paste("Formula ",formula," does not have the required form, ",
                   "e.g. ~min(T1=1,T2=2,C=0), see (examples in) help(eventTime)."))
    latentTimes <- latentTimes[-1]
    NT <- length(latentTimes)
    events <- vector(NT,mode="character")
    for (lt in seq_len(NT)){
        tmp <- strsplit(latentTimes[lt],"=")[[1]]
        stopifnot(length(tmp) %in% c(1,2))
        if (length(tmp)==1){
            events[lt] <- as.character(lt)
            latentTimes[lt] <- tmp
        }
        else{
            events[lt] <- tmp[2]
            latentTimes[lt] <- tmp[1]
        }
    }
    events <- gsub(" ","",events)
    eventnum <- char2num(events)
    if (all(!is.na(eventnum))) {
        events <- eventnum
    } else {
        events <- gsub("\"","",events)
    }

    addvar(object) <- timeName
    eventTime <- list(names=c(timeName,eventName),
                      latentTimes=gsub(" ","",latentTimes),
                      events=events
                      )

    transform(object,
              y=eventTime$names,
              x=eventTime$latentTimes) <-
        function(z) {
            idx <- apply(z,1,which.min)
            cbind(z[cbind(seq(NROW(z)),idx)],
                  eventTime$events[idx])
        }

    if (is.null(object$attributes$eventHistory)) {
        object$attributes$eventHistory <- list(eventTime)
        names(object$attributes$eventHistory) <- timeName
    } else {
        object$attributes$eventHistory[[timeName]] <- eventTime
    }
    return(object)
}

##' @export
"eventTime<-" <- function(object,...,value) {
    eventTime(object,value,...)
}

## addhook("color_eventHistory","color.hooks")
## color.eventHistory <- function(x,subset=vars(x),...) {
##   return(list(vars=intersect(subset,binary(x)),col="indianred1"))
## }

addhook("plothook_eventHistory","plot.post.hooks")
plothook_eventHistory <- function(x,...) {
    eh <- x$attributes$eventHistory
    for (f in eh) {
        x <- regression(x,to=f$names[1],from=f$latentTimes)
        latent(x) <- f$latentTimes
        kill(x) <- f$names[2]
    }
    timedep <- x$attributes$timedep
    for (i in seq_len(length(timedep))) {
        x <- regression(x,to=names(timedep)[i],from=timedep[[i]])
    }
    return(x)
}

addhook("colorhook_eventHistory","color.hooks")
colorhook_eventHistory <- function(x,subset=vars(x),...) {
    return(list(vars=intersect(subset,unlist(x$attributes$timedep)),col="lightblue4"))
}

addhook("print_eventHistory","print.hooks")
print_eventHistory <- function(x,...) {
    eh <- x$attributes$eventHistory
    timedep <- x$attributes$timedep
    if (is.null(eh) & is.null(timedep)) return(NULL)
    ehnames <- unlist(lapply(eh,function(x) x$names))
    cat("Event History Model\n")
    ff <- formula(x,char=TRUE,all=TRUE)
    R <- c()
    for (f in ff) {
        oneline <- as.character(f);
        y <- gsub(" ","",strsplit(f,"~")[[1]][1])
        if (!(y %in% ehnames)) {
            col1 <- as.character(oneline)
            D <- attributes(distribution(x)[[y]])$family
            col2 <- "Normal"
            if (!is.null(D$family)) col2 <- paste0(D$family)
            if (!is.null(D$link)) col2 <- paste0(col2,"(",D$link,")")
            if (!is.null(D$par)) col2 <- paste0(col2,"(",paste(D$par,collapse=","),")")
            R <- rbind(R,c(col1,"  ",col2))
        }
    }
    for (y in names(eh)) {
        col1 <- paste0(y, " = min(",paste(eh[[y]]$latentTimes,collapse=","),")")
        eh[[y]]$names[2]
        col2 <- paste0(eh[[y]]$names[2], " := {",paste(eh[[y]]$events,collapse=","),"}")
        R <- rbind(R,c(col1,"",col2))
    }
    rownames(R) <- rep("",nrow(R)); colnames(R) <- rep("",ncol(R))
    print(R,quote=FALSE,...)
    cat("\n")
    for (i in seq_len(length(timedep))) {
        cat("Time-dependent covariates:\n\n")
        cat(paste("",names(timedep)[i],"~", paste(timedep[[i]],collapse="+")),"\n")
    }
    TRUE
}

## addhook("simulate_eventHistory","sim.hooks")

## simulate.eventHistory <- function(x,data,...){
##   if (is.null(eventTime(x))) {
##     return(data)
##   }
##   else{
##     for (eh in eventTime(x)) {
##       if (any((found <- match(eh$latentTimes,names(data),nomatch=0))==0)){
##         warning("Cannot find latent time variable: ",
##                 eh$latentTimes[found==0],".")
##       }
##       else{
##         for (v in seq_along(eh$latentTimes)) {
##           if (v==1){ ## initialize with the first latent time and event
##             eh.time <- data[,eh$latentTimes[v]]
##             eh.event <- rep(eh$events[v],NROW(data))
##           } else{ ## now replace if next time is smaller
##             ## in case of tie keep the first event
##             eh.event[data[,eh$latentTimes[v]]<eh.time] <- eh$events[v]
##             eh.time <- pmin(eh.time,data[,eh$latentTimes[v]])
##           }
##         }
##       }
##       data[,eh$names[1]] <- eh.time
##       data[,eh$names[2]] <- eh.event
##     }
##     return(data)
##   }
## }



##' @export
coxWeibull.lvm <- function(scale=1/100,shape=2, param=1) {
    ## proportional hazard (Cox) parametrization.
    ##
    ## Here we parametrize the Weibull distribution
    ## (without covariates) as
    ##
    ## hazard(t) = scale * shape * t^(shape-1)
    ##
    ## The linear predictor (LP) enters like this
    ##
    ## hazard(t) = = scale * exp(LP) * shape * t^(shape-1)
    ##
    ## Thus, we simulate
    ##
    ## T = (scale^{-1} * exp(-LP) * -log(U))^{shape-1})
    ##
    ## The hazard is:
    ## - rising if shape > 1
    ## - declining if shape <1
    ## - constant if shape=1
    ##
    ## scale = exp(b0 + b1*X)
    ##
    ## An alternative parametrization, which can be found in ABGK and
    ## Kalbfleisch&Prentice, can be obtained with \code{param=2}. This leads to
    ## hazard(t) = scale*shape*(scale*t)^shape, and cum.hazard(t) = (scale*t)^shape.
    ## \code{param=3} gives the same parametrization as \code{rweibull}, i.e.,
    ## scale[rweibull] := 1/scale^shape
    if (tolower(param)%in%c("2","kp","abgk"))
       scale <- scale^(shape)
    if (tolower(param)%in%c("3","r","rweibull")) {
       scale <- scale^(-shape)
    }
    f <- function(n,mu=0,Scale=scale,Shape=shape,...) {
        (- log(runif(n)) / (Scale * exp(mu)))^(1/Shape)
    }
    ff <- formals(f)
    expr <- "(- log(runif(n)) / (Scale * exp(mu)))^{1/Shape}"
    if (inherits(scale,"formula")) scale <- all.vars(scale)[1]
    if (is.character(scale)) {
        names(ff)[3] <- scale
        expr <- gsub("Scale",scale,expr)
    }
    if (inherits(shape,"formula")) shape <- all.vars(shape)[1]
    if (is.character(shape)) {
        names(ff)[4] <- shape
        expr <- gsub("Shape",shape,expr)
    }
    formals(f) <- ff
    e <- parse(text=expr)
    body(f) <- as.call(c(as.name("{"), e))
    attr(f,"family") <- list(family="weibull",
                             regression="PH",
                             par=c(shape=shape,scale=scale))
    return(f)
}


##' @export
coxExponential.lvm <- function(scale=1,rate=1,timecut){
    if (missing(timecut)) {
        return(coxWeibull.lvm(shape=1, rate))
    }
    if (NROW(rate)>length(timecut))
        stop("Number of time-intervals (cuts) does not agree with number of rate parameters (beta0)")
    par <- paste(timecut,rate,sep=":")
    if (is.matrix(rate)) par <- "..."
    timecut <- c(timecut,Inf)
    f <- function(n,mu=0,...) {
        Ai <- function() {
            vals <- matrix(0,ncol=length(timecut)-1,nrow=n)
            ival <- numeric(n)
            if (is.matrix(rate)) {
                mu <- cbind(mu[,1],cbind(1,as.matrix(mu[,-1]))%*%t(rate))
                rate <- rep(1,length(timecut)-1)
            }
            for (i in seq(length(timecut)-1)) {
                u <- -log(runif(n)) ##rexp(n,1)
                if (NCOL(mu)>1) {
                    vals[,i] <-  timecut[i] + u*exp(-mu[,1]-mu[,i+1])/(rate[i])
                } else {
                    vals[,i] <-  timecut[i] + u*exp(-mu)/(rate[i])
                }
                idx <- which(vals[,i]<=timecut[i+1] & ival==0)
                ival[idx] <- vals[idx,i]
            }
            ival
        }
        Ai()
    }
    attributes(f)$family <- list(family="CoxExponential",par=par)
    return(f)
}

##' @export
aalenExponential.lvm <- function(rate=1,timecut=0){
    if (missing(timecut)==1) {
        return(coxWeibull.lvm(shape=1,rate))
    }

    if (length(rate)>length(timecut))
        stop("Number of time-intervals (cuts) does not agree with number of rate parameters (beta0)")
    par <- paste(timecut,rate,sep=":")
    if (is.matrix(rate)) par <- "..."
    timecut <- c(timecut,Inf)
    f <- function(n,mu,...) {
        Ai <- function() {
            vals <- matrix(0,ncol=length(timecut)-1,nrow=n)
            ival <- numeric(n)
            if (is.matrix(rate)) {
                mu <- cbind(mu[,1],cbind(1,as.matrix(mu[,-1]))%*%t(rate))
                rate <- rep(1,length(timecut)-1)
            }
            for (i in seq(length(timecut)-1)) {
                u <- -log(runif(n)) ##rexp(n,1)
                if (NCOL(mu)>1) {
                    vals[,i] <-  timecut[i] + u/(rate[i]+mu[,1]+mu[,i+1])
                } else {
                    vals[,i] <-  timecut[i] + u/(rate[i]+mu)
                }
                idx <- which(vals[,i]<=timecut[i+1] & ival==0)
                ival[idx] <- vals[idx,i]
            }
            ival
        }
        Ai()
    }
    attributes(f)$family <- list(family="aalenExponential",par=par)
    return(f)
}


##' @export
coxGompertz.lvm <- function(shape=1,scale) {
    f <- function(n,mu,var,...) {
        (1/shape) * log(1 - (shape/scale) * (log(runif(n)) * exp(-mu)))
    }
    attr(f,"family") <- list(family="gompertz",par=c(shape,scale))
    return(f)
}
