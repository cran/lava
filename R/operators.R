##' For matrices a block-diagonal matrix is created. For all other
##' data types he operator is a wrapper of \code{paste}.
##'
##' Concatenation operator
##' @aliases %++%
##' @rdname op_concat
##' @usage x \%++\% y
##' @title Concatenation operator
##' @param x First object
##' @param y Second object of same class
##' @author Klaus K. Holst
##' @keywords utilities misc
##' @seealso \code{blockdiag}, \code{\link{paste}}, \code{\link{cat}},
##' @examples
##' ## Block diagonal
##' matrix(rnorm(25),5)%++%matrix(rnorm(25),5)
##' ## String concatenation
##' "Hello "%++%" World"
##' ## Function composition
##' f <- log %++% exp
##' f(2)
##' @export
`%++%` <- function(x,y) UseMethod("%++%",y)

## ##' @export
## `%+%` <- function(x,y) UseMethod("%+%",y)

##' @export
`%++%.default` <- function(x,y) paste0(x,y)

##' @export
`%++%.character` <- function(x,y) paste0(x,y)

##' @export
`%++%.matrix` <- function(x,y) blockdiag(x,y)

##' @export
`%++%.function` <- function(x,y) function(...) x(y(...))


notin <- Negate(get("%in%"))
##' Matching operator (x not in y) oposed to the \code{\%in\%}-operator (x in y)
##'
##' Matching operator
##' @rdname op_match
##' @aliases %ni% %in.open% %in.closed%
##' @usage x \%ni\% y
##' @param x vector
##' @param y vector of same type as \code{x}
##' @return A logical vector.
##' @author Klaus K. Holst
##' @seealso \code{\link{match}}
##' @keywords utilities misc
##' @examples
##'
##' 1:10 %ni% c(1,5,10)
##'
##' @export
"%ni%" <- function(x,y) notin(x,y)

## function(x,y) {
##   is.na(match(x,y))
## }

##' @export
"%in.open%" <- function(x, y) {
  if (length(y) == 1) y <- c(y, y)
  if (length(y) != 2 || !is.numeric(y)) stop("rhs should be a range (numeric vector of length 2)")
  x > y[1] & x < y[2]
}

##' @export
"%in.closed%" <- function(x,y) {
  if (length(y) == 1) y <- c(y, y)
  if (length(y) != 2 || !is.numeric(y)) stop("rhs should be a range (numeric vector of length 2)")
  x >= y[1] & x <= y[2]
}
