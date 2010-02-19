### ============================================================================
###
### timelags and delay differential equations
###
### ============================================================================

## =============================================================================
## lagged values and derivates are obtained in the R-code via functions 
## lagvalue and lagderiv 
## =============================================================================

lagvalue <- function (t, nr=NULL) {
  if (is.null(nr)) nr <- 0
  out <- .Call("getLagValue", t = t, PACKAGE = "deSolve", as.integer(nr))
  return(out)
}

lagderiv <- function (t, nr=NULL) {
  if (is.null(nr)) nr <- 0
  out <- .Call("getLagDeriv", t = t, PACKAGE = "deSolve", as.integer(nr))
  return(out)
}

### ============================================================================
### solving Delay Differential Equations
### ============================================================================

dede <- function(y, times, func=NULL, parms, method = c( "lsoda", "lsode", 
    "lsodes", "lsodar", "vode", "daspk"), control=NULL,  ...) {
    if (is.null(control)) control <- list(mxhist = 1e4)

    if (is.null(method)) 
        method <- "lsoda"
    else if (is.function(method)) 
        res <- method(y, times, func, parms, lags = control, ...)
    else res <- switch(match.arg(method), 
       lsoda = lsoda(y, times, func, parms, lags = control, ...), 
       vode = vode(y, times, func, parms, lags = control, ...), 
       lsode = lsode(y, times, func, parms, lags = control, ...), 
       lsodes = lsodes(y, times, func, parms, lags = control, ...), 
       lsodar = lsodar(y, times, func, parms, lags = control, ...), 
       daspk = daspk(y, times, func, parms, lags = control, ...))
    return(res)
}

