## ========================================================================
## General functions of deSolve
## ========================================================================

timestep <- function (prev = TRUE) {
  out <- .Call("getTimestep", PACKAGE = "deSolve")
  if (prev)
    return(out[1])
  else
    return(out[2])
}


## ========================================================================
## Check solver input - livermore solvers and rk
## ========================================================================

checkInput <- function(y, times, func, rtol, atol,
  jacfunc, tcrit, hmin, hmax, hini, dllname, jacname = "jacfunc")
{
  if (!is.numeric(y))     stop("`y' must be numeric")
  n <- length(y)
  if (! is.null(times) && !is.numeric(times))
    stop("`times' must be NULL or numeric")
  if (!is.function(func) && !is.character(func))
    stop("`func' must be a function or character vector")
  if (is.character(func) && (is.null(dllname) || !is.character(dllname)))
    stop("specify the name of the dll or shared library where func can be found (without extension)")
  if (!is.numeric(rtol))  stop("`rtol' must be numeric")
  if (!is.numeric(atol))  stop("`atol' must be numeric")
  if (!is.null(tcrit) & !is.numeric(tcrit)) stop("`tcrit' must be numeric")
  if (!is.null(jacfunc) && !(is.function(jacfunc) || is.character(jacfunc)))
    stop(paste(jacname," must be a function or character vector"))
  if (length(atol) > 1 && length(atol) != n)
    stop("`atol' must either be a scalar, or as long as `y'")
  if (length(rtol) > 1 && length(rtol) != n)
    stop("`rtol' must either be a scalar, or as long as `y'")
  if (!is.numeric(hmin))   stop("`hmin' must be numeric")
  if (hmin < 0)            stop("`hmin' must be a non-negative value")
  if (is.null(hmax))
    hmax <- if (is.null(times)) 0 else max(abs(diff(times)))
  if (!is.numeric(hmax))   stop("`hmax' must be numeric")
  if (hmax < 0)            stop("`hmax' must be a non-negative value")
  if (hmax == Inf)  hmax <- 0
  if (!is.null(hini))
   if(hini < 0)            stop("`hini' must be a non-negative value")
  return(hmax)
}

## ========================================================================
## Check solver input - euler and rk4
## ========================================================================

checkInputEuler <- function (y, times, func, dllname) {
    if (!is.numeric(y))  stop("`y' must be numeric")
    n <- length(y)
    if (! is.null(times) && !is.numeric(times))
        stop("`times' must be NULL or numeric")
    if (!is.function(func) && !is.character(func))
      stop("`func' must be a function or character vector")
    if (is.character(func) && (is.null(dllname) || !is.character(dllname)))
      stop("You need to specify the name of the dll or shared library where func can be found (without extension)")
}

## ========================================================================
## Check ode function call - livermore solvers
## ========================================================================

checkFunc<- function (Func2, times, y, rho) {
    ## Call func once to figure out whether and how many "global"
    ## results it wants to return and some other safety checks
    tmp <- eval(Func2(times[1], y), rho)
    if (!is.list(tmp))
      stop("Model function must return a list\n")
    if (length(tmp[[1]]) != length(y))
      stop(paste("The number of derivatives returned by func() (",
                 length(tmp[[1]]),
                 ") must equal the length of the initial conditions vector (",
                 length(y), ")", sep = ""))

    ## use "unlist" here because some output variables are vectors/arrays
    Nglobal <- if (length(tmp) > 1)
      length(unlist(tmp[-1]))  else 0
    ## Karline: changed this:
    ##  Nmtot is now a list with names, dimensions,... for 1-D, 2-D vars
    Nmtot <- list()
    Nmtot$colnames <- attr(unlist(tmp[-1]), "names")

    Nmtot$lengthvar <- unlist(lapply(tmp, length))
    if (length(Nmtot$lengthvar) < Nglobal+1){
      Nmtot$dimvar <- lapply(tmp[-1], dim)
    }
   return(list(Nglobal = Nglobal, Nmtot = Nmtot))
}

## ========================================================================
## Check event function calls
## ========================================================================

checkEventFunc<- function (Func, times, y, rho) {
    ## Call func once
    tmp <- eval(Func(times[1], y), rho)
    if (length(tmp) != length(y))
      stop(paste("The number of values returned by events$func() (",
                 length(tmp),
                 ") must equal the length of the initial conditions vector (",
                 length(y), ")", sep = ""))
    if (!is.vector(tmp))
      stop("The event function 'events$func' must return a vector\n")
}

## ========================================================================
## Check ode function call - euler and rk solvers
## ========================================================================

checkFuncEuler<- function (Func, times, y, parms, rho, Nstates) {
      ## Call func once to figure out whether and how many "global"
      ## results it wants to return and some other safety checks
      tmp <- eval(Func(times[1], y, parms), rho)
      if (!is.list(tmp)) stop("Model function must return a list\n")
      if (length(tmp[[1]]) != Nstates)
        stop(paste("The number of derivatives returned by func() (",
                   length(tmp[[1]]),
                   "must equal the length of the initial conditions vector (",
                   Nstates, ")", sep=""))

      ## use "unlist" because output variables can be vectors/arrays
      Nglobal <- if (length(tmp) > 1)
        length(unlist(tmp[-1])) else 0
      Nmtot <- list()
      Nmtot$colnames <- attr(unlist(tmp[-1]), "names")
      Nmtot$lengthvar <- unlist(lapply(tmp, length))
      if (length(Nmtot$lengthvar) < Nglobal+1){
        Nmtot$dimvar <- lapply(tmp[-1], dim)
      }
   return(list(Nglobal = Nglobal, Nmtot = Nmtot))

}

## ========================================================================
## check ode DLL input
## ========================================================================

checkDLL <- function (func, jacfunc, dllname,
                      initfunc, verbose, nout, outnames, JT = 1) {

    if (sum(duplicated (c(func, initfunc, jacfunc))) > 0)
      stop("func, initfunc, or jacfunc cannot be the same")
    ModelInit <- NA
    if (! is.null(initfunc))  # to allow absence of initfunc
      if (inherits (initfunc, "CFunc"))
        ModelInit <- body(initfunc)[[2]]
      else if (is.loaded(initfunc, PACKAGE = dllname, type = "") ||
        is.loaded(initfunc, PACKAGE = dllname, type = "Fortran"))  {
        ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
      } else if (initfunc != dllname && ! is.null(initfunc))
        stop(paste("'initfunc' not loaded ", initfunc))

    ## Easier to deal with NA in C-code
    if (is.null(initfunc)) ModelInit <- NA

    ## copy value of func to funcname
    ## check to make sure it describes a function in a loaded dll

    funcname <- func
    ## get the pointer and put it in func

    if (inherits (func, "CFunc"))
        Func <- body(func)[[2]]
    else if(is.loaded(funcname, PACKAGE = dllname)) {
      Func <- getNativeSymbolInfo(funcname, PACKAGE = dllname)$address
    } else stop(paste("dyn function 'func' not loaded", funcname))

    ## Finally, is there a Jacobian?
    if (!is.null(jacfunc)) {
      if (!is.character(jacfunc))
        switch (JT,
          stop("If 'func' is dynloaded, so must 'jacfunc' be"),
          stop("If 'func' is dynloaded, so must 'jacvec' be")
        )
      jacfuncname <- jacfunc
      if (inherits (jacfunc, "CFunc"))
        JacFunc <- body(jacfunc)[[2]]

      else if(is.loaded(jacfuncname, PACKAGE = dllname))  {
        JacFunc <- getNativeSymbolInfo(jacfuncname, PACKAGE = dllname)$address
      } else stop(paste("cannot integrate: jac function not loaded ", jacfunc))
    } else JacFunc <- NULL
    Nglobal <- nout
    Nmtot <- list()
    if (is.null(outnames))
      { Nmtot$colnames   <- NULL} else
    if (length(outnames) == nout)
      { Nmtot$colnames   <- outnames} else
    if (length(outnames) > nout)
      Nmtot$colnames <- outnames[1:nout] else
    Nmtot$colnames <- c(outnames,(length(outnames)+1):nout)

    cnames <- outnames
    unames <- unique(outnames)
    if (length(cnames) > length(unames))
      Nmtot$lengthvar <- c(NA,
        sapply (unames, FUN = function(x) length(which(cnames == x))))

    ## thpe: set extra class to be checked by solvers
    ret <- list(ModelInit = ModelInit, Func = Func, JacFunc = JacFunc,
                Nglobal = Nglobal, Nmtot = Nmtot)
    class(ret) <- c("deSolve.symbols", "list")
   return(ret)
}

## =============================================================================
## print integration task
## =============================================================================
printtask <- function(itask, func, jacfunc) {
    printM("\n--------------------")
    printM("Time settings")
    printM("--------------------\n")
    if (itask==1) printM("  Normal computation of output values of y(t) at t = TOUT") else
    if (itask==2) printM("  Take one step only and return.")                          else
    if (itask==3) printM("  istop at the first internal mesh point at or beyond t = TOUT and return. ")  else
    if (itask==4) printM("  Normal computation of output values of y(t) at t = TOUT but without overshooting t = TCRIT.") else
    if (itask==5) printM("  Take one step, without passing TCRIT, and return.")
    printM("\n--------------------")
    printM("Integration settings")
    printM("--------------------\n")
    if (is.character(func)) printM(paste("  Model function a DLL: ", func)) else
    printM("  Model function an R-function: ")
    if (is.character(jacfunc)) printM(paste ("  Jacobian specified as a DLL: ", jacfunc)) else
    if (!is.null(jacfunc))     printM("  Jacobian specified as an R-function: ") else
    printM("  Jacobian not specified")
    cat("\n")
}

## =============================================================================
## Make Istate vector similar for all solvers.
## =============================================================================

setIstate <- function(istate, iin, iout)
{
  IstateOut <- rep(NA, 21)
  IstateOut[iout] <- istate[iin]
  IstateOut
}


## =============================================================================
## Output cleanup  - for the Livermore solvers
## =============================================================================

saveOut <- function (out, y, n, Nglobal, Nmtot, func, Func2,
  iin, iout, nr = 4) {
  troot  <- attr(out, "troot")
  istate <- attr(out, "istate")
  istate <- setIstate(istate,iin,iout)

  valroot  <- attr(out, "valroot")
  indroot <- attr(out, "indroot")

  Rstate <- attr(out, "rstate")
  rstate <- rep(NA,5)
  rstate[1:nr] <- Rstate[1:nr]

  nm <- c("time",
          if (!is.null(attr(y, "names"))) names(y) else as.character(1:n))
  if (Nglobal > 0) {
    nm <- c(nm,
            if (!is.null(Nmtot$colnames))
              Nmtot$colnames else as.character((n+1) : (n + Nglobal)))
  }
  attr(out,"istate") <- istate
  attr(out, "rstate") <- rstate
  if (! is.null(Nmtot$lengthvar))
    if (is.na(Nmtot$lengthvar[1]))Nmtot$lengthvar[1] <- length(y)
  attr(out, "lengthvar") <- Nmtot$lengthvar
  if (! is.null(troot)) attr(out, "troot") <-  troot
  if (! is.null(valroot)) attr(out, "valroot") <- matrix(nrow = n, valroot)
  if (! is.null(indroot)) attr(out, "indroot") <- indroot

  ii <- if (is.null(Nmtot$dimvar))
    NULL else !(unlist(lapply(Nmtot$dimvar, is.null))) # variables with dimension
  if (sum(ii) >0)
    attr(out, "dimvar") <- Nmtot$dimvar[ii]     # dimensions that are not null
  class(out) <- c("deSolve", "matrix")          # a differential equation
  dimnames(out) <- list(nm, NULL)
  return (t(out))
}

## =============================================================================
## Output cleanup  - for the Runge-Kutta solvers
## =============================================================================

saveOutrk <- function(out, y, n, Nglobal, Nmtot, iin, iout, transpose = FALSE)  {
  ## Names for the outputs
  nm <- c("time",
    if (!is.null(attr(y, "names"))) names(y) else as.character(1:n)
   )
  ## Global outputs
  if (Nglobal > 0) {
    nm  <- c(nm,
      if (!is.null(Nmtot$colnames))
        Nmtot$colnames else as.character((n + 1) : (n + Nglobal))
    )
  }
  ## Column names and state information
  dimnames(out) <- list(NULL, nm)
  istate <- attr(out, "istate")
  istate <- setIstate(istate, iin, iout)
  attr(out,"istate") <- istate
  if (! is.null(Nmtot$lengthvar))
    if (is.na(Nmtot$lengthvar[1])) Nmtot$lengthvar[1] <- length(y)
  attr(out, "lengthvar") <- Nmtot$lengthvar

  ii <- if (is.null(Nmtot$dimvar))
    NULL else !(unlist(lapply(Nmtot$dimvar, is.null))) # variables with dimension
  if (sum(ii) >0)
    attr(out, "dimvar") <- Nmtot$dimvar[ii] # only those which are not null
  class(out) <- c("deSolve", "matrix")      # output of a differential equation
  if (transpose)
    return(t(out))
  else
    return(out)
}

