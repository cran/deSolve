### ============================================================================
### Interface to C code for Euler's ODE solver
### with fixed step size and without interpolation, see helpfile for details.
### ============================================================================

euler <- function(y, times, func, parms, verbose = FALSE, ynames = TRUE,
  dllname = NULL, initfunc = dllname, initpar = parms,
  rpar = NULL,  ipar = NULL, nout = 0, outnames = NULL, forcings = NULL,
  initforc = NULL, fcontrol = NULL, ...) {

    ## check input
    checkInputEuler(y, times, func, dllname)
    n <- length(y)

    ## Model as shared object (DLL)?
    Ynames   <- attr(y, "names")
    Initfunc <- NULL
    flist    <-list(fmat = 0, tmat = 0, imat = 0, ModelForc = NULL)
    Nstates <- length(y) # assume length of states is correct

    if (is.character(func)) {
      DLL <- checkDLL(func, NULL, dllname,
                    initfunc, verbose, nout, outnames)

      Initfunc <- DLL$ModelInit
      Func     <- DLL$Func
      Nglobal  <- DLL$Nglobal
      Nmtot    <- DLL$Nmtot

      if (! is.null(forcings))
        flist <- checkforcings(forcings, times, dllname, initforc, verbose, fcontrol)

      rho <- NULL
      if (is.null(ipar)) ipar <- 0
      if (is.null(rpar)) rpar <- 0

    } else {
      initpar <- NULL # parameter initialisation not needed if function is not a DLL
      rho <- environment(func)
      ## func and jac are overruled, either including ynames, or not
      ## This allows to pass the "..." arguments and the parameters
      if(ynames) {
        Func   <- function(time, state, parms) {
          attr(state, "names") <- Ynames
          func   (time, state, parms, ...)
        }
      } else { # no ynames ...
        Func   <- function(time, state, parms)
          func   (time, state, parms, ...)
      }

      ## Call func once to figure out whether and how many "global"
      ## results it wants to return and some other safety checks
      FF <- checkFuncEuler(Func, times, y, parms, rho, Nstates)
      Nglobal <- FF$Nglobal
      Nmtot   <- FF$Nmtot

    }

    ## the CALL to the integrator
    out <- .Call("call_euler", as.double(y), as.double(times),
                 Func, Initfunc, parms, as.integer(Nglobal), rho, as.integer(verbose),
                 as.double(rpar), as.integer(ipar), flist, PACKAGE = "deSolve")

    ## saving results
    out <- saveOutrk(out, y, n, Nglobal, Nmtot,
                     iin = c(1, 12, 13, 15), iout = c(1:3, 18))
    ## testing code: 'call_euler_t' is a version with transposed data structure in memory
    #    out <- .Call("call_euler_t", as.double(y), as.double(times),
    #                 Func, Initfunc, parms, as.integer(Nglobal), rho, as.integer(verbose),
    #                 as.double(rpar), as.integer(ipar), flist, PACKAGE = "deSolve")
    #    out <- saveOutrk(out, y, n, Nglobal, Nmtot,
    #                 iin = c(1, 12, 13, 15), iout = c(1:3, 18), transpose = TRUE)
    if (verbose) diagnostics(out)
    attr(out, "type") <- "rk"
    out
}
