diagnostics <- function(obj) {
  Attr <- attributes(obj)
  type <- Attr$type
  if (is.null(type))
    stop("cannot print ODE characteristics; output not of correct type")
  istate <- Attr$istate
  rstate <- Attr$rstate
    cat("\n--------------------\n")
  if (type == "vode")
    cat("vode return code")
  else if (type == "lsoda")
    cat("lsoda return code")
  else   if (type == "lsodes")
    cat("lsodes return code")
  else if (type == "lsode")
    cat("lsode return code")
  else if (type=="lsodar")
    cat("lsodar return code")
  else if (type == "daspk")
    cat("daspk return code")
  else if (type == "rk")
    cat("rk return code")

    cat("\n--------------------\n")

  idid <- istate[1]
  cat(paste("\n  idid = ", idid), "\n")

  if (type %in% c("vode", "lsoda", "lsodes", "lsode", "lsodar")) {
    if (type == "lsodar" && idid ==2 )
      cat("  Integration was successful but no root was found.\n") else
    if (idid == 2)  cat("  Integration was successful.\n") else
    if (idid == 3)  cat("  Integration was successful and a root was found before reaching the end.\n") else
    if (idid == -1) cat("  Excess work done on this call. (Perhaps wrong Jacobian type MF.)\n") else
    if (idid == -2) cat("  Excess accuracy requested. (Tolerances too small.)\n") else
    if (idid == -3) cat("  Illegal input detected. (See printed message.)\n") else
    if (idid == -4) cat("  Repeated error test failures. (Check all input.)\n") else
    if (idid == -5) cat("  Repeated convergence failures. (Perhaps bad Jacobian supplied or wrong choice of MF or tolerances.)\n") else
    if (idid == -6) cat("  Error weight became zero during problem. (Solution component i vanished, and ATOL or ATOL(i) = 0.)\n") else
    if (type == lsodes && idid == -7)
       cat("  A fatal error came from sparse solver CDRV by way of DPRJS or DSOLSS.\n") else
    if (idid == -7) cat("  Work space insufficient to finish (see messages).\n")
  } else if (type == "daspk") {
    if (idid > 0)   {
      cat ("  *** TASK COMPLETED ***\n")
      if (idid == 1) cat("  A step was successfully taken in the intermediate-output mode.  The code has not yet reached TOUT.\n")
      if (idid == 2) cat("  The integration to TSTOP was successfully completed (T = TSTOP) by stepping exactly to TSTOP.\n")
      if (idid == 3) cat("  The integration to TOUT was successfully completed (T = TOUT) by stepping past TOUT. Y(*) and YPRIME(*) are obtained by interpolation.\n")
      if (idid == 4) cat("  The initial condition calculation, with INFO(11) > 0, was successful, and INFO(14) = 1. No integration steps were taken, and the solution is not considered to have been started.\n")
    } else if (idid < 0 & idid > -33)  {
      cat ("  *** TASK INTERRUPTED ***\n")
      if (idid == -1) cat("  A large amount of work has been expended (about 500 steps).\n") else
      if (idid == -2) cat("  The error tolerances are too stringent.\n") else
      if (idid == -3) cat("  The local error test cannot be satisfied because a zero component in ATOL was specified and the corresponding computed solution component is zero.  Thus, a pure relative error test is impossible for this component.\n") else
      if (idid == -5) cat("  There were repeated failures in the evaluation or processing of the preconditioner (in jacfunc).\n") else
      if (idid == -6) cat("  DDASPK had repeated error test failures on the last attempted step.\n") else
      if (idid == -7) cat("  The nonlinear system solver in the time integration could not converge.\n") else
      if (idid == -8) cat("  The matrix of partial derivatives appears to be singular (direct method).\n") else
      if (idid == -9) cat("  The nonlinear system solver in the time integration failed to achieve convergence, and there were repeated error test failures in this step.\n") else
      if (idid == -10) cat("  The nonlinear system solver in the time integration failed to achieve convergence because IRES was equal to -1.\n") else
      if (idid == -11) cat("  IRES = -2 was encountered and control is being returned to the calling program.\n") else
      if (idid == -12) cat("  DDASPK failed to compute the initial Y, YPRIME.\n") else
      if (idid == -13) cat("  Unrecoverable error encountered inside user's PSOL routine, and control is being returned to the calling program.\n") else
      if (idid == -14) cat("  The Krylov linear system solver could not achieve convergence.\n")
    } else if (idid ==-33)  {
      cat ("  *** TASK TERMINATED ***\n")
      cat("  The code has encountered trouble from which it cannot recover.  A message is printed explaining the trouble and control is returned to the calling program.\n")
    }
  } else if (type == "rk") {
      if (idid == 0)  cat("  Integration was successful.\n") else
      if (idid == -1) cat("  A large amount of work has been expended. Increase maxsteps.\n") else
      if (idid == -2) cat("  Excess accuracy requested. Tolerances too small.\n") else
        cat("  rk returned with undefined return code.\n")
  } else {
    warning("Unknown return type.")
  }

#### istate

  cat("\n--------------------\n")
  cat("ISTATE values\n")
  cat("--------------------\n")
  df <- c( "The return code (istate):",
           "The number of steps taken for the problem so far:",
           "The number of function evaluations for the problem so far:",
           "The number of Jacobian evaluations so far:",
           "The method order last used (successfully):",
           "The order to be attempted on the next step:",
           "If istate=-4,-5: the largest comp in error vector",
           "The length of rwork actually required:",
           "The length of iwork actually required:",
           "The number of matrix LU decompositions so far:",
           "The number of nonlinear (Newton) iterations so far:",
           "The number of convergence failures of the solver so far ",
           "The number of error test failures of the integrator so far:")
  if (type == "vode")  {
    ii <- c(1,12:23)
  } else if (type %in% c("lsoda", "lsodar")) {
    df[4] <- "The number of Jacobian evaluations and LU decompositions so far:"
    df[10]<- "The method indicator for the last succesful step, 1=adams (nonstiff), 2= bdf (stiff):"
    df[11]<- "The current method indicator to be attempted on th next step, 1=adams (nonstiff), 2= bdf (stiff):"
    ii <- c(1,12:21)

  } else if (type == "lsodes") {
    df[4] <- "The number of Jacobian evaluations and LU decompositions so far:"
    df[10]<- "The number of nonzero elements in the sparse Jacobian:"
    ii <- c(1,12:20)

  } else if (type == "lsode") {
    df[4] <- "The number of Jacobian evaluations and LU decompositions so far:"
    ii <- c(1,12:19)

  } else if (type == "daspk") {
    df <- c( "The return code (idid):",
             "The order of the method to be attempted on th next step:",
             "The order of the method used on the last step:",
             "The number of steps taken for the problem so far:",
             "The number of res evaluations for the problem so far:",
             "The number of Jacobian evaluations for the problem so far:",
             "The total number of error test failures so far:",
             "The total number of nonlinear convergence failures so far:",
             "The number of convergence failures of the linear iteration so far:",
             "The length of iwork actually required:",
             "The length of rwork actually required:",
             "The number of nonlinear iterations so far:",
             "The number of linear (Krylov) iterations so far ",
             "The number of psol calls so far:")
    ii <- c(1,8:9,12:22)
  } else if (type == "rk") {
    df <- c("The return code:",
            "The number of steps taken for the problem so far:",
            "The number of function evaluations for the problem so far:",
            "The order of the method:")
    ii <- c(1, 12, 13, 15)
  }

### rstate
  printmessage(df[1:length(ii)], istate[ii])

  if (type != "rk") {
    cat("--------------------\n")
    cat("RSTATE values\n")
    cat("--------------------\n")
    ii <- 1:4
    df <- c( "The step size in t last used (successfully):",
      "The step size to be attempted on the next step:",
      "The current value of the independent variable which the solver has actually reached:",
      "Tolerance scale factor > 1.0 computed when a request for too much accuracy was detected:")

    if (type %in% c("lsoda", "lsodar")) {
      df <- c(df,"The value of t at the time of the last method switch, if any:")
      ii <- 1:5
    }
    printmessage(df[ii], rstate[ii])
  }
  invisible(list(istate=istate, rstate=rstate)) # return something more useful; may change in the future
}
