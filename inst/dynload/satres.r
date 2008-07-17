library(deSolve)
dyn.load("satres.dll")

## Dose is the Dose in mg/kg
## Doseint is NA for single dose, interval
## between doses in hours for repeated dosing.
## Other parms as in satres.c
defParms <- 
  c(Vc=0.0027, Vt=0.0545, kd=0.00059/0.0027, ka=0.537,
    Tm=860.9, KT=0.0015, kfil=0.6830/0.0027,
    Vfil=0.01, free=0.02, BW=.025, Dose=NA, Doseint=NA,
    Qd=NA, TDose=NA)

## initparms is called as, for example
## P <- initparms(list(Dose=60, Doseint=24, Vc=0.0030))
## Gives a parameter list that the model can use, for 60 mg/kg
## every 24 hours dosing, setting Vc to 0.003 L

## newParms is a list with parameter names
initparms <- function(newParms=NULL){
  if (!is.null(newParms)) {
    ldots <- as.list(newParms)
    if (!all(names(ldots) %in% names(defParms)[1:12])) {
      stop("illegal parameter name")
    }
  }
  Parms <- defParms
  if (!is.null(newParms)) Parms[names(ldots)] <- unlist(ldots)
  lParms <- as.list(Parms)
  Parms["Qd"] <- with(lParms, kd * Vc)
  Parms["TDose"] <- with(lParms, Dose * BW)
  Parms
}

## initState returns the initialized state vector.  If new is TRUE,
## the vector state is ignored, and the return value is the initial
## value of state.  If new is FALSE, then on input, state is the
## current state of the system, and initState increments the
## fourth compartment (the gut) by the dose value
initState <- function(Parms,state,new=TRUE){
  if (new){
    state <- c(rep(0,3),Parms["TDose"], 0)
  } else {
    state[4] <- state[4] + Parms["TDose"]
  }
  state
}

## pfoasat runs the model.  On input,
##   - Times is a vector of time values
##     at which model results are desired.
##   - newParms is a list like the input
##     to initparms, above.
##   - method is a string giving the solution method to use
##     see the documentation for deSolve::ode for details
##     there.  the elipsis (...) is for additional arguments
##     to the odesolver (see ode and the individual methods
##     for details).
## The return value is a matrix of values.  Column 1 is the
## time vector, Columns 2 - 5 are the concentrations in
## compartments 1 - 4 (just before dosing, in the case of repeated
## dosing).
##
## Example: to match the 7 and 17 day 20 mg/kg repeated dosing
## using lsode:
## out <- pfoasat(24 * c(0,7,17), newParms=list(Dose=20, Doseint=24))
## when finished, you can unload the dll with
## dyn.unload("satres")

pfoasat <- function(Times, newParms, method="lsode", ...){
  ### Initialize parameters
  Parms <- initparms(newParms)
  ### Single or repeate dosing?
  if (is.na(Parms["Doseint"])){ ## single dose
    y <- initState(Parms,new=TRUE)
    ode(y, Times, "derivs", Parms, method=method,
        dllname="satres", initfunc="initmod",
        nout=1, outnames="Total", ...)
  } else { ## repeated dose
    ## restart the solver at the elements of breaks
    breaks <- seq(min(Times),max(Times),by=Parms["Doseint"])
    segments <- split(Times, cut(Times, breaks))
    ## put the start and stop time for each segment in the segment,
    ## then reduce redundant values
    tmp <- matrix(nrow=0, ncol=7)
    for (i in seq(along=segments)) {
      tx <- unique(c(breaks[i],segments[[i]],breaks[i+1]))
      y <- initState(Parms,tmp[nrow(tmp),2:6], new=(i==1))
      tmp2 <- ode(y, tx, "derivs", Parms, method=method,
                  dllname="satres", initfunc="initmod",
                  nout=1, outnames="Total", ...)
      tmp <- rbind(tmp,tmp2)
    }
    tmp <- tmp[!duplicated(tmp[,1]),]
    tmp[tmp[,1] %in% Times,]
  }
}
#out <- pfoasat(24*c(0,7,17), newParms=list(Dose=20, Doseint=24))
out <- pfoasat(24*seq(0,17,0.1), newParms=list(Dose=20, Doseint=24))

plot(out,type="l")
