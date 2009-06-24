aquaphy <- function(times, y, parms, ...) {

  if (length(y) != 4)
    stop ("length of state variable vector should be 4")
  if (length(parms) != 19)
    stop ("length of parameter vector should be 19")

  names(y) <- c("DIN","PROTEIN","RESERVE","LMW")
  outnames <- c("PAR","TotalN","PhotoSynthesis",
                "NCratio","ChlCratio","Chlorophyll")
  ode(y,times,dllname="deSolve",
      func="aquaphy",initfunc="iniaqua",
      parms=parms,nout=6,outnames=outnames,...)
}

