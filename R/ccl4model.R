ccl4model <- function(times, # time sequence for which the model has to be integrated
                      y,     # intial values of the state variables
                      parms, # parameter values
                    ...)
{ 
  if (length(y) != 7) stop ("length of state variable vector should be 7")
  if (length(parms) != 21) stop ("length of parameter vector should be 21")

  names(y) <- c("AI","AAM","AT","AF","AL","CLT","AM")

  ode(y=y,dllname="deSolve",func="derivsccl4",
        initfunc = "initccl4",parms=parms,
        times=times,nout=3,outnames=c("DOSE","MASS","CP"),...)
}

