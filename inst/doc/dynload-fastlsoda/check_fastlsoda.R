## Check of a feture, introduced in deSolve 1.35, where
## shared object symbol identification is separated from integration
##
## The function \code{checkDLL} is normally called internally by the solver
## functions. It can be used to avoid overhead, when a small compiled
## model with a low number of integration steps is repeatedly called.
## The feature is currently only available for the \code{lsoda} solver and
## not yet implemented for symbols of roots and forcings.

library("deSolve")

parameters <- c(a = -8/3, b = -10, c =  28)
state <- c(X = 1, Y = 1, Z = 1)
times <- seq(0, 100, by = 0.01)  # many time steps for plotting

system("R CMD SHLIB lorenzc.c")

dll <- paste0("lorenzc", .Platform$dynlib.ext)
dyn.load(dll)

## Comparison
out1 <- lsoda(state, times, func = "derivs", parms = parameters,
          dllname = "lorenzc", initfunc = "initmod")


## pre-identified symbols of DLL functions.
symbols <- checkDLL(func = "derivs", jacfunc = NULL, dllname = "lorenzc",
             initfunc = "initmod", verbose = TRUE, nout = 0,
             outnames = NULL, JT = 1)

## has now a new class deSolve.symbols
class(symbols)

out2 <- lsoda(state, times, func = symbols, parms = parameters,
      dllname = "lorenzc", initfunc = "initmod")

## dllname and initfunc not needed
out3 <- lsoda(state, times, func = symbols, parms = parameters)

plot(out1, out2, out3, mfrow=c(3,1))

all.equal(out1, out3)

## Benchmark ===================================================================
times <- seq(0, 1)    # one integration step for large number of runs
Nruns <- 5000         # or use 200 for large number of time steps

system.time(
  for(i in 1:Nruns)
    out1 <- lsoda(state, times, func = "derivs", parms = parameters,
      dllname = "lorenzc", initfunc = "initmod")
)


system.time(
  for(i in 1:Nruns) {
    out2 <- lsoda(state, times, func = symbols, parms = parameters)
  }
)

## setting hmax saves another 10%
system.time(
  for(i in 1:Nruns) {
    out3 <- lsoda(state, times, func = symbols, parms = parameters, hmax=1)
  }
)

plot(out1, out2, out3)


dyn.unload(dll)

