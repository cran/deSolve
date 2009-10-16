##------------------------------------------------------------------------------
##  This demo demonstrates different options and possibilities
##    of the Runge-Kutta-Solvers.
##  The external signal ("input") contains a discontinuity
##    to challenge stepsize algorithm and interpolation.
##  See help files for more options and examples.
##------------------------------------------------------------------------------

oask   <- devAskNewPage(dev.interactive(orNone = TRUE))
oldpar <- par(mfrow = c(2,2))

## a helper function for plotting
plotIt <- function(col = "red") {

  plot (out1$time, out1$S,  type = "l",   ylim = c(0,3))
  lines(out2$time, out2$S, col = col,   lty = "dotted", lwd = 3)
  lines(out3$time, out3$S, col = "green", lty = "dotted")

  plot (out1$time, out1$P, type = "l",    ylim = c(0,3))
  lines(out2$time, out2$P, col = col,   lty = "dotted", lwd = 3)
  lines(out3$time, out3$P, col = "green", lty = "dotted")

  plot (out1$time, out1$C, type = "l",    ylim = c(0,3))
  lines(out2$time, out2$C, col = col,   lty = "dotted", lwd = 3)
  lines(out3$time, out3$C, col = "green", lty = "dotted")

  plot (out1$P, out1$C, type = "l")
  lines(out2$P, out2$C, col = col,   lty = "dotted", lwd = 3)
  lines(out3$P, out3$C, col = "green", lty = "dotted")
}


## A simple resource limited Lotka-Volterra-Model
## The parameters
parms  <- c(b = 0.0, c = 0.1, d = 0.1, e = 0.1, f = 0.1, g = 0.0)

## The model
SPCmod <- function(t, x, parms, input)  {
  with(as.list(c(parms, x)),  {
    import <- input(t)
    dS <- import - b*S*P + g*C
    dP <- c*S*P  - d*C*P
    dC <- e*P*C  - f*C
    res<-c(dS, dP, dC)
    list(res)
  })
}

## vector of timesteps
times  <- seq(0, 100, length = 1001)

## external signal with rectangle impulse
signal <- as.data.frame(list(times = times,
                          import = rep(0,length(times))))

signal$import[signal$times > 10 & signal$times <= 11] <- 0.2

sigimp <- approxfun(signal$times, signal$import, rule = 2)


## Start values for steady state
xstart <- c(S = 1, P = 1, C = 1)


## Classical RK4 with fixed time step
system.time(out1  <- as.data.frame(rk4(xstart, times, SPCmod, parms, input = sigimp)))

## LSODA
system.time(out3 <- as.data.frame(lsoda(xstart, times, SPCmod,
  hmax = 1, parms, input = sigimp)))

## same: rk4 fixed step, generalized implementation
system.time(out2 <- as.data.frame(rk(xstart, times, SPCmod, parms, hini = 1,
  method = rkMethod("rk4"), input = sigimp)))

plotIt("blue")

## rk2, note smaller time step !!
system.time(out2 <- as.data.frame(rk(xstart, times, SPCmod, parms, hini = .5,
  method = rkMethod("rk2"), input = sigimp)))

plotIt()

## Euler, note smaller time step !!
system.time(out2 <- as.data.frame(rk(xstart, times, SPCmod, parms, hini = .1,
  method = rkMethod("euler"), input = sigimp)))

plotIt("blue")

# rk23bs
system.time(out2 <- as.data.frame(rk(xstart, times, SPCmod, parms, hini = 1,
  method = rkMethod("rk23bs"), input = sigimp)))

plotIt()

## Runge-Kutta-Fehlberg  45
## allows larger internal steps, but discontinuity disturbs
##  default 5th-order polynomial interpolation
system.time(out2 <- as.data.frame(rk(xstart, times, SPCmod, parms, hini = 1,
  hmax = 1, method = rkMethod("rk45f"), input = sigimp)))

plotIt("blue")

## Runge-Kutta-Fehlberg  45
##  3rd order polynomials are smoother
system.time(out2 <- as.data.frame(rk(xstart, times, SPCmod, parms, hini = 1,
  hmax = 1, method = rkMethod("rk45f", nknots=3), input = sigimp)))

plotIt()

## Prince-Dormand  5(4)7m
##  has built-in polynomial interpolation (dense output)
##  that is much better
system.time(out2 <- as.data.frame(rk(xstart, times, SPCmod, parms, hini = 1,
  hmax = 1, method = rkMethod("rk45dp7"), input = sigimp)))

plotIt("blue")

##other syntax
system.time(out2 <- as.data.frame(rk(xstart, times, SPCmod, parms, hini = 1,
  hmax = 1, method = "rk45dp7", input = sigimp)))

plotIt()

## tolerance values can also be set for single state variables
## if all tolerance values are zero, hmax is taken as fixed step
## Prince-Dormand  5(4)7m
system.time(out2 <- as.data.frame(rk(xstart, times, SPCmod, parms, hini = 1,
  hmax = 1, method = rkMethod("rk45dp7"), input = sigimp,
  atol = c(0, 0, 0), rtol = c(0, 0, 0), maxsteps = 1e5)))

plotIt("blue")

## clean up
par(oldpar)
devAskNewPage(oask)
