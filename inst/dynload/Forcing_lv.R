###############################################################################
# Implements the lv test model, as given in Forcing_lv.c
# A model in C-code and comprising a forcing function
# before trying this code, c program has to be compiled
# this can be done in R:
# system("R CMD SHLIB Forcing_lv.c")
# do make sure that these files are in the working directory...
# (if not, use setwd() )
###############################################################################

library(deSolve)

dyn.load(paste("Forcing_lv", .Platform$dynlib.ext, sep = ""))


#===============================================================================
# The R-code
#===============================================================================

SPCmod <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    import <- sigimp(t)
    dS <- import - b*S*P + g*C     # substrate
    dP <- c*S*P  - d*C*P           # producer
    dC <- e*P*C  - f*C             # consumer
    res <- c(dS, dP, dC)
    list(res,signal=import)
  })
}

## define states, time steps and parameters
init  <- c(S = 1, P = 1, C = 1)      # initial conditions
times  <- seq(0, 800, by=0.1)        # output times
parms  <- c(b = 0.1, c = 0.1, d = 0.1, e = 0.1, f = 0.1, g = 0.0)

## external input signal with rectangle impulse
signal <- as.data.frame(list(times = times,
                            import = rep(0, length(times))))
signal$import[signal$times >= 10 & signal$times <= 11] <- 0.2

signal$import <- ifelse((trunc(signal$times) %% 2 == 0), 0, 1)
ftime  <- seq(0, 900, 0.1)
sigimp <- approxfun(signal$times, signal$import, rule = 2)

Sigimp <- approx(signal$times, signal$import, xout=ftime ,rule = 2)$y
forcings <- cbind(ftime, Sigimp)

## Start values for steady state
xstart <- y <- c(S = 1, P = 1, C = 1)

## solve R version of the model
print(system.time(
  Out <- ode(xstart, times, SPCmod, parms))
)

## solve C version of the model
print(system.time(
  out <-  ode(y = y, times, func = "derivsc",
            parms = parms, dllname = "Forcing_lv", initforc="forcc",
            forcings = forcings, initfunc = "parmsc", nout = 2,
            outnames = c("Sum", "signal"))
))


## Plotting
out <- as.data.frame(out) # R version
#out <- as.data.frame(Out) # C version

mf <- par(mfrow = c(2,2))
plot(out$time, out$S,  type = "l", ylab = "substrate")
plot(out$time, out$P,  type = "l", ylab = "producer")
plot(out$time, out$C,  type = "l", ylab = "consumer")
plot(out$P, out$C,     type = "l", xlab = "producer", ylab = "consumer")
#points(Out$P,Out$C)

par(mfrow = mf)
tail(out)

dyn.unload(paste("Forcing_lv", .Platform$dynlib.ext, sep = ""))
