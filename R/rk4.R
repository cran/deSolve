rk4 <- function(y, times, func, parms, hini=min(diff(times)), ...) {
  rk(y, times, func, parms, hmin = 0, hini = hini, 
    method = rkMethod("rk4", ...), maxsteps = Inf, ...)
}

euler <- function(y, times, func, parms, hini=min(diff(times)), ...) {
  rk(y, times, func, parms, hmin = 0, hini = hini, 
    method = rkMethod("euler", ...), maxsteps = Inf, ...)
}
