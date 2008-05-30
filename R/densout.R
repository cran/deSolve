## polynomial interpolation for dense output
denspar <- function(FF, y0, y1, dt, stage, d) {
  r <- matrix(0, nrow = length(y0), ncol = 5)
  r[,1] <- y0
  ydiff <- y1 - y0
  r[,2] <- ydiff
  bspl  <- dt * FF[, 1] - ydiff
  r[,3] <- bspl
  r[,4] <- ydiff - dt * FF[, stage] - bspl
  for (i in 1:stage) r[,5]  <- r[,5] + d[i] * FF[,i]
  #r[,5] <- colSums(t(FF) * d) # this is slower
  r[,5] <- r[,5] * dt
  return(r)
}

## is there a simpler and faster way??
densout <- function(r, t0, t, dt) {
   ntimes  <- length(t)
   nstates <- nrow(r)
   t   <- rep(t, each = nstates)
   s   <- (t - t0) /dt
   s1  <- 1 - s
   res <- r[,1] + s * (r[,2] + s1 * (r[,3] + s * (r[,4] + s1 * r[,5])))
   matrix(res, nrow = ntimes, ncol = nstates, byrow = TRUE)
}
