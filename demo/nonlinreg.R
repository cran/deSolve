### Another nonlinear regression example (from Jim Lindsey)

if (require(gnlm, quietly=TRUE)) {
  t <- 1:10

  ## known initial conditions
  ## explicit regression formula with gamma variability
  fn1 <- ~exp(-a*t)
  y <- rgamma(10, 2, finterp(fn1)(0.1)/2)
  print(gnlr(y, dist="gamma", mu=fn1, pmu=0.2, pshape=0))

  ## regression formula defined by differential equation
  dfn <- function(t, y, p) list(-p*y)
  fn2 <- ~lsoda(1, c(0,t), dfn, p)[2:11,2]
  print(finterp(fn1)(0.1))
  print(finterp(fn2)(0.1))
  print(gnlr(y, dist="gamma", mu=fn2, pmu=0.2, pshape=0))

  ## estimated initial conditions
  ## explicit regression formula with gamma variability
  fn3 <- ~b*exp(-a*t)
  y <- rgamma(10, 2, finterp(fn3)(c(2,0.1))/2)
  print(gnlr(y, dist="gamma", mu=fn3, pmu=c(1,0.2), pshape=0))

  ## regression formula defined by differential equation
  fn4 <- ~lsoda(b, c(0,t), dfn, a)[2:11,2]
  print(finterp(fn3)(c(1,0.1)))
  print(finterp(fn4)(c(1,0.1)))
  print(gnlr(y, dist="gamma", mu=fn4, pmu=c(1,0.2), pshape=0))
} else {
  cat("This example requires Jim Lindsay's gnlm package\n")
}
