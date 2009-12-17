### ============================================================================
### Butcher tables for selected explicit ODE solvers of Runge-Kutta type
### Note that for fixed step methods A is a vector (the subdiagonal of matrix A)
###   For variable time step methods, A must be strictly lower triangular.
###   The underlying rk code is currently restricted to explicit methods.
### ============================================================================

rkMethod <- function(method = NULL, ...) {
  methods <- list(
    euler = list(ID = "euler",
        varstep = FALSE,
          A      = c(0),
          b1     = c(1),
          c      = c(0),
          stage  = 1,
          Qerr   = 1
    ),
    ## Heun's method
    rk2 = list(ID = "rk2",
        varstep = FALSE,
          A      = c(0, 1),
          b1     = c(0.5, 0.5),
          c      = c(0, 1),
          stage  = 2,
          Qerr   = 1
    ),
    ## classical Runge-Kutta 4th order method
    rk4 = list(ID = "rk4",
        varstep = FALSE,
          A      = c(0, .5, .5, 1),
          b1     = c(1/6, 1/3, 1/3, 1/6),
          c      = c(0, .5, .5, 1),
          stage  = 4,
          Qerr   = 4
    ),
    ## One of the numerous RK23 formulae
    rk23 = list(ID = "rk23",
      varstep = TRUE,
      FSAL    = FALSE,
      A  = matrix(c(0, 0, 0,
                  1/2, 0, 0,
                  -1, 2, 0), 3, 3, byrow = TRUE),
      b1 = c(0, 1, 0),
      b2 = c(1/6, 2/3, 1/6),
      c  = c(0, 1/2, 2),
      stage = 3,
      Qerr  = 2
    ),
    ## Bogacki & Shampine
    rk23bs = list(ID = "rk23bs",
      varstep = TRUE,
      FSAL    = TRUE,
      A  = matrix(c(0, 0, 0, 0,
                  1/2, 0, 0, 0,
                  0, 3/4, 0, 0,
                  2/9, 1/3, 4/9, 0), 4, 4, byrow = TRUE),
      b1 = c(2/9, 1/3, 4/9, 0),
      b2 = c(7/24, 1/4, 1/3, 1/8),
      c  = c(0, 1/2, 3/4, 1),
      stage = 4,
      Qerr  = 2
    ),
    ## RK-Fehlberg 34
    rk34f = list(ID = "rk34f",
         varstep = TRUE,
         A  = matrix(c(0, 0, 0, 0,
                      2/7, 0, 0, 0,
                      77/900, 343/900, 0, 0,
                      805/1444, -77175/54872, 97125/54872, 0,
                      79/490, 0, 2175/3626, 2166/9065),
                      5, 4, byrow = TRUE),
         b1 = c(79/490, 	0, 2175/3626, 2166/9065, 	0),
         b2 = c(229/1470, 0, 1125/1813, 13718/81585, 1/18),
         c  = c(0,	2/7, 	7/15, 35/38, 	1),
         stage = 5,
         Qerr  = 3
    ),
    ## RK-Fehlberg Method 45
    rk45f = list(ID = "rk45f",
         varstep = TRUE,
         A  = matrix(c(0, 0, 0, 0, 0,
                      1/4, 0, 0, 0, 0,
                      3/32, 9/32, 0, 0, 0,
                      1932/2197, -7200/2197, 7296/2197, 0, 0,
                      439/216, -8, 3680/513, -845/4104, 0,
                      -8/27, 2, -3544/2565, 1859/4104, -11/40),
                      6, 5, byrow = TRUE),
         b1 = c(16/135, 	0, 	6656/12825, 	28561/56430, 	-9/50, 	2/55),
         b2 = c(25/216, 	0, 	1408/2565, 	2197/4104, 	-1/5, 	0),
         c  = c(0,	1/4, 	3/8, 	12/13, 	1, 	1/2),
         stage = 6,
         Qerr  = 4
    ),
    ## Cash-Karp method
    rk45ck = list(ID = "rk45ck",
            varstep = TRUE,
            A = matrix(c(0,    0,       0,         0,            0,
                         1/5,  0,       0,         0,            0,
                         3/40, 9/40,    0,         0,            0,
                         3/10, -9/10,   6/5,       0,            0,
                       -11/54, 5/2,    -70/27,     35/27,        0,
                   1631/55296, 175/512, 575/13824, 44275/110592, 253/4096),
                         6, 5, byrow = TRUE),
            b1 = c(37/378, 0, 250/621, 125/594, 0, 512/1771),
            b2 = c(2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 1/4),
            c = c(0, 1/5, 3/10, 3/5,  1, 7/8),
            stage = 6,
            Qerr = 4),
    ## England Method
    rk45e = list(ID = "rk45e",
         varstep = TRUE,
         A  = matrix(c(0, 0, 0, 0, 0,
                      1/2, 0, 0, 0, 0,
                      1/4, 1/4, 0, 0, 0,
                      0, -1, 2, 0, 0,
                      7/27, 10/27, 0, 1/27, 0,
                      28/625, -125/625, 546/625, 54/625, -378/625),
                      6, 5, byrow = TRUE),
         b1 = c(1/6, 	0, 4/6, 1/6, 	0, 	0),
         b2 = c(14/336, 0, 0,	35/336, 162/336, 125/336),
         c  = c(0,	1/2, 	1/2, 	1, 	2/3, 	1/5),
         stage = 6,
         Qerr  = 4
    ),
    ## Prince-Dormand 5(4)6m
    rk45dp6 = list(ID = "rk45dp6",
         varstep = TRUE,
         A  = matrix(c(0, 0, 0, 0, 0,
                      1/5, 0, 0, 0, 0,
                      3/40, 9/40, 0, 0, 0,
                      3/10, -9/10, 6/5, 0, 0,
                      226/729, -25/27, 880/729, 55/729, 0,
                      -181/270, 5/2, -266/297, -91/27, 189/55),
                      6, 5, byrow = TRUE),
         b1 = c(31/540, 	0, 190/297, -145/108, 351/220, 1/20),
         b2 = c(19/216, 0, 1000/2079,	-125/216, 81/88, 5/56),
         c  = c(0,	1/5, 	3/10, 3/5, 	2/3, 	1),
         stage = 6,
         Qerr  = 4
    ),
    ## Prince-Dormand 5(4)7m -- recommended by the Octave developers
    rk45dp7 = list(ID = "rk45dp7",
         varstep = TRUE,
         FSAL    = TRUE,
         A  = matrix(c(0, 0, 0, 0, 0, 0,
                      1/5, 0, 0, 0, 0, 0,
                      3/40, 9/40, 0, 0, 0, 0,
                      44/45, -56/15, 32/9, 0, 0, 0,
                      19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,
                      9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,
                      35/384, 0, 500/1113, 125/192, -2187/6784, 11/84),
                      7, 6, byrow = TRUE),
         b1 = c(5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40),
         b2 = c(35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0),
         c  = c(0, 1/5, 3/10, 4/5, 8/9, 1, 1),
         d  = c(-12715105075.0/11282082432.0, 0, 87487479700.0/32700410799.0,
                -10690763975.0/1880347072.0, 701980252875.0/199316789632.0,
                -1453857185.0/822651844.0, 69997945.0/29380423.0),
         stage = 7,
         Qerr  = 4
    ),
    ## Runge-Kutta-Fehlberg 78 method
    rk78f = list(ID = "rk78f",
        varstep = TRUE,
        FSAL    = FALSE,
        A  = matrix(
         c(rep(0,12),
         2/27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1/36, 1/12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1/24, 0, 1/8, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         5/12, 0, -25/16, 25/16, 0, 0, 0, 0, 0, 0, 0, 0,
         0.05, 0, 0, 0.25, 0.2, 0, 0, 0, 0, 0, 0, 0,
         -25/108, 0, 0, 125/108, -65/27, 125/54, 0, 0, 0, 0, 0, 0,
         31/300, 0, 0, 0, 61/225, -2/9, 13/900, 0, 0, 0, 0, 0,
         2, 0, 0, -53/6, 704/45, -107/9, 67/90, 3, 0, 0, 0, 0,
         -91/108, 0, 0, 23/108, -976/135, 311/54, -19/60, 17/6, -1/12, 0, 0, 0,
         2383/4100, 0, 0, -341/164, 4496/1025, -301/82, 2133/4100, 45/82, 45/164, 18/41, 0, 0,
         3/205, 0, 0, 0, 0, -6/41, -3/205, -3/41, 3/41, 6/41, 0, 0,
         -1777/4100, 0, 0, -341/164, 4496/1025, -289/82, 2193/4100, 51/82, 33/164, 12/41, 0, 1
        ), nrow=13, ncol=12, byrow = TRUE),
        b1 = c(41/840, 0,0,0,0, 34/105, 9/35, 9/35, 9/280, 9/280, 41/840, 0, 0),
        b2 = c(0, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 0, 41/840, 41/840),
        c  = c(0, 2./27., 1/9, 1/6, 5/12, 0.5, 5/6, 1/6, 2/3, 1/3, 1, 0, 1),
        stage = 7,
        Qerr  = 8
    )
  )
  ## look if the method is known; ode23 and ode45 are used as synonyms
  knownMethods <- c(lapply(methods,"[[", "ID"), "ode23", "ode45")

  if (!is.null(method)) {
    method <- unlist(match.arg(method, knownMethods))
    if (method == "ode23")
      method <- "rk23bs"
    else if (method == "ode45")
      method <- "rk45dp7"

    out <- methods[[method]]
  } else {
    out <- vector("list", 0)
  }

  ## modify a known or add a completely new method)
  ldots <- list(...)
  out[names(ldots)] <- ldots

  ## return the IDs of the methods if called with an empty argument list
  if (is.null(method) & length(ldots) == 0) {
    out <- as.vector(unlist(knownMethods))
  } else {
    ## check size consistency of parameter sets
    sl    <- lapply(out, length)
    stage <- out$stage
    if (is.matrix(out$A)) {
      if (nrow(out$A) != stage | ncol(out$A)  < stage -1 | ncol(out$A) > stage)
        stop("Size of matrix A does not match stage")
    } else {
      if (length(out$A) != stage) stop("Size of A does not match stage")
    }
    if (stage != sl$b1 | stage != sl$c)
      stop("Wrong rkMethod, length of parameters do not match")
    if (out$varstep & is.null(out$b2))
      stop("Variable stepsize method needs non-empty b2")
    if (!is.null(out$b2))
      if (sl$b2 != stage)
        stop("Wrong rkMethod, length of b2 must be empty or equal to stage")
    if (!is.null(out$d))
      if (sl$d != stage)
        stop("Wrong rkMethod, length of d must be empty or equal to stage")
    class(out) <- c("list", "rkMethod")
  }

  out
}
