# require (shape)

#===============================================================================
# R-examples from CHAPTER 3
# chapter 3.1 - the basic lotka-volterra predator-prey model.
#===============================================================================

LVmod <- function(Time,State,Pars)
 {
   with(as.list(c(State,Pars)),
    {
    Ingestion    <- rIng  * Prey*Predator
    GrowthPrey   <- rGrow * Prey*(1-Prey/K)
    MortPredator <- rMort * Predator

    dPrey        <- GrowthPrey - Ingestion
    dPredator    <- Ingestion*assEff -MortPredator

    return(list(c( dPrey, dPredator)))
    })
 }


  nrun <- 1  # set to 10 for timing

  pars    <- c(rIng   =0.2,    # /day, rate of ingestion
               rGrow  =1.0,    # /day, growth rate of prey
               rMort  =0.2 ,   # /day, mortality rate of predator
               assEff =0.5,    # -, assimilation efficiency
               K      =10  )   # mmol/m3, carrying capacity

  yini    <- c(Prey=1,Predator=2)
  times   <- seq(0,200,by=1)

print(system.time(
 for (i in 1:nrun)
  out     <- as.data.frame(lsoda(func= LVmod, y=yini,
                           parms=pars, times=times))
  )/nrun)

print(system.time(
 for (i in 1:nrun)
  out     <- as.data.frame(lsode(func= LVmod, y=yini,
                           parms=pars, times=times))
  )/nrun)

print(system.time(
 for (i in 1:nrun)
  out     <- as.data.frame(vode(func= LVmod, y=yini,
                           parms=pars, times=times))
  )/nrun)

print(system.time(
 for (i in 1:nrun)
  out     <- as.data.frame(daspk(func= LVmod, y=yini,
                           parms=pars, times=times))
  )/nrun)

print(system.time(
 for (i in 1:nrun)
  out     <- as.data.frame(lsodes(func= LVmod, y=yini,
                           parms=pars, times=times))
  )/nrun)
  

  matplot(out$time,out[,2:3],type="l",xlab="time",ylab="Conc",
          main="Lotka-Volterra",lwd=2)
  legend("topright",c("prey", "predator"),col=1:2, lty=1:2)

#===============================================================================
# chapter 3.2 - predator-prey model with stopping criterium.
#===============================================================================

rootfun <- function(Time,State,Pars)
{
  dstate <- unlist(LVmod(Time,State,Pars))
  root   <- sum(abs(dstate))- 1e-4
}

print(system.time(
 for (i in 1:nrun)
out <- as.data.frame(lsodar(func=LVmod,y=yini,parms=pars,
                            times=times,rootfun=rootfun))
  )/nrun)
  matplot(out$time,out[,2:3],type="l",xlab="time",ylab="Conc",
          main="Lotka-Volterra with root",lwd=2)

#===============================================================================
# chapter 3.3 - predator-prey model in 1-D.
#===============================================================================
lvmod1D <- function (time, state, parms, N, Da, dx)

{
 with (as.list(parms),{

  Prey <- state[1:N]
  Pred <- state[(N+1):(2*N)]

# Dispersive fluxes; zero-gradient boundaries

  FluxPrey <- -Da * diff(c(Prey[1],Prey,Prey[N]))/dx
  FluxPred <- -Da * diff(c(Pred[1],Pred,Pred[N]))/dx

# Biology: Lotka-Volterra dynamics
  Ingestion     <- rIng * Prey *Pred
  GrowthPrey    <- rGrow* Prey *(1- Prey/K)
  MortPredator  <- rMort* Pred

# Rate of change = Flux gradient + Biology
  dPrey   <- -diff(FluxPrey)/dx  + GrowthPrey - Ingestion
  dPred   <- -diff(FluxPred)/dx  + Ingestion*assEff -MortPredator

  return (list(c(dPrey,dPred)))
  })
}

  R  <- 20                    # total length of surface, m
  N  <- 1000                  # number of boxes
  dx <- R/N                   # size of box in x-direction
  Da <- 0.05                  # m2/d, dispersion coefficient

  yini    <- rep(0,2*N)
  yini[500:501] <- yini[1500:1501] <- 10

  times  <-seq(0,200,by=1)   # output wanted at these time intervals

# based on lsode
print(system.time(
 for (i in 1:nrun)
  out    <- ode.1D(y=yini,times=times,func=lvmod1D,parms=pars,nspec=2,
                    N=N,dx=dx,Da=Da)
)/nrun)

print(system.time(
 for (i in 1:nrun)
  out    <- ode.1D(y=yini,times=times,func=lvmod1D,parms=pars,nspec=2,
                    N=N,dx=dx,Da=Da,method="vode")
)/nrun)

print(system.time(
 for (i in 1:nrun)
  out    <- ode.1D(y=yini,times=times,func=lvmod1D,parms=pars,nspec=2,
                    N=N,dx=dx,Da=Da,method="lsoda")
)/nrun)

print(system.time(
 for (i in 1:nrun)
  out    <- ode.1D(y=yini,times=times,func=lvmod1D,parms=pars,nspec=2,
                    N=N,dx=dx,Da=Da,method="lsodes")
)/nrun)


  PREY   <- out[,2:(N  +1)]

  filled.contour(x=times,z=PREY,color= topo.colors,
                 xlab="time, days", ylab= "Distance, m",
                 main="Prey concentration, 1-D")

#===============================================================================
# chapter 3.4 - predator-prey model in 2-D.
#===============================================================================

lvmod2D <- function (time, state, parms, N, Da, dx, dy)

{
  NN <- N*N
  Prey <- matrix(nr=N,nc=N,state[1:NN])
  Pred <- matrix(nr=N,nc=N,state[(NN+1):(2*NN)])

  with (as.list(parms),
  {
   dPrey   <- rGrow* Prey *(1- Prey/K) - rIng* Prey *Pred
   dPred   <- rIng* Prey *Pred*assEff -rMort* Pred

   zero <- rep(0,N)

   # Fluxes in x-direction; zero fluxes near boundaries
   FluxPrey <- -Da * rbind(zero,(Prey[2:N,]-Prey[1:(N-1),]),zero)/dx
   FluxPred <- -Da * rbind(zero,(Pred[2:N,]-Pred[1:(N-1),]),zero)/dx

   dPrey    <- dPrey - (FluxPrey[2:(N+1),]-FluxPrey[1:N,])/dx
   dPred    <- dPred - (FluxPred[2:(N+1),]-FluxPred[1:N,])/dx


   # Fluxes in y-direction
   FluxPrey <- -Da * cbind(zero,(Prey[,2:N]-Prey[,1:(N-1)]),zero)/dy
   FluxPred <- -Da * cbind(zero,(Pred[,2:N]-Pred[,1:(N-1)]),zero)/dy

   dPrey    <- dPrey - (FluxPrey[,2:(N+1)]-FluxPrey[,1:N])/dy
   dPred    <- dPred - (FluxPred[,2:(N+1)]-FluxPred[,1:N])/dy

  return (list(c(as.vector(dPrey),as.vector(dPred))))
 })}


  R  <- 20                    # total length of surface, m
  N  <- 50                    # number of boxes
  dx <- R/N                   # size of box in x-direction
  dy <- R/N                   # size of box in y-direction
  Da <- 0.05                  # m2/d, dispersion coefficient
  NN <- N*N
  
  yini     <- rep(0,2*N*N)
  cc       <- c((NN/2):(NN/2+1)+N/2,(NN/2):(NN/2+1)-N/2)
  yini[cc] <- yini[NN+cc] <- 10

  times  <-seq(0,200,by=1)   # output wanted at these time intervals

print(system.time(
 for (i in 1:nrun)

  out <- ode.2D(y=yini,times=times,func=lvmod2D,parms=pars,dimens=c(N,N),
                N=N,dx=dx,dy=dy,Da=Da,ynames=FALSE,lrw=440000)

                )/nrun)


 Col<- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
#                 topo.colors
 par(mfrow=c(2,2))
 par(oma=c(0,0,2,0))
 image(matrix(nr=N,nc=N,out[1,-1]), col=Col(100),main="initial",xlab="x",ylab="y")
 image(matrix(nr=N,nc=N,out[21,-1]), col=Col(100),main="20 days",xlab="x",ylab="y") 
 image(matrix(nr=N,nc=N,out[31,-1]), col=Col(100),main="30 days",xlab="x",ylab="y") 
 image(matrix(nr=N,nc=N,out[41,-1]), col=Col(100),main="40 days",xlab="x",ylab="y") 
 mtext(side=3,outer=TRUE,cex=1.25,"Lotka-Volterra Prey concentration on 2-D grid")
# filled.contour(matrix(nr=N,nc=N,out[20,-1]), color.palette=topo.colors,main="2-D grid") 
                
                
## DAE example
  Res_DAE <- function (t,y,yprime, pars, K)
      {
      with (as.list(c(y,yprime,pars)), {
  
        # residuals of lumped rates of changes
        res1 = -dD - dA + prod
        res2 = -dB + dA - r*B
        
        # and the equilibrium equation
        eq = K*D - A*B
  
        return(list(c(res1,res2,eq),
                    CONC=A+B+D))
  })
  }
  times <- seq(0,100,by=2)

  pars <- c(r     = 1,
            prod  = 0.1) 
  
  K <- 1

  # Initial conc; D is in equilibrium with A,B
  yini  <- c(A=2,B=3,D=2*3/K)
  
  # Initial rate of change
  dyini <- c(dA=0, dB=0, dD=0) 
  
  # DAE model solved with daspk
  DAE <- as.data.frame(daspk(y=yini,dy=dyini,times=times,
           res=Res_DAE,parms=pars,atol=1e-10,rtol=1e-10, K=1))
  
  #------------------------------------------------------
  # plotting output
  #------------------------------------------------------
  opa <- par(mfrow=c(2,2))
  for (i in 2:5) 
  {
  nm <- paste("[",names(DAE[i]),"]")
  if (i==5) nm <- "total conc"
  plot(DAE$time,DAE[,i],xlab="time", lwd=2,
       ylab="conc",main=nm,type="l")
#  writelabel(nr=i-1)
  }

mtext(outer=TRUE, side=3, "DAE chemical model",cex=1.25)

               