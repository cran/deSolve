\name{aquaphy}
\alias{aquaphy}
\title{a physiological model of unbalanced algal growth}
\description{A phytoplankton model with uncoupled carbon and nitrogen assimilation as a function of light and Dissolved Inorganic Nitrogen (DIN) concentration.
  \cr Algal biomass is described via 3 different state variables:
  \item low molecular weight carbohydrates (LMW), the product of photosynthesis,
  \item storage molecules (RESERVE) and
  \item the biosynthetic and photosynthetic apparatus (PROTEINS).
  
  All algal state variables are expressed in mmol C m-3.
  \cr Only proteins contain nitrogen and chlorophyll, with a fixed stoichiometric ratio.
  \cr As the relative amount of proteins changes in the algae, so does the N:C and the Chl:C ratio.
  
  An additional state variable, dissolved inorganic nitrogen (DIN) has units of mmol N m-3.
  
  The algae grow in a dilution culture (chemostat): there is constant inflow of DIN and outflow of culture water, including DIN and algae, at the same rate. 
  
  There is a day-night illumination regime, i.e. the light is switched on and off at fixed times (where the sum of illuminated+dark period=24 hours).}
\usage{aquaphy(times,y,parms,...)}
\arguments{
  \item{times }{time sequence for which output is wanted; the first value of times must be the initial time}
  \item{y }{the initial (state) values ("DIN","PROTEIN","RESERVE","LMW"), in that order.}
  \item{parms }{vector or list with the aquaphy model parameters; see the example for the order in which these have to be defined}
  \item{... }{ any other parameters passed to the integrator ode (which solves the model)}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{#-----------------------#
  # the model parameters: #
  #-----------------------#
  
  parameters<-c(maxPhotoSynt   =0.125,      #molC/molC/hr
                rMortPHY       =0.001,      #/hr
                alpha          =-0.125/150, #uEinst/m2/s/hr
                pExudation     =0.0,        #-
                maxProteinSynt =0.136,      #molC/molC/hr
                ksDIN          =1.0,        #mmolN/m3
                minpLMW        =0.05,       #molC/molC
                maxpLMW        =0.15,       #molC/molC
                minQuotum      =0.075,      #molC/molC
                maxStorage     =0.23,       #/h
                respirationRate=0.0001,     #/h
                pResp          =0.4,        #-
                catabolismRate =0.06,       #/h
                dilutionRate   =0.01,       #/h
                rNCProtein     =0.2,        #molN/molC
                inputDIN       =10.0,       #mmolN/m3
                rChlN          =1,          #gChl/molN
                parMean        =250.,       #umolPhot/m2/s
                dayLength      =15.         #hours
                )
  
  #-------------------------#
  # the initial conditions: #
  #-------------------------#
  
  state     <-c(DIN     =6.,     #mmolN/m3
                PROTEIN =20.0,   #mmolC/m3
                RESERVE =5.0,    #mmolC/m3
                LMW     =1.0)    #mmolC/m3
  
  #----------------------#
  # RUNNING the model:   #
  #----------------------#
  
  times <-seq(0,24*20,1)
  
  out <- as.data.frame(aquaphy(times,state,parameters))
  
  #------------------------#
  # PLOTTING model output: #
  #------------------------#
  
  par(mfrow=c(2,2), oma=c(0,0,3,0))    
  col <- grey(0.9)
  ii <- 1:length(out$PAR)              
  
  plot (times[ii],out$Chlorophyll[ii],type="l",
        main="Chlorophyll",xlab="time, hours",ylab="ug/l")
  polygon(times[ii],out$PAR[ii]-10,col=col,border=NA);box()
  lines (times[ii],out$Chlorophyll[ii]  ,lwd=2 )
  
  
  plot (times[ii],out$DIN[ii]        ,type="l",main="DIN",
        xlab="time, hours",ylab="mmolN/m3")
  polygon(times[ii],out$PAR[ii]-10,col=col,border=NA);box()
  lines (times[ii],out$DIN[ii]  ,lwd=2 )
  
  
  plot (times[ii],out$NCratio[ii]    ,type="n",main="NCratio",
        xlab="time, hours",ylab="molN/molC")
  polygon(times[ii],out$PAR[ii]-10,col=col,border=NA);box()
  lines (times[ii],out$NCratio[ii]  ,lwd=2 )
  
  
  plot (times[ii],out$PhotoSynthesis[ii],type="l",
         main="PhotoSynthesis",xlab="time, hours",
         ylab="mmolC/m3/hr")
  polygon(times[ii],out$PAR[ii]-10,col=col,border=NA);box()
  lines (times[ii],out$PhotoSynthesis[ii]  ,lwd=2 )
  
  mtext(outer=TRUE,side=3,"AQUAPHY",cex=1.5)
  
  #------------------------#
  # SUMMARY  model output: #
  #------------------------#
  t(summary(out))}
\references{Lancelot, C., Veth, C. and Mathot, S. (1991). 
  Modelling ice-edge phytoplankton bloom in the Scotia-Weddel sea sector of the Southern Ocean during spring 1988. 
  Journal of Marine Systems 2, 333..346.
  
  Soetaert, K. and Herman, P. (2008). A practical guide to ecological modelling. 
  Using R as a simulation platform. Springer.
}
\details{The model is implemented primarily to demonstrate the linking of FORTRAN with R-code.
  \cr The source can be found in the \file{dynload} subdirectory of the package.}
\seealso{\code{\link{ccl4model}}, the CCl4 inhalation model}
\keyword{models}

