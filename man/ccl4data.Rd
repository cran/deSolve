\name{ccl4data}
\docType{data}
\alias{ccl4data}
\title{Closed Chamber Study of CCl4 Metabolism by Rats.}
\description{The results of a closed chamber experiment to determine
  metabolic parameters for CCl4 (carbon tetrachloride) in rats.
}
\usage{data(ccl4data)}
\format{This data frame contains the following columns:
  \describe{
    \item{time}{the time (in hours after starting the experiment).}
    \item{initconc}{initial chamber concentration (ppm).}
    \item{animal}{this is a repeated measures design; this variable
      indicates which animal the observation pertains to.
    }
    \item{ChamberConc}{chamber concentration at \code{time}, in ppm.}
  }
}
\source{
  Evans, et al. 1994 Applications of sensitivity analysis to a
  physiologically based pharmacokinetic model for carbon tetrachloride
  in rats.  Toxicology and Applied Pharmacology \bold{128}: 36 -- 44.
}
\examples{
plot(ChamberConc ~ time, data = ccl4data, xlab = "Time (hours)",
  xlim = range(c(0, ccl4data$time)),
  ylab = "Chamber Concentration (ppm)", log = "y")
ccl4data.avg <- aggregate(ccl4data$ChamberConc,
                  by = ccl4data[c("time", "initconc")], mean)
points(x ~ time, data = ccl4data.avg, pch = 16)                        
}
\keyword{datasets}

