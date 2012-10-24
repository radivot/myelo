\name{epo}
\alias{epo}
\docType{data}
\title{Erythropoietin (Epo) Dynamic Response Data of Becker et al. Science 2010}
\description{
  This is the portion of the Becker et al. Science 2010 data that is used 
  in the profile likelihood analysis Raue et al. Chaos 2010.
}

\usage{data(epo)}
\format{
  A data frame with 3 variables.	
  \describe{
    \item{\code{epo.e}}{The counts per minute (cpm) of Epo external to the cell, degraded or intact, but free.}
    \item{\code{epo.m}}{The cpm of Epo bound to the Epo Receptor (EpoR) on the cell's plasma membrane.}
    \item{\code{epo.i}}{The cpm of Epo internal to the cell, either on an internal EpoR or degraded and free.}
   }
}

\details{
  The values involve a scale factor. 
  }

\source{
  The data were taken from the excel file in the zip download of "Epo receptor model" under 
  \bold{Applications} at 
  \url{http://www.fdmold.uni-freiburg.de/~araue/pmwiki.php/Main/Software}
}

\references{ Becker, V., Schilling, M., Bachmann, J., Baumann, U., Raue, A., Maiwald,
T., Timmer, J., and Klingmüller, U., Covering a broad dynamic range:
Information processing at the erythropoietin receptor, \emph{Science} \bold{328},
1404–1408, 2010. }
\examples{
library(myelo)
epo
with(epo,matplot(time,epo[,-1]))
}
\keyword{datasets}
