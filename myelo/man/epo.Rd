\name{epo}
\alias{epo}
\docType{data}
\title{Erythropoietin Data}
\description{
  This is the Epo data of the Epo Model of Raue et al. 
}
\usage{data(epo)}
\format{
  A data frame with 3 variables.
  \describe{
    \item{\code{epo.e}}{The counts per minute (cpm) of epo external to the cell, degraded or intact, but free.}
    \item{\code{epo.m}}{The cpm of epo bound to epoR on the cell's plasma membrane.}
    \item{\code{epo.i}}{The cpm of epo internal to the cell, either on an internal epoR or degraded.}
   }
}

\details{
  The values involve a scale factor. 
  }
\source{
   An excel file in the downloaded zip obtained under "epo receptor model" 
  \url{http://www.fdmold.uni-freiburg.de/~araue/pmwiki.php/Main/Software}
 
}
\references{ Becker, V., Schilling, M., Bachmann, J., Baumann, U., Raue, A., Maiwald,
T., Timmer, J., and Klingmüller, U., Covering a broad dynamic range:
Information processing at the erythropoietin receptor, \emph{Science} \bold{328},
1404–1408, 2010. }
\examples{
library(myelo)
epo
with(epo,matplot(time,epo[,-1])
}
\keyword{datasets}
