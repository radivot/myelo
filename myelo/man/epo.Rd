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
  A data frame with time in minutes and 3 variables.	
  \describe{
    \item{\code{time}}{The time point in minutes after Epo addition.}
    \item{\code{epo.e}}{The counts per minute (cpm) of Epo external to the cell, degraded or intact, but free, see \link{raue10}.}
    \item{\code{epo.m}}{The cpm of Epo bound to the Epo Receptor (EpoR) on the cell's plasma membrane.}
    \item{\code{epo.i}}{The cpm of Epo internal to the cell, either on an internal EpoR or degraded and free.}
   }
}

\details{
  The values involve a scale factor. 
  }

\source{
  The data were taken from the excel file in the zip download of Epo receptor model under 
  \bold{Applications} at 
  \url{http://www.fdmold.uni-freiburg.de/~araue/pmwiki.php/Main/Software}
}

\note{  This coding was supported by the National Cancer Institute and Tufts Integrative 
        Cancer Biology Program under U54CA149233-029689. }

\seealso{\code{\link{myelo-package}, \link{rcj12}, \link{fokas91}, \link{raue10} }}


\references{ Becker, V., Schilling, M., Bachmann, J., Baumann, U., Raue, A., Maiwald,
T., Timmer, J., and Klingmuller, U., Covering a broad dynamic range:
Information processing at the erythropoietin receptor, \emph{Science} \bold{328},
1404-1408, 2010. }
\examples{
library(myelo)
epo
with(epo,matplot(time,epo[,-1]))
}
\keyword{datasets}
