

#'Erythropoietin (Epo) Dynamic Response Data of Becker et al. Science 2010
#'
#'This is the portion of the Becker et al. Science 2010 data that is used in
#'the profile likelihood analysis Raue et al. Chaos 2010.
#'
#'The values involve a scale factor.
#'
#'@name epo
#'@docType data
#'@format A data frame that includes time and 3 variables.  \describe{
#'\item{time}{The time point in minutes after Epo addition.}
#'\item{epo.e}{The counts per minute (cpm) of Epo external to the cell,
#'degraded or intact, but free, see \link{raue10}.} 
#'\item{epo.m}{The
#'cpm of Epo bound to the Epo Receptor (EpoR) on the cell's plasma membrane.}
#'\item{epo.i}{The cpm of Epo internal to the cell, either on an
#'internal EpoR or degraded and free.} }
#'@note Coding supported by the National Cancer Institute and Tufts
#'Integrative Cancer Biology Program U54CA149233-029689.
#'@seealso \code{\link{myelo-package}, \link{rcj12}, \link{fokas91},
#'\link{raue10} }
#'@references Becker, V., Schilling, M., Bachmann, J., Baumann, U., Raue, A.,
#'Maiwald, T., Timmer, J., and Klingmuller, U., Covering a broad dynamic range:
#'Information processing at the erythropoietin receptor, \emph{Science}
#'\bold{328}, 1404-1408, 2010.
#'@source The data were taken from the excel file in the zip download of Epo
#'receptor model under \bold{Applications} at
#'\url{http://www.fdmold.uni-freiburg.de/~araue/pmwiki.php/Main/Software}
#'@keywords datasets
#'@examples
#'
#'library(myelo)
#'epo
#'with(epo,matplot(time,epo[,-1]))
#'
NULL





#'Models of myeloid hematopoiesis
#'
#'This package contains functions that compute the right hand sides of ODEs of
#'models relevant to myeloid hematopoiesis.
#'
#'\tabular{ll}{ Package: \tab \pkg{myelo}\cr Type: \tab Package\cr Depends:
#'\tab deSolve\cr Suggests: \tab bbmle \cr License: \tab GPL-2\cr LazyLoad:
#'\tab yes\cr LazyData: \tab yes\cr URL: \tab
#'\url{http://epbi-radivot.cwru.edu/myelo/myelo.html}\cr }
#'
#'@name myelo-package
#'@docType package
#'@note This work was supported by the National Cancer Institute and Tufts
#'Integrative Cancer Biology Program under U54CA149233-029689.
#'@author Tom Radivoyevitch (\email{txr24@@case.edu})
#'@seealso \code{\link{fokas91},\link{rcj12}}
#'@keywords package
NULL



