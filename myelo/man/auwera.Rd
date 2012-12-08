\docType{data}
\name{auwera}
\alias{auwera}
\title{G-CSF Dynamic Response Data of Van Der Auwera 2001}
\format{A data frame that includes time in hours and 8 variables.  \describe{
\item{hours}{Time points in hours after subcutanteous filgrastim injection.}
\item{gcsf5}{GCSF time course after a 5 ug/kg injection.}
\item{gcsf10}{GCSF time course after a 10 ug/kg injection.}
\item{cd34p5}{The number of CD34+ cells per cubic mm after a 5 ug/kg injection.}
\item{cd34p10}{The number of CD34+ cells per cubic mm after a 10 ug/kg injection.}
\item{cfugm5}{The number of CFUGM per 10^5 cells after 7 daily 5 ug/kg injections.}
\item{cfugm10}{The number of CFUGM per 10^5 cells after 7 daily 10 ug/kg injections.}
\item{anc5}{Absolute Neutrophil Count (ANC) time course after a 5 ug/kg injection.}
\item{anc10}{ANC time course after a 10 ug/kg injection.} }}
\source{
  The data were obtained from Figures 1,2,3, and 5 of
  vander Auwera et al. by first copying into powerpoint,
  then saving as a gif, and then using PlotDigitizer.
}
\description{
  In addition to the single injection data used in Scholz
  2012 this data.frame includes CFU-GM and CD34 cells
  counts in response to 7 daily filgrastim injections.
}
\note{
  Coding supported by the National Cancer Institute and
  Tufts Integrative Cancer Biology Program
  U54CA149233-029689.
}
\examples{
library(myelo)
auwera
par(mfrow=c(2,2),oma=c(0,0,2,0))
with(auwera,{matplot(hours,auwera[,2:3],ylab="GCSF")
matplot(hours,auwera[,4:5],ylab="CD34")
matplot(hours,auwera[,6:7],ylab="CFU-GM")
matplot(hours,auwera[,8:9],ylab="ANC")}
)
title("1 = 5 ug/kg           2 = 10 ug/kg",outer=TRUE)
}
\references{
  Phillippe van der Auwera, Erich Platzer, Z.-X. Xu, Rainer
  Schulz, Olivier Feugeas, Renaud Capdeville, and David J.
  Edwards, Pharmacodynamics and Pharmacokinetics of Single
  Doses of Subcutaneous Pegylated Human G-CSF Mutant (Ro
  25-8315) in Healthy Volunteers: Comparison With Single
  and Multiple Daily Doses of Filgrastim, \emph{American
  Journal of Hematology} \bold{66}, 245?251, 2001.
}
\seealso{
  \code{\link{myelo-package}, \link{rcj12}, \link{fokas91},
  \link{raue10} }
}
\keyword{datasets}

