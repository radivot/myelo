\name{rs13}
\alias{rs13}
\title{Radivoyevitch and Sachs 2013 model of
Auwera et al. 2001 data}
\usage{
  rs13(Time, State, Pars)
}
\arguments{
  \item{Time}{The current time.}

  \item{State}{Vector of current states. The elements of
  this vector must be named because within this function's
  definition is \code{with(as.list(State,Pars), {XXX} )}
  and code in XXX presupposes that it can call named
  elements as variables with those names.}

  \item{Pars}{Vector of parameter values with \emph{named}
  (see above) elements.}
}
\value{
  A list of length 1 where the first list element is the
  vector of derivatives, i.e. the right hand side of the
  ODEs
}
\description{
  This function returns the right hand side of a very
  simple ordinary differential equation model of neutrophil
  dynamics after G-CSF. The intended use of this function
  is as an argument to \code{ode()} of the \code{deSolve}
  package.
}
\details{
  In this model CD34 = CD34+ cells, Ni = immature cells in
  marrow that have finished dividing, Nm = mature
  neutrophils in marrow, and Nb = neutrophils in blood. If
  G-CSF mobilizes CD34, it likely mobilizes Ni as well, but
  since Auwere does not have a separate read-out of these
  in the blood, we will assume that the numbers in the
  blood are negligible relative to the numbers of mature
  cells pulled from the marrow reserve.
}
\note{
  This work was supported by the National Cancer Institute
  and Tufts Integrative Cancer Biology Program under
  U54CA149233-029689.
}
\examples{
\dontrun{
library(myelo)
# This is all still under construction ... don't expect it to work.
times <- c(0:200)
yout <- ode(c(Ni=1e12,Nm=1e11,CD34b=0,Nb=5000,Gt=0,Gb=0),
		times = times, func = rs13,parms = c(CD34=1e8,kgain0=1e4/240, # 240 hours = 10 days to amplify
      tau=120, # 5 days to mature
      kmobN0=2.5e10/24,   # 5000/uL = 25e9/5L and assume this many every 24 hours
      kmobNSlope=0,kgainSlope=0, # to be estimated from G-CSF data
      kt=0,kb=0,   #          to be estimated from G-CSF data
      Vb=5e6, # flux of mature cell numbers in marrow drop into 5 liters (5e6 uL) of blood as they mobilize
par(mfrow=c(1,1),mar=c(4.7,6,3.3,1),oma=c(0,0,0,0),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
plot(times,yout[,3]/1e8,type="l",xlab="days",ylab="Neutrophils")


}
}
\author{
  Tom Radivoyevitch (\email{txr24@case.edu})
}
\references{
  Phillippe van der Auwera, Erich Platzer, Z.-X. Xu,1
  Rainer Schulz, Olivier Feugeas, Renaud Capdeville, and
  David J. Edwards, Pharmacodynamics and Pharmacokinetics
  of Single Doses of Subcutaneous Pegylated Human G-CSF
  Mutant (Ro 25-8315) in Healthy Volunteers: Comparison
  With Single and Multiple Daily Doses of Filgrastim
  \emph{American Journal of Hematology} \bold{66} 245?251
  (2001).
}
\seealso{
  \code{\link{myelo-package}, \link{brooks12},
  \link{zhuge12}}
}
\keyword{IO}

