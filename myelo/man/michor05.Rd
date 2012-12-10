\name{michor05}
\alias{michor05}
\title{CML model of Michor et al. 2005}
\usage{
  michor05(Time, State, Pars)
}
\arguments{
  \item{Time}{The parameters of this model do not depend on
  time: this argument exists here as a dummy only because
  deSolve expects it to exist.}

  \item{State}{Vector of current states. The elements of
  this vector must be named because within this function's
  definition is \code{with(as.list(State,Pars), {XXX} )}
  and code in XXX presupposes that it can call named
  elements as variables with those names.}

  \item{Pars}{Vector of parameter values with \emph{named}
  (see above) elements.}
}
\value{
  A list of length 2 (as expected by deSolve) where the
  first list element is the vector of derivatives, i.e. the
  right hand side of the ODEs, and the second element of
  the list is a vector of auxiliary variables that one may
  wish to track over time.
}
\description{
  This function returns the right hand side of an ordinary
  differential equation model of CML that was published by
  Michor et al. in 2005. The intended use of this function
  is as an argument to \code{ode()} of the \code{deSolve}
  package.
}
\note{
  This work was supported by the National Cancer Institute
  and Tufts Integrative Cancer Biology Program under
  U54CA149233-029689.
}
\examples{
\dontrun{
library(myelo)
pars=c(d0=0.003,d1=0.008,d2=0.05,d3=1,ax=0.8,bx=5,cx=100,ry=0.008,
ay=1.6, byy=10, cy=100, rz=0.023, az=1.6, bz=10, cz=100,trt=0,u=0,X0o=2e6)
(y0<-c(X0=2e6,X1=0,X2=0,X3=0,Y0=0,Y1=0,Y2=0,Y3=0,Z0=0,Z1=0,Z2=0,Z3=0))
out=ode(y=y0,times=seq(0,1000,1),michor05, parms=pars, rtol=1e-4, atol= rep(1e-4,12))
plot(out)
}
}
\author{
  Tom Radivoyevitch (\email{txr24@case.edu})
}
\references{
  Franziska Michor, Timothy P. Hughes, Yoh Iwasa, Susan
  Branford, Neil P. Shah, Charles L. Sawyers & Martin A.
  Nowak Dynamics of chronic myeloid leukaemia,
  \emph{Nature}, \bold{435}, 1267-1270 (2005)
}
\seealso{
  \code{\link{myelo-package}, \link{rcj12}}
}
\keyword{IO}

