\name{scholz12}
\alias{scholz12}
\title{PKPD G-CSF model of Scholz et al. TBMM 2012}
\usage{
  scholz12(Time, State, Pars)
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
  right hand side of the ODEs and the second element of the
  list is a vector of auxiliary variables that one may wish
  to track over time.
}
\description{
  This function returns the right hand side of an ordinary
  differential equation model of G-CSF published by M.
  Scholz, S. Schirm, M. Wetzler, C. Engel and M. Loeffler
  in Theoretical Biology and Medical Modelling 2012, 9:32.
  The intended use of this function is as an argument to
  \code{ode()} of the \code{deSolve} package.
}
\details{
  This model is depicted below.
  \if{html}{\figure{scholz.svg}}
  \if{latex}{\figure{scholz.png}{options: width=5in}} Here
  S = stem cells, CG = colony forming units of granulocytes
  and macrophages, PGB = proliferating granulopoietic
  blasts, MGB = maturing granulopoietic blasts subdivided
  into metamyelocytes (G4),banded granulocytes (G5) and
  segmented granulocytes (G6) and GRA = circulating
  granulocytes. The system is regulated by feedback loops.
  A major loop is mediated by G-CSF which is produced
  endogenously but can also be applied subcutaneously.
}
\note{
  This work was supported by the National Cancer Institute
  and Tufts Integrative Cancer Biology Program under
  U54CA149233-029689.
}
\examples{
\dontrun{
library(myelo)
library(deSolve)
pars=scholzPars

X0=rep(0,14)
names(X0)<-c(paste("A",1:7,sep=""),paste("M",1:7,sep=""))
times <- seq(1, 3000, by = 1)
out   <- ode(X0, times, scholz12, pars)
tail(out)

}
}
\author{
  Tom Radivoyevitch (\email{txr24@case.edu})
}
\references{
  M. Scholz, S. Schirm, M. Wetzler, C. Engel and M.
  Loeffler, Pharmacokinetic and -dynamic modelling of G-CSF
  derivatives in humans, \emph{Theoretical Biology and
  Medical Modelling} \bold{9} 32 (2012).
}
\seealso{
  \code{\link{myelo-package}, \link{rcj12}}
}
\keyword{IO}

