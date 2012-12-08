\name{zhuge12}
\alias{zhuge12}
\title{Neutrophil dynamics model of Zhuge et al 2012}
\usage{
  zhuge12(Time, State, Pars)
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
  This function returns the right hand side of an ordinary
  differential delay equation model of neutrophil dynamics
  published by C. Zhuge et al. in 2012 The intended use of
  this function is as an argument to \code{dede()} of the
  \code{deSolve} package.
}
\details{
  In this model Q = quiescent stem cells and N =
  neutrophils
}
\note{
  This work was supported by the National Cancer Institute
  and Tufts Integrative Cancer Biology Program under
  U54CA149233-029689.
}
\examples{
\dontrun{
library(myelo)
times <- c(-zhugePars[["tauN"]]:600)
yout <- dede(c(Q=zhugePars[["Qss"]],N=zhugePars[["Nss"]]), times = times, func = zhuge12,	parms = zhugePars)
plot(yout)

}
}
\author{
  Tom Radivoyevitch (\email{txr24@case.edu})
}
\references{
  Changjing Zhuge, Jinzhi Lei, and Michael C. Mackey,
  Neutrophil dynamics in response to chemotherapy and
  G-CSF, \emph{Journal of Theoretical Biology} \bold{293}
  111–120 (2012).
}
\seealso{
  \code{\link{myelo-package}, \link{rcj12}}
}
\keyword{IO}

