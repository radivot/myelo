% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colijn05b.R
\name{colijn05b}
\alias{colijn05b}
\title{Cyclic Neutropenia model of Colijn and Mackey (JTB 2005)}
\usage{
colijn05b(Time, State, Pars)
}
\arguments{
\item{Time}{The current time.}

\item{State}{Vector of current states. The elements of this vector must be
named because within this function's definition is
\code{with(as.list(State,Pars), {XXX} )} and code in XXX presupposes that it
can call named elements as variables with those names.}

\item{Pars}{Vector of parameter values with \emph{named} (see above)
elements.}
}
\value{
A list of length 1 where the first list
element is the vector of derivatives, i.e. the right hand side of the ODEs
}
\description{
This function returns the right hand side of an ordinary differential delay
equation model. 
The intended use of this function is as an argument to \code{dede()} of the \code{deSolve} package.
}
\details{
In this model Q = quiescent stem cells, N = neutrophils, R = red blood cells, and P = platelets
}
\note{
This work was supported by the Cleveland Clinic Foundation.
}
\examples{

\dontrun{
library(deSolve)
library(myelo)
colijnPars05b
times <- seq(-200,200,by=0.1)
yout <- dede(c(Q=colijnPars05b[["Qss"]],N=colijnPars05b[["Nss"]],
   R=colijnPars05b[["Rss"]],P=colijnPars05b[["Pss"]]),
		times = times, func = colijn05b,	parms = colijnPars05b)
plot(yout)
}

}
\references{
Caroline Colijn and Michael C. Mackey,
A mathematical model of hematopoiesis: II. Cyclical neutropenia, 
\emph{Journal of Theoretical Biology} 
\bold{237} 133-146 (2005).
}
\seealso{
\code{\link{myelo-package}, \link{zhuge12}}
}
\author{
Tom Radivoyevitch (\email{txr24@case.edu})
}
\keyword{IO}
