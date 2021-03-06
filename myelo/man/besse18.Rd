% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/besse18.R
\name{besse18}
\alias{besse18}
\title{CML-T-Cell Model of Besse et al 2018}
\usage{
besse18(Time, State, Pars)
}
\arguments{
\item{Time}{The parameters of this model do not depend on time: this argument
exists here as a dummy only because deSolve expects it to exist.}

\item{State}{Vector of current states. The elements of this vector must be
named because within this function's definition is
\code{with(as.list(State,Pars), {XXX} )} and code in XXX presupposes that it
can call named elements as variables with those names.}

\item{Pars}{Vector of parameter values with \emph{named} (see above) elements.}
}
\value{
A list of length 2 (as expected by deSolve) where the first list
 element is the vector of derivatives, i.e. the right hand side of the ODEs,
 and the second element of the list is a vector of auxiliary variables that
 one may wish to track over time.
}
\description{
This function returns the right hand side of an ordinary differential equation
model of T cells killing CML cells.  The intended use of this function is as
an argument to \code{ode()} of the \code{deSolve} package.
}
\references{
Apollos Besse, Geoffrey D. Clapp, Samuel Bernard, Franck E.
 Nicolini, Doron Levy, Thomas Lepoutre, Stability Analysis of a Model of
 Interaction Between the Immune System and Cancer Cells in Chronic
 Myelogenous Leukemia, \emph{Bull Math Biol}, \bold{80}, 1084-1110 (2018)
}
\seealso{
\code{\link{myelo-package}, \link{clapp15}, \link{michor05}}
}
\author{
Tom Radivoyevitch
}
\keyword{IO}
