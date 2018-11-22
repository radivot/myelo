#'CML-T-Cell Model of Besse et al 2018
#'
#'
#'This function returns the right hand side of an ordinary differential equation
#'model of T cells killing CML cells.  The intended use of this function is as
#'an argument to \code{ode()} of the \code{deSolve} package.
#'
#'
#'@param Time The parameters of this model do not depend on time: this argument
#'  exists here as a dummy only because deSolve expects it to exist.
#'@param State Vector of current states. The elements of this vector must be
#'  named because within this function's definition is
#'  \code{with(as.list(State,Pars), {XXX} )} and code in XXX presupposes that it
#'  can call named elements as variables with those names.
#'@param Pars Vector of parameter values with \emph{named} (see above) elements.
#'@return A list of length 2 (as expected by deSolve) where the first list
#'  element is the vector of derivatives, i.e. the right hand side of the ODEs,
#'  and the second element of the list is a vector of auxiliary variables that
#'  one may wish to track over time.
#'@author Tom Radivoyevitch
#'@seealso \code{\link{myelo-package}, \link{clapp15}, \link{michor05}}
#'@references Apollos Besse, Geoffrey D. Clapp, Samuel Bernard, Franck E.
#'  Nicolini, Doron Levy, Thomas Lepoutre, Stability Analysis of a Model of
#'  Interaction Between the Immune System and Cancer Cells in Chronic
#'  Myelogenous Leukemia, \emph{Bull Math Biol}, \bold{80}, 1084-1110 (2018)
#'@keywords IO
#'@name besse18
#'@export
#'@import deSolve

besse18<-function(Time, State, Pars) {
  with(as.list(c(Time, State, Pars)),{
    if (Time>0) a1=a1/inh1
    dY1 = r*Y1*(1-Y1/K)                         -  mu*Y1*Z
    dY2 = a1*Y1         - d2*Y2                 -  mu*Y2*Z
    dZ  = s             -  d*Z                  + alf*Y2*Z/(1+eps*Y2^2)
    list(c(dY1,dY2,dZ),c(ratio=log10(100*beta*Y2/(Y2+2*x))))
  }) 
}
