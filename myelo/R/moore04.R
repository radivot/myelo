#'CML-T-Cell Model of Moore and Li 2004
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
#'@seealso \code{\link{myelo-package}, \link{fokas91}}
#'@references Helen Moorea and Natasha K. Li  
#'A mathematical model for chronic myelogenous leukemia (CML) and T cell interaction, 
#' \emph{Journal of Theoretical Biology},  \bold{227},  513–523 (2004)
#'@keywords IO
#'@name moore04
#'@export
#'@import deSolve

moore04<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)),{
    dTn = sn - dn*Tn - kn*Tn*C/(C+eta) 
    dTe = alfn*kn*Tn*C/(C+eta) + alfe*Te*C/(C+eta) - de*Te - ge*C*Te
    dC  = rc*C*log(Cmax/C) - dc*C - gc*C*Te 
    list( c(dTn,dTe,dC), c(T=Tn+Te) )
  }) 
}
