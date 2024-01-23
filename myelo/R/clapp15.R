#'CML-T-Cell Model of Clapp et al 2015
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
#'@seealso \code{\link{myelo-package}, \link{michor05}}
#'@references Geoffrey D. Clapp, Thomas Lepoutre, Raouf El Cheikh, Samuel Bernard, Jeremy Ruby, 
#'  Helene Labussie Wallet, Franck E. Nicolini, and Doron Levy, Implication of the Autologous
#'  Immune System in BCR-ABL Transcript Variations in Chronic Myelogenous
#'  Leukemia Patients Treated with Imatinib, \emph{Cancer Research},  \bold{75}, 4053-4062 (2015)
#'@keywords IO
#'@name clapp15
#'@export
#'@import deSolve

clapp15<-function(Time, State, Pars) {
  with(as.list(c(Time, State, Pars)),{
    if (Time>0) {a1=a1/inh1; a2=a2/inh2}
    dY0 = b1*Y1 - a0*Y0                         -  mu*Y0*Z/(1+eps*Y3^2)
    dY1 = a0*Y0 - b1*Y1 - d1*Y1 + r*Y1*(1-Y1/K) -  mu*Y1*Z/(1+eps*Y3^2)
    dY2 = a1*Y1         - d2*Y2                 -  mu*Y2*Z/(1+eps*Y3^2)
    dY3 = a2*Y2         - d3*Y3                 -  mu*Y3*Z/(1+eps*Y3^2)
    dZ  = sz            - dz*Z                  + alf*Y3*Z/(1+eps*Y3^2)
    list(c(dY0,dY1,dY2,dY3,dZ),c(ratio=log10(100*beta*Y3/(Y3+2*x))))
  })  
}
