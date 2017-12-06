#'LSC-HSC Competition Model of Stiehl et al. 2015
#'
#'
#'This function returns the right hand side of an ordinary differential equation
#'model of LSC competing with HSC.  The intended use of this function is as
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
#'@references Thomas Stiehl, Natalia Baran, Anthony D. Ho, and Anna Marciniak-Czochra,
#'  Cell Division Patterns in Acute Myeloid Leukemia Stem-like Cells Determine
#'  Clinical Course: A Model to Predict Patient Survival,  \emph{Cancer Research},
#'  \bold{75}, 6 (2015)
#'@keywords IO
#'@name stiehl15
#'@export
#'@import deSolve

stiehl15<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)),{
    s=1/(1+kc*C3+kl*L3)
    dC1=(2*a1c*s-1)*p1c*C1 #parent is always removed (-1) and both daughters either renew or go to C2 (first term below)
    dC2=2*(1-a1c*s)*p1c*C1 + (2*a2c*s-1)*p2c*C2 #2nd term is like one above
    dC3=2*(1-a2c*s)*p2c*C2 - d3c*C3 #leave marrow if become C3, only C3 can die
    dL1=(2*a1l*s-1)*p1l*L1
    dL2=2*(1-a1l*s)*p1l*L1 + (2*a2l*s-1)*p2l*L2
    dL3=2*(1-a2l*s)*p2l*L2 - d3l*L3
    list(c(dC1,dC2,dC3,dL1,dL2,dL3),
         c(s=s,B=L3/(C1+C2+C3+L1+L2+L3))
    )
  }) 
}
