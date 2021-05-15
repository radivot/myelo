#'CML model of Michor et al. 2005
#'
#'This function returns the right hand side of an ordinary differential
#'equation model of CML that was published by Michor et al. in 2005.
#'The intended use of this
#'function is as an argument to \code{ode()} of the \code{deSolve} package.
#'
#'
#'@param Time The parameters of this model do not depend on time: this argument
#'exists here as a dummy only because deSolve expects it to exist.
#'@param State Vector of current states. The elements of this vector must be
#'named because within this function's definition is
#'\code{with(as.list(State,Pars), {XXX} )} and code in XXX presupposes that it
#'can call named elements as variables with those names.
#'@param Pars Vector of parameter values with \emph{named} (see above)
#'elements.
#'@return A list of length 2 (as expected by deSolve) where the first list
#'element is the vector of derivatives, i.e. the right hand side of the ODEs, 
#'and the second element of the list is a vector of auxiliary
#'variables that one may wish to track over time.
#'@note This work was supported by the National Cancer Institute and Tufts
#'Integrative Cancer Biology Program under U54CA149233-029689.
#'@author Tom Radivoyevitch (\email{txr24@@case.edu})
#'@seealso \code{\link{myelo-package}, \link{rcj12}}
#'@references Franziska Michor, Timothy P. Hughes, Yoh Iwasa, Susan Branford,
#' Neil P. Shah, Charles L. Sawyers & Martin A. Nowak
#'Dynamics of chronic myeloid leukaemia, \emph{Nature},
#'\bold{435}, 1267-1270 (2005)
#'@keywords IO
#'@export
#'@examples
#'
#'\dontrun{
#'library(myelo)  
#'
#'pars=c(d0=0.003, d1=0.008, d2=0.05, d3=1,ax=0.8, bx=5, cx=100, ry=0.008, 
#'D0a=45/99, D0b= 45/749,  
#'rz=0.023, az=1.6, bz=10,
#'u=4e-8, sp=2e6)
#'(y0<-c(X0=2e6,X1=0,X2=0,X3=0,
#' Y0=0,Y1=0,Y2=0,Y3=0,
#' Z0=0,Z1=0,Z2=0,Z3=0,D=0))
#'out=ode(y=y0,times=seq(0,1000,1),michor05, parms=pars, rtol=1e-4, atol= rep(1e-4,12))
#'plot(out)
#' }
#'
michor05<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    # lambda is a controller that manipulates inflow
    # to move a level (x) toward setpoint (sp)
    lambda= -0.5*(X0-sp)
    dX0 = (lambda-d0)*X0 
    dX1 = ax*X0-d1*X1
    dX2 = bx*X1-d2*X2
    dX3 = cx*X2-d3*X3
    dY0 = (ry*(1-u)-d0)*Y0
    dY1 = Y0*az/(1+D/D0a) - d1*Y1
    dY2 = Y1*bz/(1+D/D0b) - d2*Y2
    dY3 = cx*Y2 - d3*Y3
    dZ0 = (rz-d0)*Z0+ry*u*Y0
    dZ1 = az*Z0-d1*Z1
    dZ2 = bz*Z1-d2*Z2
    dZ3 = cx*Z2-d3*Z3
    dD  = 0
    return(list(c(dX0,dX1, dX2, dX3, 
                  dY0,dY1, dY2, dY3, 
                  dZ0,dZ1, dZ2, dZ3,dD),
                c(ratio=(Y3+Z3)/(Y3+Z3+X3))))
  })
}
