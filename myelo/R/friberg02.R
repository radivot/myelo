#'Friberg's myelosuppression model (JCO 2002)
#'
#'
#'This function returns the right hand side of an ordinary differential equation
#'model of docetaxel killing proliferating neutrophils.  The intended use of this function is as
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
#'@references Lena E. Friberg, Anja Henningsson, Hugo Maas, Laurent Nguyen, and Mats O. Karlsson
#'Model of Chemotherapy-Induced Myelosuppression With Parameter Consistency Across Drugs, 
#' \emph{Journal of Clinical Oncology},  \bold{20},  4713-4721 (2002)
#'@keywords IO
#'@name friberg02
#'@export
#'@import deSolve

friberg02<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dC1=-(k12+k13+k10)*C1 + k21*C2 + k31*C3
    dC2=k12*C1 - k21*C2
    dC3=k13*C1 - k31*C3
    Cp=C1/V1/mw  #mw in mg/umole => uM
    Edrug=slope*Cp
    dProl = ktr*Prol*(1-Edrug)*(Circ0/Circ)^gam - ktr*Prol
    dTrans1=ktr*Prol-ktr*Trans1
    dTrans2=ktr*Trans1-ktr*Trans2
    dTrans3=ktr*Trans2-ktr*Trans3
    dCirc=ktr*Trans3-ktr*Circ
    return(list(c(dC1,dC2,dC3,dProl,dTrans1,dTrans2,dTrans3,dCirc)))
  })
}

