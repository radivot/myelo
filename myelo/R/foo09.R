#'CML model of Foo et al. 2009
#'
#'This function returns the right hand side of the CML ODE
#'model of Foo et al 2009. The intended use of this
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
#'element is the vector of derivatives and the second element of the list is a vector of auxiliary
#'variables that one may wish to track over time.
#'@note This work was supported by the National Cancer Institute and Tufts
#'Integrative Cancer Biology Program under U54CA149233-029689.
#'@author Tom Radivoyevitch (\email{txr24@@case.edu})
#'@seealso \code{\link{myelo-package}, \link{rcj12}}
#'@references Jasmine Foo, Mark W. Drummond, Bayard Clarkson, Tessa Holyoake, Franziska Michor,
#'Eradication of Chronic Myeloid Leukemia Stem Cells: 
#'A Novel Mathematical Model Predicts No Therapeutic
#'Benefit of Adding G-CSF to Imatinib, \emph{PLOS Computational Biology},
#'\bold{5(9)}, (2009)
#'@keywords IO
#'@export
#'@examples
#'
#'\dontrun{
#'
#'
#'library(myelo)  
#'pars=c(d0=0.003,d1=0.008,d2=0.05,d3=1,ax=0.8,bx=5,cx=100,rx=0.005,ry=0.008,
#'ay=1.6, byy=10, cy=100, trt=0,X0o=1e6,Y0o=1e7,alf=0.0001,bet=0.0009)
#'(y0<-c(X0=9e5,Xq=1e5,X1=1e8,X2=1e10,X3=1e12,Y0=9e6,Yq=1e6,Y1=0,Y2=0,Y3=0))
#'out=ode(y=y0,times=seq(0,4000,10),foo09, parms=pars, rtol=1e-4, atol= 1e-4)
#'plot(out)
#' }
#'
foo09<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    if (trt==1) {   #During therapy 
      ay=ay/100
      byy=byy/750
      ry=ry/4
    }
    px=(rx/d0-1)/X0o
    py=(ry/d0-1)/Y0o
    phix=1/(1+px*(X0+Y0))
    phiy=1/(1+py*(X0+Y0))
    dX0 = (rx*phix-d0-alf)*X0+bet*Xq
    dXq = bet*X0
    dX1 = ax*X0-d1*X1
    dX2 = bx*X1-d2*X2
    dX3 = cx*X2-d3*X3
    dY0 = (ry*phiy-d0-alf)*Y0+bet*Yq
    dYq = alf*Y0-bet*Yq
    dY1 = ay*Y0-d1*Y1
    dY2 = byy*Y1-d2*Y2
    dY3 = cy*Y2-d3*Y3
    return(list(c(dX0, dXq, dX1, dX2, dX3, dY0, dYq, dY1, dY2, dY3),c(ratio=(Y3)/(Y3+X3))))
  })
}

