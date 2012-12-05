#'PKPD G-CSF model of Scholz et al. TBMM 2012
#'
#'This function returns the right hand side of an ordinary differential
#'equation model of G-CSF published by M. Scholz, S. Schirm, M. Wetzler, C. Engel
#'and M. Loeffler in Theoretical Biology and Medical Modelling 2012, 9:32.  The intended use of this
#'function is as an argument to \code{ode()} of the \code{deSolve} package.
#'
#'This model is depicted below.  
#'\if{html}{\figure{scholz.svg}} \if{latex}{\figure{scholz.png}{options: width=5in}}
#' Here S = stem cells, CG = colony forming units of granulocytes
#' and macrophages, PGB = proliferating granulopoietic blasts, MGB = maturing granulopoietic blasts
#' subdivided into metamyelocytes (G4),banded granulocytes (G5) and segmented granulocytes (G6) and GRA =
#' circulating granulocytes. The system is regulated by feedback loops. A major loop is mediated by G-CSF
#' which is produced endogenously but can also be applied subcutaneously.
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
#'element is the vector of derivatives, i.e. the right hand side of the ODEs 
#' and the second element of the list is a vector of auxiliary
#'variables that one may wish to track over time.
#'@note This work was supported by the National Cancer Institute and Tufts
#'Integrative Cancer Biology Program under U54CA149233-029689.
#'@author Tom Radivoyevitch (\email{txr24@@case.edu})
#'@seealso \code{\link{myelo-package}, \link{rcj12}}
#'@references  M. Scholz, S. Schirm, M. Wetzler, C. Engel
#'and M. Loeffler, Pharmacokinetic and -dynamic modelling of
#'G-CSF derivatives in humans, \emph{Theoretical Biology and Medical Modelling} 
#' \bold{9} 32 (2012).
#'@keywords IO
#'@export
#'@examples
#'
#'\dontrun{
#'library(myelo)
#'library(deSolve)
#'pars=scholzPars 
#' 
#'X0=rep(0,14)
#'names(X0)<-c(paste("A",1:7,sep=""),paste("M",1:7,sep=""))
#'times <- seq(1, 3000, by = 1)
#'out   <- ode(X0, times, scholz12, pars)
#'tail(out)  
#'
#'}
#'
scholz12<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				dA1 = fb*q -  A1/T
				dA2 = 2*fb*A1/T - A2/T
				dA3 = 2*fb*A2/T - A3/T
				dA4 = 2*fp*A3/T - A4/T
				dA5 = 2*fp*A4/T - A5/T
				dA6 = 2*fm*A5/T - A6/T
				dA7 = 2*fm*A6/T - A7/T
				dM1 = (1-fb)*q -  M1/T 
				dM2 = 2*(1-fb)*A1/T - M2/T  + M1/T
				dM3 = 2*(1-fb)*A2/T - M3/T  + M2/T
				dM4 = 2*(1-fp)*A3/T - M4/T  + M3/T
				dM5 = 2*(1-fp)*A4/T - M5/T  + M4/T
				dM6 = 2*(1-fm)*A5/T - M6/T  + M5/T
				dM7 = 2*(1-fm)*A6/T - M7/T  + M6/T
				Nb=A1+A2+A3+M1+M2+M3
				Np=A4+A5+M4+M5
				Nm=A6+A7+M6+M7
				Ntot=Nb+Np+Nm
				Q=(2*A7+M7)/T
				return(list(c(dA1,dA2,dA3,dA4,dA5,dA6,dA7,dM1,dM2,dM3,dM4,dM5,dM6,dM7),
							c(Nb=Nb,Np=Np,Nm=Nm,povb=Np/Nb,movb=Nm/Np,Qovm=Q/Nm,Qovq=Q/q,Ntot=Ntot,Q=Q)))
			})
}

