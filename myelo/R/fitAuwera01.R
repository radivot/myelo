#'Inspired by Brooks et al 2012, a model that fits the G-CSF data of Auwera 2001
#'
#'This function returns the right hand side of a barebone ordinary differential 
#'equation model of neutrophil dynamics after G-CSF.  
#'The intended use of this function is as an argument to \code{ode()} of the \code{deSolve} package.
#'
#'In this model CD34 = CD34+ cells, Nm = mature neutrophils in marrow, and Nb = neutrophils in blood
#'
#'@param Time The current time.
#'@param State Vector of current states. The elements of this vector must be
#'named because within this function's definition is
#'\code{with(as.list(State,Pars), {XXX} )} and code in XXX presupposes that it
#'can call named elements as variables with those names.
#'@param Pars Vector of parameter values with \emph{named} (see above)
#'elements.
#'@return A list of length 1 where the first list
#'element is the vector of derivatives, i.e. the right hand side of the ODEs 
#'@note This work was supported by the National Cancer Institute and Tufts
#'Integrative Cancer Biology Program under U54CA149233-029689.
#'@author Tom Radivoyevitch (\email{txr24@@case.edu})
#'@seealso \code{\link{myelo-package}, \link{brooks12}, \link{zhuge12}}
#'@references  Phillippe van der Auwera, Erich Platzer, Z.-X. Xu,1 Rainer Schulz, Olivier Feugeas,
#' Renaud Capdeville, and David J. Edwards,
#' Pharmacodynamics and Pharmacokinetics of Single
#' Doses of Subcutaneous Pegylated Human G-CSF Mutant
#' (Ro 25-8315) in Healthy Volunteers: Comparison With
#' Single and Multiple Daily Doses of Filgrastim
#' \emph{American Journal of Hematology} 
#' \bold{66} 245–251 (2001).
#'@keywords IO
#'@export
#'@examples
#'
#'\dontrun{
#'library(myelo)
#' 
#' times <- c(-20:-1,seq(-.1,1,by=0.02),2:14,seq(14.01,15,by=0.02),16:200)
#' yout <- ode(c(CD34=1e6,Nm1=3e12,Nm2=1e10,Nm3=9e11,Nm=1e11,CD34b=3,Nb=5000,Gt=0,Gb=0),
#' 		times = times, func = fitAuwera01,	parms = c(Gain=3e6,tau=12)
#' par(mfrow=c(1,1),mar=c(4.7,6,3.3,1),oma=c(0,0,0,0),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
#' plot(times,yout[,3]/1e8,type="l",xlab="days",ylab="Neutrophils")
#' 
#'
#'}
#'

fitAuwera01<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				Gain=3e12+eta*Gb
				tau= 12 - tauSlope*Gb
				I0g=0
				if (Time < 0) {
					delEtaNP=etaNP*tauNP
					delGam0=gam0*tauNMmax
					delGamS=gamS*tauS
					An=exp(delEtaNP - delGam0)
					Aq=2*exp(-delGamS)
					dCD34=-(k0/(1+(Q/the2)^s2) + f0/(1+(N/the1)^s1)+kdel)*Q + Aq*k0/(1+(Qss/the2)^s2)*Qss
					dN=-gamN*N + An*f0/(1+(Nss/the1)^s1)*Qss
					dC=0;	dGt=0;	dGb=0
				}	else {
					if ((Time%%T > T1)&&(Time%%T < T1+Delg) ) I0g=Dg/Delg #in G-CSF injection period
					dCD34=-kmob*Gb*CD34
					dNb=-gamN*N + An*f0/(1+(Ntn/the1)^s1)*Qtn
					dGt= I0g -kt*Gt
					dG=  kt*Gt/VB -kb*Gb - sig*Nb*G^2/(kG+G^2)
				}
				dEtaNP=etaNPv
				dGam0=gam0v
				dGamS=gamSv
				list(c(dQ,dN,dEtaNP,dGam0,dGamS,dG,dX,dC),
						c(delEtaNP=delEtaNP,delGam0=delGam0,delGamS=delGamS,
								etaNPv=etaNPv,gam0v=gam0v,gamSv=gamSv,trt=trt,tauNMv=tauNMv,tauNv=tauNv,An=An,Aq=Aq))
			})
}



