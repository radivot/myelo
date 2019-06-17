#'Neutrophil-GCSF-Chemo model of Craig et al 2015
#'
#'This function returns the right hand side of an ordinary differential delay
#'equation model of Craig et al. JTB 2015. 
#'The intended use of this function is as an argument to \code{dede()} of the \code{deSolve} package.
#'
#' In this model Q = stem cells, Nr = marrow reserve neutrophils, and N = circulating neutrophils
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
#'@note This work was supported by the Cleveland Clinic Foundation.
#'@author Tom Radivoyevitch (\email{txr24@@case.edu})
#'@seealso \code{\link{myelo-package}, \link{craig15}}


#'@references  # Morgan Craig ... Michael C. Mackey, Neutrophil dynamics during
#'  concurrent chemotherapy and G-CSF administration: Mathematical modelling
#'  guides dose optimisation to minimise neutropenia, \emph{Journal of
#'  Theoretical Biology} \bold{385} 77-89  (2015).
#'
#'@keywords IO
#'
#'@export
#'@examples
#'
#'\dontrun{
#'library(deSolve)
#'library(myelo)
#' craigPars15
#' times <- seq(-30,200,by=0.1)
#' yout <- dede(c(Q=craigPars15[["Qss"]],Nr=craigPars15[["Nrss"]],N=craigPars15[["Nss"]]),
#' 		times = times, func = craig15,	parms = craigPars15)
#' plot(yout)
#'}
#'
craig15<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
	      As=2*exp(-gamS*tauS)
				if (Time < 0) {
				  etaN=etaNbar*(N/the1script)^nu1/(1+(N/the1script)^nu1)
				  etaP=etaPbar*(P/the4script)^nu4/(1+(P/the4script)^nu4)
				  An=Anmax*exp(-etaN*tauNM)
				  Ap=Apmax*exp(-etaP*tauPM)
				  beta=k0/(1+(S/the2)^s2)
				  kN=f0/(1+(N/the1)^s1)
				  kP=Kpbar/(1+Kp*P^s4)
				  dS=-(beta+kN+kP+kR)*S + As*beta*S
				  dN=-gamN*N + (An/1e2)*kN*S  # 1e2 map units of S to N
				  dP=-gamP*P + (Ap/1e4)*(kP*S-exp(-gamP*tauPS)*kP*S)# 1e4 maps units of S to P
				}	else {
				  St=lagvalue(Time - tauS,1)
				  Snt=lagvalue(Time - tauNM,1)
				  Spt=lagvalue(Time - tauPM,1)
				  Sptsum=lagvalue(Time - (tauPM+tauPS),1)
				  Nt=lagvalue(Time - tauNM,2)
				  Pt=lagvalue(Time - tauPM,3)
				  Ptsum=lagvalue(Time - (tauPM+tauPS),3)
				  etaNt=etaNbar*(Nt/the1script)^nu1/(1+(Nt/the1script)^nu1)
				  etaPt=etaPbar*(Pt/the4script)^nu4/(1+(Pt/the4script)^nu4)
				  etaPtsum=etaPbar*(Ptsum/the4script)^nu4/(1+(Ptsum/the4script)^nu4)
				  Ant=Anmax*exp(-etaNt*tauNM)
				  Apt=Apmax*exp(-etaPt*tauPM)
				  Aptsum=Apmax*exp(-etaPtsum*tauPM)
				  kN=f0/(1+(N/the1)^s1)
				  kNt=f0/(1+(Nt/the1)^s1)
				  kP=Kpbar/(1+Kp*P^s4)
				  kPt=Kpbar/(1+Kp*Pt^s4)
				  kPtsum=Kpbar/(1+Kp*Ptsum^s4)
				  beta=k0/(1+(S/the2)^s2)
				  betat=k0/(1+(St/the2)^s2)
				  dS=-(beta+kN+kP+kR)*S + As*betat*St
				  dN=-gamN*N + (Ant/1e2)*kNt*Snt
				  dP=-gamP*P + (Apt/1e4)*(kPt*Spt-exp(-gamP*tauPS)*kPtsum*Sptsum) # in paper
				  # dP=-gamP*P + (Apt/1e4)*kPt*Spt-(Aptsum/1e4)*exp(-gamP*tauPS)*kPtsum*Sptsum) # should be?
					}
				list(c(dS,dN,dP))
			})
}



