#'Neutrophil and Platelet model of Zhuge et al 2019
#'
#'This function returns the right hand side of an ordinary differential delay
#'equation model of Zhuge et al. JTB 2019. 
#'The intended use of this function is as an argument to \code{dede()} of the \code{deSolve} package.
#'
#' In this model S = stem cells, N = neutrophils, and P = platlets
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
#'@seealso \code{\link{myelo-package}, \link{zhuge12}}

#'@references  Changjing Zhuge, Michael C. Mackey, and Jinzhi Lei,
#'Origins of oscillation patterns in cyclical thrombocytopenia, 
#'\emph{Journal of Theoretical Biology} 
#' \bold{462} 432-445 (2019).
#'
#'@keywords IO
#'@export

zhuge19<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
	      Aq=2*exp(-gamS*tauS)
				if (Time < 0) {
				  etaN=etaNbar/(1+(Nss/the1script)^nu1)
				  etaP=etaPbar/(1+(Pss/the4script)^nu4)
				  An=Anmax*exp(-etaN*tauNM)
				  Ap=Apmax*exp(-etaP*tauPM)
					beta=k0/(1+(Sss/the2)^s2)
					Kn=f0/(1+(Nss/the1)^s1)
					Kp=Kpbar/(1+Kp*Pss^s4)
					dS=-(beta+Kn+Kp+kR)*S + Aq*beta*Sss
					dN=-gamN*N + An*Kn*Sss
					dP=-gamP*P + Ap*(1-exp(-gamP*tauPS))*Kp*Sss 
				}	else {
				  St=lagvalue(Time - tauS,1)
				  Snt=lagvalue(Time - tauNM,1)
				  Spt=lagvalue(Time - tauPM,1)
				  Sptsum=lagvalue(Time - (tauPM+tauPS),1)
				  Nt=lagvalue(Time - tauNM,2)
				  Pt=lagvalue(Time - tauPM,3)
				  Ptsum=lagvalue(Time - (tauPM+tauPS),3)
				  etaN=etaNbar/(1+(Nt/the1script)^nu1)
				  etaP=etaPbar/(1+(Pt/the4script)^nu4)
				  Ant=Anmax*exp(-etaN*tauNM)
				  Apt=Apmax*exp(-etaP*tauPM)
				  Kn=f0/(1+(N/the1)^s1)
				  Knt=f0/(1+(Nt/the1)^s1)
          Kp=Kpbar/(1+Kp*P^s4)
          Kpt=Kpbar/(1+Kp*Pt^s4)
          Kptsum=Kpbar/(1+Kp*Ptsum^s4)
          beta=k0/(1+(S/the2)^s2)
          betat=k0/(1+(St/the2)^s2)
          dS=-(beta+Kn+Kp+kR)*S + Aq*betat*St
				  dN=-gamN*N + Ant*Knt*Snt
				  dP=-gamP*P + Apt*(Kpt*Spt-exp(-gamP*tauPS)*Kptsum*Sptsum) 
					}
				list(c(dS,dN,dP))
			})
}



