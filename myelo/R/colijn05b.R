#'Cyclic Neutropenia model of Colijn and Mackey (JTB 2005)
#'
#'This function returns the right hand side of an ordinary differential delay
#'equation model. 
#'The intended use of this function is as an argument to \code{dede()} of the \code{deSolve} package.
#'
#' In this model Q = quiescent stem cells, N = neutrophils, R = red blood cells, and P = platelets
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


#'@references  Caroline Colijn and Michael C. Mackey,
#'A mathematical model of hematopoiesis: II. Cyclical neutropenia, 
#'\emph{Journal of Theoretical Biology} 
#'\bold{237} 133-146 (2005).
#'
#'@keywords IO
#'@export
#'@examples
#'
#'\dontrun{
#'library(deSolve)
#'library(myelo)
#' colijnPars05b
#' times <- seq(-200,200,by=0.1)
#' yout <- dede(c(Q=colijnPars05b[["Qss"]],N=colijnPars05b[["Nss"]],
#'    R=colijnPars05b[["Rss"]],P=colijnPars05b[["Pss"]]),
#' 		times = times, func = colijn05b,	parms = colijnPars05b)
#' plot(yout)
#'}
#'
colijn05b<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				if (Time < 0) {
				  dQ=0
				  dN=0
				  dR=0
				  # dL=0   #since Retics L do not feedback on anything, we can skip them in this model
				  dP=0
				}	else {
				  tauRS=tauSum-tauRM-tauRet
				  
				  Qt=lagvalue(Time - tauS,1)
				  Qnt=lagvalue(Time - tauNM,1)
				  Qrt=lagvalue(Time - tauRM,1)
				  Qpt=lagvalue(Time - tauPM,1)
				  Qptsum=lagvalue(Time - (tauPM+tauPS),1)
				  Qrtsum=lagvalue(Time - (tauRM+tauRS),1)
				  # Qrtsum2=lagvalue(Time - (tauSum),1)
				  
				  Nt=lagvalue(Time - tauNM,2)
				  
				  Rt=lagvalue(Time - tauRM,3)
				  Rtsum=lagvalue(Time - (tauRM+tauRS),3)
				  # Rtsum2=lagvalue(Time - tauSum,3)
				  
				  Pt=lagvalue(Time - tauPM,4)
				  Ptsum=lagvalue(Time - (tauPM+tauPS),4)
				  
				  # etaN=etaNbar/(1+(Nt/the1script)^nu1)
				  # etaP=etaPbar/(1+(Pt/the4script)^nu4)
   	      Aq=2*exp(-gamS*tauS)
				  # Ant=Anmax*exp(-etaN*tauNM)
				  # Apt=Apmax*exp(-etaP*tauPM)
				  
				  kN=f0/(1+(N/the1)^n)
				  kNt=f0/(1+(Nt/the1)^n)
				  
				  kR=Krbar/(1+Kr*R^m)
				  kRt=Krbar/(1+Kr*Pt^m)
				  kRtsum=Krbar/(1+Kr*Rtsum^m)
				  # kRtsum2=Krbar/(1+Kr*Rtsum2^m) #bigger tau sum
				  
          kP=Kpbar/(1+Kp*P^r)
          kPt=Kpbar/(1+Kp*Pt^r)
          kPtsum=Kpbar/(1+Kp*Ptsum^r)
          beta=k0/(1+(Q/the2)^s)
          betat=k0/(1+(Qt/the2)^s)
          dQ=-(beta+kN+kP+kR)*Q + Aq*betat*Qt
				  dN=-gamN*N + An*kNt*Qnt
				  dR=-gamR*R + Ar*(kRt*Qrt-exp(-gamR*tauRS)*kRtsum*Qrtsum) 
				  dP=-gamP*P + Ap*(kPt*Qpt-exp(-gamP*tauPS)*kPtsum*Qptsum) 
					}
				list(c(dQ,dN,dR,dP))
			})
}



