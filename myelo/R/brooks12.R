#'Neutrophil dynamics model of Brooks et al 2012
#'
#'This function returns the right hand side of an ordinary differential delay
#'equation model of neutrophil dynamics published by Brooks et al. in 2012 
#'The intended use of this function is as an argument to \code{dede()} of the \code{deSolve} package.
#'
#' In this model Q = quiescent stem cells and N = neutrophils
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
#'@seealso \code{\link{myelo-package}, \link{rcj12}}
#'@references  Grace Brooks, Gabriel Provencher, Jinzhi Lei, Michael C. Mackey,
#'Neutrophil dynamics after chemotherapy and G-CSF: The role of pharmacokinetics 
#' in shaping the response, \emph{Journal of Theoretical Biology} 
#' \bold{315} 97-109 (2012).
#'@keywords IO
#'@export
#'@examples
#'
#'\dontrun{
#'library(myelo)
#' 
#' times <- c(-20:-1,seq(-.1,1,by=0.02),2:14,seq(14.01,15,by=0.02),16:100)
#' yout <- dede(c(Q=brooksPars[["Qss"]],N=brooksPars[["Nss"]],EtaNP=0,Gam0=0,GamS=0,G=0,X=0,C=0),
#' 		times = times, func = brooks12,	parms = brooksPars)
#' # this should make Figure 3A. The response appears to be off by a dose-amplitude parameter
#' par(mfrow=c(1,1),mar=c(4.7,6,3.3,1),oma=c(0,0,0,0),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
#' plot(times,yout[,3]/1e8,type="l",xlab="days",ylab="Neutrophils")
#' 
#'
#'}
#'

brooks12<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				gamSChemo=gamS+hS*C
				etaNPChemo=etaNP-hNP*C
				gamSv=gamMinS+(gamSChemo-gamMinS)*bS/(bS+G)
				gam0v=gamMin0+(gam0-gamMin0)*bn/(bn+G)
				etaNPv=etaNPChemo+(etaMaxNP-etaNPChemo)*G/(cn+G)
				tauNMv=tauNMmax/(1+(Vmax-1)*G/(bv+G))
				tauNv=tauNMv+tauNP
				I0c=0  # these here since both for negative times and when not infusing
				I0g=0
				if (Time < 0) {
					delEtaNP=etaNP*tauNP
					delGam0=gam0*tauNMmax
					delGamS=gamS*tauS
					An=exp(delEtaNP - delGam0)
					Aq=2*exp(-delGamS)
					dQ=-(k0/(1+(Q/the2)^s2) + f0/(1+(N/the1)^s1)+kdel)*Q + Aq*k0/(1+(Qss/the2)^s2)*Qss
					dN=-gamN*N + An*f0/(1+(Nss/the1)^s1)*Qss
					dC=0;	dX=0;	dG=0; trt=FALSE
				}	else {
					trt=(Time%/%T  < cycles)
					if (trt&&(Time%%T < Delc)) I0c=Dc/Delc # in chemo infusion period
					if (trt&&(Time%%T > T1)&&(Time%%T < T1+Delg) ) I0g=Dg/Delg #in G-CSF injection period
					lagM<-lagvalue(Time - tauNMv)
					lagN<-lagvalue(Time - tauNv)
					lagS<-lagvalue(Time - tauS)
					delEtaNP=lagM[3]-lagN[3]
					delGam0=Gam0 -lagM[4]
					delGamS=GamS -lagS[5]
					An=exp(delEtaNP - delGam0)
					Aq=2*exp(-delGamS)
					Qts=lagS[1]
					Qtn=lagN[1]
					Ntn=lagN[2]
					dQ=-(k0/(1+(Q/the2)^s2) + f0/(1+(N/the1)^s1)+kdel)*Q + Aq*k0/(1+(Qts/the2)^s2)*Qts
					dN=-gamN*N + An*f0/(1+(Ntn/the1)^s1)*Qtn
					dC= I0c/phi - del*C
					dX= I0g + kT*VB*G -kB*X
					dG= Gprod - kT*G + kB*X/VB -gamG*G - sig*N*G^2/(kG+G^2)
				}
				dEtaNP=etaNPv
				dGam0=gam0v
				dGamS=gamSv
				list(c(dQ,dN,dEtaNP,dGam0,dGamS,dG,dX,dC),
						c(delEtaNP=delEtaNP,delGam0=delGam0,delGamS=delGamS,
								etaNPv=etaNPv,gam0v=gam0v,gamSv=gamSv,trt=trt,tauNMv=tauNMv,tauNv=tauNv,An=An,Aq=Aq))
			})
}



