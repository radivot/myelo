#'Neutrophil dynamics model of Zhuge et al 2012
#'
#'This function returns the right hand side of an ordinary differential delay
#'equation model of neutrophil dynamics published by C. Zhuge et al. in 2012 
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
#'@references  Changjing Zhuge, Jinzhi Lei, and Michael C. Mackey,
#'Neutrophil dynamics in response to chemotherapy and G-CSF, \emph{Journal of Theoretical Biology} 
#' \bold{293} 111-120 (2012).
#'@keywords IO
#'@export
#'@examples
#'
#'\dontrun{
#'library(myelo)
#' 
#' myplot=function(times, y1,y2) {
#' 	graphics.off()
#' 	windows(width=6,height=6)
#' 	par(mfrow=c(1,1),mar=c(4.7,6,3.3,1),oma=c(0,0,0,0),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
#' 	plot(times,y1/1e8,type="l",lty=2,log="y",yaxt="n",ylim=c(1e-2,1e4),ylab="",xlab="days")
#' 	lines(times,y2/1e8,lty=1)
#' 	abline(h=0.63)
#' 	axis(side=2,las=1, at=c(1e-2,1,1e2,1e4),labels=expression(10^-2,10^0,10^2,10^4))
#' 	mtext(expression(paste("Neutrophils in ",10^8," per kg")),side=2,line=3.5,cex=2)
#' }
#' 
#' times <- seq(-zhugePars[["tauN"]],200,by=0.1)
#' zhugePars["T"]=21
#' zhugePars["T1"]=4
#' yout4 <- dede(c(Q=zhugePars[["Qss"]],N=zhugePars[["Nss"]],Eta=0,Gam=0,GamS=0),
#' 		times = times, func = zhuge12,	parms = zhugePars)
#' 
#' zhugePars["T1"]=14
#' yout14 <- dede(c(Q=zhugePars[["Qss"]],N=zhugePars[["Nss"]],Eta=0,Gam=0,GamS=0),
#' 		times = times, func = zhuge12,	parms = zhugePars)
#' 
#' # this makes Figure 6B, though with a slight delay in T1=14 destabilization 
#' myplot(yout14[,1],yout14[,3],yout4[,3])
#' 
#'
#'}
#'

zhuge12<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				dEta=etaNP
				dGam=gam0
				dGamS=gamS
				tauNMv=tauNM
				tauNv=tauN
				if (Time < 0) {
					Aq=2*exp(-gamS*tauS)
					An=exp(etaNP*tauNP-gam0*tauNM)
					dQ=-(k0/(1+(Q/the2)^s2) + f0/(1+(N/the1)^s1))*Q + Aq*k0/(1+(Qss/the2)^s2)*Qss
					dN=-gamN*N + An*f0/(1+(Nss/the1)^s1)*Qss
				}	else {
					if (!is.na(T)&(Time%%T < 1)) { # in chemo
						dEta=etaMinNP
						dGamS=gamMaxS
					} 
					if ( !is.na(T1)&&(Time%%T -T1 > 0)&(Time%%T -T1 < 1) ) { # in G-CSF exposure period
						dEta=etaMaxNP
						dGam=gamMin0
						dGamS=gamMinS
						tauNMv=tauNMgcsf
						tauNv=tauNMv+tauNP
					}
					delEta=lagvalue(Time - tauNMv)[3]-lagvalue(Time - tauNv)[3]
					delGam=Gam -lagvalue(Time - tauNMv)[4]
					delGamS=GamS -lagvalue(Time - tauS)[5]
					An=exp(delEta - delGam)
					Aq=2*exp(-delGamS)
					Qts=lagvalue(Time - tauS)[1]
					Qtn=lagvalue(Time - tauNv)[1]
					Ntn=lagvalue(Time - tauNv)[2]
					dQ=-(k0/(1+(Q/the2)^s2) + f0/(1+(N/the1)^s1))*Q + Aq*k0/(1+(Qts/the2)^s2)*Qts
					dN=-gamN*N + An*f0/(1+(Ntn/the1)^s1)*Qtn
				}
				list(c(dQ,dN,dEta,dGam,dGamS))
			})
}



