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
#' \bold{293} 111–120 (2012).
#'@keywords IO
#'@export
#'@examples
#'
#'\dontrun{
#'library(myelo)
#' times <- c(-zhugePars[["tauN"]]:600)
#'yout <- dede(c(Q=zhugePars[["Qss"]],N=zhugePars[["Nss"]]), times = times, func = zhuge12,	parms = zhugePars)
#'plot(yout)
#'
#'}
#'

zhuge12<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				Aq=2*exp(-gamS*tauS)
				An=exp(etaNP*tauNP-gam0*tauNM)
				if (Time < 0) {
					dQ=-(k0/(1+(Q/the2)^s2) + f0/(1+(N/the1)^s1))*Q + Aq*k0/(1+(Qss/the2)^s2)*Qss
					dN=-gamN*N + An*f0/(1+(Nss/the1)^s1)*Qss
				}
				else{
					Qts=lagvalue(Time - tauS)[1]
					Qtn=lagvalue(Time - tauN)[1]
					Ntn=lagvalue(Time - tauN)[2]
					dQ=-(k0/(1+(Q/the2)^s2) + f0/(1+(N/the1)^s1))*Q + Aq*k0/(1+(Qts/the2)^s2)*Qts
					dN=-gamN*N + An*f0/(1+(Ntn/the1)^s1)*Qtn
				}
				list(c(dQ,dN))
			})
}


