#'Neutrophil-GCSF-Chemo model of Craig et al 2015
#'
#'This function returns the right hand side of an ordinary differential delay
#'equation model of Craig et al. JTB 2015. 
#'The intended use of this function is as an argument to \code{dede()} of the \code{deSolve} package.
#'
#'In this model Q = stem cells, Nr = marrow reserve neutrophils, and N = circulating neutrophils
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
    # # etaNPChemo=etaNP-hNP*C
    # tauNMv=tauNMmax/(1+(Vmax-1)*G/(bv+G))
    # tauNv=tauNMv+tauNP
    # I0c=0  # these here since both for negative times and when not infusing
    # I0g=0
    # An=Anmax*exp(-etaN*tauNM)
    
    fbeta=function(Q) fQ/(1+(Q/the2)^s2)
    fkapN=function(N) fN/(1+(N/the1)^s1)
    beta=fbeta(Q) 
    kapN=fkapN(N) 
    
    
    
    # gamSv=gamSmin+(gamSChemo-gamSmin)*bS/(bS+G)
    # gam0v=gamMin0+(gam0-gamMin0)*bn/(bn+G)
    # etaNPv=etaNPChemo+(etaMaxNP-etaNPChemo)*G/(cn+G)
    fluxin= ANss*kapN*Q/1e3  
    fluxout=ftransss*Nr  
    fluxdeath=gamNr*Nr  
  
    aux=c(fluxin=fluxin,fluxout=fluxout,fluxdeath=fluxdeath)  
    
    
    if (Time < 0) {
      # As=2*exp(-gamS*tauS)
      # An=exp(etaNPss*tauNP-gamNr*tauNr)
      dQ=-(beta+kapN+kapDel)*Q + AQss*beta*Q
      dNr=-gamNr*Nr-ftransss*Nr + (ANss/1e3)*kapN*Q  # 1e3 maps units of Q to N
      # dNr=-gamNr*Nr-ftransss*Nr + ANss*kapN*Q  # 1e3 maps units of Q to N
      dN= ftransss*Nr-gamNss*N 
      dGs=0
      dG=0   #5
      dTn=0  #6
      dCp=0  #7
      dCf=0
      dCs1=0
      dCs2=0
      dAq=0
      dAn=0
    }	else {
      Qst=lagvalue(Time - tauS,1)
      Qnt=lagvalue(Time - Tn,1)
      Nnt=lagvalue(Time - Tn,3)
      Gst=lagvalue(Time - tauS,5)
      Cpst=lagvalue(Time - tauS,7)
      tauNM=Tn-tauNr-tauNP
      Gnmt=lagvalue(Time - tauNM,5)
      Cpnmt=lagvalue(Time - tauNM,7)
      Gnt=lagvalue(Time - Tn,5)
      Cpnt=lagvalue(Time - Tn,7)
      
      kapNnt=fkapN(Nnt)
      betast=fbeta(Qst)
    
      Vn=function(G) 1 + (Vmax-1)*(G-Gss)/(G-Gss+bV)
      Vrat=Vn(G)/Vn(Gnmt)
      dTn=1-Vrat
      fgamS=function(Cp,G) {
        gamSchemo=gamSss+hS*Cp
        gamSmin-(gamSmin-gamSchemo)*bS/(G-Gss+bS)
      }
      gamS=fgamS(Cp,G)
      gamSt=fgamS(Cpst,Gst)
      dAq=Aq*(gamSt-gamS)
      
      fetaNP=function(Cp,G) {
        etaNPChemo=etaNPss/(1+(Cp/EC50)^h)
        etaNPChemo + (etaNPmax-etaNPChemo)*(G-Gss)/(G-Gss+bNP)
      }
      etaNPnmt=fetaNP(Cpnmt,Gnmt) 
      etaNPnt=fetaNP(Cpnt,Gnt) 
      fgamNM=function(G)  gamNMmin - (gamNMmin-gamNMss)*bNM/(G-Gss+bNM)
      gamNM=fgamNM(G) 
      gamNMnmt=fgamNM(Gnmt) 
      dAn=An*((1-dTn)*(etaNPnmt+gamNMnmt-etaNPnt)-gamNM)
      
      
      dQ=-(beta+kapN+kapDel)*Q + Aq*betast*Qst
      dNr=-gamNr*Nr-ftransss*Nr + (An/1e3)*kapNnt*Qnt*Vrat # 1e3 maps units of S to N
      # dNr=-gamNr*Nr-ftransss*Nr + An*kapNnt*Qnt*Vrat # 
      transRatio=transmax/ftransss
      ftranst=ftransss*(transRatio*(G-Gss)+bG)/(G-Gss+bG)
      dN= ftranst*Nr-gamNss*N + An*kapNnt*Qnt
      dGs=-ka*Gs # G in skin, add F*Dg/Vd to Gs at each injection
      # chi=Gss/Nss  # already in pars
      dG=ka*Gs + Gprod - kren*G - chi*kint*(G/Kd)^2/(1+(G/Kd)^2)
      dCp=k21*Cf+k31*Cs1-(k12+k13+kelC)*Cp  # add Dc/Vd 
      dCf=k12*Cp+k42*Cs1-(k21+k24)*Cf
      dCs1=k13*Cp-k31*Cs1
      dCs2=k24*Cf-k42*Cs2
    }
    list(c(dQ,dNr,dN,dGs,dG,dTn,dCp,dCf,dCs1,dCs2,dAq,dAn),aux)
  })
}



