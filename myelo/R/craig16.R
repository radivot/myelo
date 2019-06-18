#'Neutrophil-GCSF-Chemo model of Craig et al 2016
#'
#'This function returns the right hand side of an ordinary differential delay
#'equation model of Craig et al. Bull Math Biol 2016. 
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


#'@references  # M. Craig, A. R. Humphries, M. C. Mackey, A Mathematical Model of 
#'Granulopoiesis Incorporating the Negative Feedback Dynamics and Kinetics of 
#'G-CSF/Neutrophil Binding and Internalization, \emph{Bull Math Biol} \bold{78} 2304-2357  (2016). 
#'
#'@keywords IO
#'
#'@export
#'@examples
#'
#'\dontrun{
#'library(deSolve)
#'library(myelo)
#' craigPars16
#' times <- seq(-30,200,by=.2)
#' yout <- dede(c(Q=craigPars15[["Qss"]],Nr=craigPars15[["Nrss"]],N=craigPars15[["Nss"]]),
#' 		times = times, func = craig16,	parms = craigPars16)
#' plot(yout)
#'}
#'
craig16<-function(Time, State, Pars) { # no chemo, no external GCSF given subQ
  with(as.list(c(State, Pars)), {
    fbeta=function(Q) fQ/(1+(Q/the2)^s2)
    fkapN=function(G1) kapss + (kapss-kapmin)*(G1^s1-G1ss^s1)/(G1^s1+G1ss^s1)
    fetaNP=function(G1) etaNPss  + (etaNPss -etaNPmin)*(bNP/G1ss)*(G1-G1ss)/(G1+bNP)
    Vn=function(G1) 1 + (Vmax-1)*(G1-G1ss)/(G1-G1ss+bV)
    GBF=G2/(V*(Nr+N))
    GBFss=G2ss/(V*(Nrss+Nss))
    phiNr=function(GBF) phiNrss  + (phiNrmax-phiNrss)*(GBF-GBFss)/(GBF-GBFss+bG)
    beta=fbeta(Q)
    kapG1=fkapN(G1)
    if (Time < 0) {
      dQ=-(beta+kapG1+kapDel)*Q + AQss*beta*Q
      dNr=-gamNr*Nr-phiNr(GBF)*Nr + (An/1e3)*kapG1*Q  # 1e3 maps units of Q to N
      dN= phiNr(GBF)*Nr-gamNss*N
      # dG1=0  # 4  unbound circulating GCSF
      # dG2=0   #5  #bound GCSF
      dG1=Gprod - kren*G1 - k12g*((Nr+N)*V+G2)*G1^Pow+k21g*G2
      dG2=k12g*((Nr+N)*V+G2)*G1^Pow -kint*G2 - k21g*G2
      dTn=0  #6
      dAn=0
    }	else {
      Qst=lagvalue(Time - tauS,1)
      betast=fbeta(Qst)
      dQ=-(beta+kapG1+kapDel)*Q + AQss*betast*Qst
      tauNM=Tn-tauNP  # Tn now time to reserves
      G1nmt=lagvalue(Time - tauNM,4)
      Vrat=Vn(G1)/Vn(G1nmt)
      dTn=1-Vrat
      G1nt=lagvalue(Time - Tn,4)
      Qnt=lagvalue(Time - Tn,1)
      kapG1nt=fkapN(G1nt)
      dNr=(An/1e3)*kapG1nt*Qnt*Vrat-gamNr*Nr-phiNr(GBF)*Nr # 1e3 maps units of S to N
      dN= phiNr(GBF)*Nr-gamNss*N
      dG1=Gprod - kren*G1 - k12g*((Nr+N)*V+G2)*G1^Pow+k21g*G2
      dG2=k12g*((Nr+N)*V+G2)*G1^Pow -kint*G2 - k21g*G2
      etaNPnmt=fetaNP(G1nmt)
      etaNPnt=fetaNP(G1nt)
      dAn=An*((1-dTn)*(etaNPnmt-etaNPnt)-gamNMss*dTn)
    }
    list(c(dQ,dNr,dN,dG1,dG2,dTn,dAn))
  })
}

# craig16<-function(Time, State, Pars) {
#   with(as.list(c(State, Pars)), {
#     fbeta=function(Q) fQ/(1+(Q/the2)^s2)
#     fkapN=function(G1) kapss + (kapss-kapmin)*(G1^s1-G1ss^s1)/(G1^s1+G1ss^s1)
#     fetaNP=function(G1) etaNPss  + (etaNPss -etaNPmin)*(bNP/G1ss)*(G1-G1ss)/(G1+bNP)
#     # next func reduces to above when Cp=0 (no chemo)
#     fetaNPchemo=function(G1,Cp) etaNPInf  + (fetaNP(G1)-etaNPInf)/(1+(Cp/EC50)^Sc)
#     Vn=function(G1) 1 + (Vmax-1)*(G1-G1ss)/(G1-G1ss+bV)
#     GBF=G2/(V*(Nr+N))
#     GBFss=G2ss/(V*(Nrss+Nss))
#     phiNr=function(GBF) phiNrss  + (phiNrmax-phiNrss)*(GBF-GBFss)/(GBF-GBFss+bG)
#     
#     beta=fbeta(Q) 
#     kapG1=fkapN(G1) 
#     if (Time < 0) {
#       dQ=-(beta+kapG1+kapDel)*Q + AQss*beta*Q
#       dNr=-gamNr*Nr-phiNr(GBF)*Nr + (An/1e3)*kapG1*Q  # 1e3 maps units of Q to N
#       dN= phiNr(GBF)*Nr-gamNss*N 
#       dG1=0  # 4  unbound circulating GCSF
#       dG2=0   #5  #bound GCSF
#       dTn=0  #6
#       dAn=0  #7
#       dAq=0  #8
#       dCp=0  #9
#       dCf=0
#       dCs1=0
#       dCs2=0
#       # dGs=0 #skin pool feed into G1
#     }	else {
#       
#       Qst=lagvalue(Time - tauS,1)
#       betast=fbeta(Qst)
#       dQ=-(beta+kapG1+kapDel)*Q + AQss*betast*Qst
#       
#       tauNM=Tn-tauNP  # Tn now time to reserves
#       G1nmt=lagvalue(Time - tauNM,4)
#       Vrat=Vn(G1)/Vn(G1nmt)
#       dTn=1-Vrat
#       
#       G1nt=lagvalue(Time - Tn,4)
#       Qnt=lagvalue(Time - Tn,1)
#       kapG1nt=fkapN(G1nt)
#       dNr=(An/1e3)*kapG1nt*Qnt*Vrat-gamNr*Nr-phiNr(GBF)*Nr # 1e3 maps units of S to N
#       
#       dN= phiNr(GBF)*Nr-gamNss*N 
#       
#       dG1=Gprod - kren*G1 - k12g*((Nr+N)*V+G2)*G1^Pow+k21g*G2 
#       dG2=k12g*((Nr+N)*V+G2)*G1^Pow+k21g*G2 -kint*G2 - k21g*G2
#       
#       Cpnt=lagvalue(Time - Tn,9)
#       Cpnmt=lagvalue(Time - tauNM,9)
#       
#       etaNPnt=fetaNPchemo(G1nt,Cpnt)   # etaNPnt=fetaNP(G1nt) 
#       etaNPnmt=fetaNPchemo(G1nmt,Cpnmt)     # etaNPnmt=fetaNP(G1nmt) 
#       dAn=An*((1-dTn)*(etaNPnmt-etaNPnt)-gamNMss*dTn)
# 
#       Cpst=lagvalue(Time - tauS,9)
#       dAq=Aq*hQ*(Cpst-Cp) 
#       
#       
#       
#       
#       # 
#       # 
#       # 
#       # Nnt=lagvalue(Time - Tn,3)
#       # G1st=lagvalue(Time - tauS,5)
#       # 
#       # 
#       # 
#       # 
#       # 
#       # 
#       # fgamS=function(Cp,G) {
#       #   gamSchemo=gamSss+hS*Cp
#       #   gamSmin-(gamSmin-gamSchemo)*bS/(G-Gss+bS)
#       # }
#       # gamS=fgamS(Cp,G)
#       # gamSt=fgamS(Cpst,Gst)
#       # dAq=Aq*(gamSt-gamS)
#       # 
#       # # fetaNP=function(Cp,G) {
#       # #   etaNPChemo=etaNPss/(1+(Cp/EC50)^h)
#       # #   etaNPChemo + (etaNPmax-etaNPChemo)*(G-Gss)/(G-Gss+bNP)
#       # # }
#       # etaNPnmt=fetaNP(Cpnmt,Gnmt) 
#       # etaNPnt=fetaNP(Cpnt,Gnt) 
#       # fgamNM=function(G)  gamNMmin - (gamNMmin-gamNMss)*bNM/(G-Gss+bNM)
#       # gamNM=fgamNM(G) 
#       # # gamNMnmt=fgamNM(Gnmt) 
#       
#       # dGs=0
#       
#       # transRatio=transmax/ftransss
#       # ftranst=ftransss*(transRatio*(G-Gss)+bG)/(G-Gss+bG)
#       # dN= ftranst*Nr-gamNss*N + An*kapNnt*Qnt
#       # dGs=-ka*Gs # G in skin, add F*Dg/Vd to Gs at each injection
#       # # chi=Gss/Nss  # already in pars
#       # dCp=k21*Cf+k31*Cs1-(k12+k13+kelC)*Cp  # add Dc/Vd 
#       # dCf=k12*Cp+k42*Cs1-(k21+k24)*Cf
#       # dCs1=k13*Cp-k31*Cs1
#       # dCs2=k24*Cf-k42*Cs2
#     }
#     # fluxin=(An/1e3)*kapG1nt*Qnt*Vrat
#     # fluxout=phiNr(GBF)*Nr
#     # fluxdeath=gamNr*Nr
#     
#     # aux=c(fluxin=fluxin,fluxout=fluxout,fluxdeath=fluxdeath)  
#     
#     # list(c(dQ,dNr,dN,dGs,dG1,dG2,dTn,dCp,dCf,dCs1,dCs2,dAq,dAn))
#     list(c(dQ,dNr,dN,dG1,dG2,dTn,dAn))
#     # list(c(dQ,dNr,dN,dG1,dG2,dTn,dAn),aux)
#   })
# }
# 
# 
# 
