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
#' yout <- dede(c(Q=craigPars16[["Qss"]],Nr=craigPars16[["Nrss"]],N=craigPars16[["Nss"]]),
#' 		times = times, func = craig16,	parms = craigPars16)
#' plot(yout)
#'}
#'
# craig16<-function(Time, State, Pars) { # no chemo, no  GCSF subQ
#   with(as.list(c(State, Pars)), {
#     fbeta=function(Q) fQ/(1+(Q/the2)^s2)
#     fkapN=function(G1) kapss + (kapss-kapmin)*(G1^s1-G1ss^s1)/(G1^s1+G1ss^s1)
#     fetaNP=function(G1) etaNPss  + (etaNPss -etaNPmin)*(bNP/G1ss)*(G1-G1ss)/(G1+bNP)
#     Vn=function(G1) 1 + (Vmax-1)*(G1-G1ss)/(G1-G1ss+bV)
#     GBF=G2/(V*(Nr+N))
#     GBFss=G2ss/(V*(Nrss+Nss))
#     phiNr=function(GBF) phiNrss  + (phiNrmax-phiNrss)*(GBF-GBFss)/(GBF-GBFss+bG)
#     beta=fbeta(Q)
#     kapG1=fkapN(G1)
#     if (Time < 0) {
#       dQ=-(beta+kapG1+kapDel)*Q + AQss*beta*Q
#       dNr=-gamNr*Nr-phiNr(GBF)*Nr + (An/1e3)*kapG1*Q  # 1e3 maps units of Q to N
#       dN= phiNr(GBF)*Nr-gamNss*N
#       # dG1=0  # 4  unbound circulating GCSF
#       # dG2=0   #5  #bound GCSF
#       dG1=Gprod - kren*G1 - k12g*((Nr+N)*V+G2)*G1^Pow+k21g*G2
#       dG2=k12g*((Nr+N)*V+G2)*G1^Pow -kint*G2 - k21g*G2
#       dTn=0  #6
#       dAn=0
#     }	else {
#       Qst=lagvalue(Time - tauS,1)
#       betast=fbeta(Qst)
#       dQ=-(beta+kapG1+kapDel)*Q + AQss*betast*Qst
#       tauNM=Tn-tauNP  # Tn now time to reserves
#       G1nmt=lagvalue(Time - tauNM,4)
#       Vrat=Vn(G1)/Vn(G1nmt)
#       dTn=1-Vrat
#       G1nt=lagvalue(Time - Tn,4)
#       Qnt=lagvalue(Time - Tn,1)
#       kapG1nt=fkapN(G1nt)
#       dNr=(An/1e3)*kapG1nt*Qnt*Vrat-gamNr*Nr-phiNr(GBF)*Nr # 1e3 maps units of S to N
#       dN= phiNr(GBF)*Nr-gamNss*N
#       dG1=Gprod - kren*G1 - k12g*((Nr+N)*V+G2)*G1^Pow+k21g*G2
#       dG2=k12g*((Nr+N)*V+G2)*G1^Pow -kint*G2 - k21g*G2
#       etaNPnmt=fetaNP(G1nmt)
#       etaNPnt=fetaNP(G1nt)
#       dAn=An*((1-dTn)*(etaNPnmt-etaNPnt)-gamNMss*dTn)
#     }
#     list(c(dQ,dNr,dN,dG1,dG2,dTn,dAn))
#   })
# }

############################################

craig16<-function(Time, State, Pars) {  # with Chemo and GCSF subQ
  with(as.list(c(State, Pars)), {
    fbeta=function(Q) fQ/(1+(Q/the2)^s2)
    # fkapN=function(G1) kapss + (kapss-kapmin)*(G1^s1-G1ss^s1)/(G1^s1+G1ss^s1)
    fkapN=function(G1) kapss #+ (kapss-kapmin)*(G1^s1-G1ss^s1)/(G1^s1+G1ss^s1)
    fetaNP=function(G1) etaNPss  + (etaNPss -etaNPmin)*(bNP/G1ss)*(G1-G1ss)/(G1+bNP)
    # next func reduces to above when Cp=0 (no chemo)
    # fetaNPchemo=function(G1,Cp) etaNPInf  + (fetaNP(G1)-etaNPInf)/(1+(Cp/EC50)^Sc)
    # fetaNPchemo=function(G1,Cp) etaNPInf  + (fetaNP(G1)-etaNPInf)/(1+(Cp/EM50)^Sc)
    fetaNPchemo=function(G1,Cp) fetaNP(G1)/(1+(Cp/EM50)^Sc)
    Vn=function(G1) 1 + (Vmax-1)*(G1-G1ss)/(G1-G1ss+bV)
    GBF=G2/(V*(Nr+N))
    GBFss=G2ss/(V*(Nrss+Nss))
    phiNr=function(GBF) phiNrss  + (phiNrmax-phiNrss)*(GBF-GBFss)/(GBF-GBFss+bG)
    beta=fbeta(Q)
    kapG1=fkapN(G1)
    eta=fetaNPchemo(G1,0)
    
    dIc=0 #Constant chemo infusion rate   
    dIg=0 # same for GCSF.  Changes only by events
    
    if (Time < 0) {
      dQ=-(beta+kapG1+kapDel)*Q + Aqss*beta*Q
      dNr=-gamNrss*Nr-phiNr(GBF)*Nr + (An/1e3)*kapG1*Q  # 1e3 maps units of Q to N
      dN= phiNr(GBF)*Nr-gamNss*N
      # dG1=0  # 4  unbound circulating GCSF
      # dG2=0   #5  #bound GCSF
      dG1=Gprod - kren*G1 - k12g*((Nr+N)*V-G2)*G1^Pow+k21g*G2
      dG2=k12g*((Nr+N)*V-G2)*G1^Pow -kint*G2 - k21g*G2
      dTn=0  #6
      dAn=0  #7
      dAq=0  #8
      dCp=0  #9
      dCf=0
      dCs1=0
      dCs2=0
      dGs=0 #skin pool that feeds into G1
      
      Qtn=Qts=Qss
      Cpts=Cptn=Cptnm=0
      G1tn=G1tnm=tgam=0
      etaTnm=etaTn=eta
    }	else {
      
      # Qst=lagvalue(Time - tauS,1)
      # betast=fbeta(Qst)
      Qts=lagvalue(Time - tauS,1)
      betats=fbeta(Qts)
      # dQ=-(beta+kapG1+kapDel)*Q + Aq*betast*Qst
      dQ=-(beta+kapG1+kapDel)*Q + Aq*betats*Qts
      
      tauNM=Tn-tauNP  # Tn now time to reserves
      # G1nmt=lagvalue(Time - tauNM,4)
      G1tnm=lagvalue(Time - tauNM,4)
      # Vrat=Vn(G1)/Vn(G1nmt)
      Vrat=Vn(G1)/Vn(G1tnm)
      dTn=1-Vrat
      
      # G1nt=lagvalue(Time - Tn,4)
      # Qnt=lagvalue(Time - Tn,1)
      G1tn=lagvalue(Time - Tn,4)
      Qtn=lagvalue(Time - Tn,1)
      # kapG1nt=fkapN(G1nt)
      kapG1tn=fkapN(G1tn)
      dNr=(An/1e3)*kapG1tn*Qtn*Vrat-gamNrss*Nr-phiNr(GBF)*Nr # 1e3 maps units of S to N
      # dNr=(An/1e3)*kapG1nt*Qnt*Vrat-gamNrss*Nr-phiNr(GBF)*Nr # 1e3 maps units of S to N
      
      dN = phiNr(GBF)*Nr-gamNss*N
      
      dG1=Gprod + ka*Gs + Ig - kren*G1 - k12g*((Nr+N)*V-G2)*G1^Pow + k21g*G2
      dG2=k12g*((Nr+N)*V-G2)*G1^Pow - kint*G2 - k21g*G2
      
      # Cpnt=lagvalue(Time - Tn,9)
      # Cpnmt=lagvalue(Time - tauNM,9)
      Cptn=lagvalue(Time - Tn,9)
      Cptnm=lagvalue(Time - tauNM,9)
      
      # etaNPnt=fetaNPchemo(G1nt,Cpnt)   # etaNPnt=fetaNP(G1nt)
      # etaNPnmt=fetaNPchemo(G1nmt,Cpnmt)     # etaNPnmt=fetaNP(G1nmt)
      etaTn=fetaNPchemo(G1tn,Cptn)   # etaNPnt=fetaNP(G1nt)
      etaTnm=fetaNPchemo(G1tnm,Cptnm)     # etaNPnmt=fetaNP(G1nmt)
      # dAn=An*((1-dTn)*(etaNPnmt-etaNPnt)-gamNMss*dTn)
      dAn=An*((1-dTn)*(etaTnm-etaTn)-gamNMss*dTn)
      
      # Cpst=lagvalue(Time - tauS,9)
      Cpts=lagvalue(Time - tauS,9)
      dAq=Aq*hQ*(Cpts-Cp)
      
      dCp=k21*Cf+k31*Cs1-(k12+k13+kelC)*Cp  + Ic
      dCf=k12*Cp+k42*Cs2-(k21+k24)*Cf
      dCs1=k13*Cp-k31*Cs1
      dCs2=k24*Cf-k42*Cs2
      dGs=-ka*Gs # G in skin, add F*Dg/Vd to Gs at each injection
    }
    
    list(c(dQ,dNr,dN,dG1,dG2,dTn,dAn,dAq,dCp,dCf,dCs1,dCs2,dGs,dIc,dIg),
         c(ANC=N*8.19,Qts=Qts,Cpts=Cpts,Qtn=Qtn,G1tn=G1tn,Cptn=Cptn,G1tnm=G1tnm,
           Cptnm=Cptnm,dAn=dAn,dTn=dTn,eta=eta,etaTnm=etaTnm,etaTn=etaTn,tGam=dTn*gamNMss))
  })
}


