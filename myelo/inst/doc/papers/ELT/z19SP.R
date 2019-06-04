############################ ZHUGE19 SP Model ##################
library(tidyverse)
library(deSolve)
library(myelo)  
(old=zhugePars19)

zhuge19SP<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    if (Time < 0) {
      # etaN=etaNbar/(1+(N/the1script)^nu1)
      # etaP=etaPbar/(1+(P/the4script)^nu4)
      # An=Anmax*exp(-etaN*tauNM)
      # Ap=Apmax*exp(-etaP*tauPM)
      # beta=k0/(1+(S/the2)^s2)
      # Kn=f0/(1+(N/the1)^s1)
      # Kp=Kpbar/(1+Kp*P^s4)
      # dS=-(beta+Kn+Kp+kR)*S + Aq*beta*Sss
      # dN=-gamN*N + An*Kn*Sss
      # dP=-gamP*P + Ap*(1-exp(-gamP*tauPS))*Kp*Sss 
      # dP=-gamP*P + Ap*(Kp*Sss-exp(-gamP*tauPS)*Kp*Sss)
      dS=0
      # dN=0
      dP=0
      # dP=Ap*(Kp*Sss-exp(-gamP*tauPS)*Kp*Sss)
    }	else {
      St=lagvalue(Time - tauS,1)
      Snt=lagvalue(Time - tauNM,1)
      Spt=lagvalue(Time - tauPM,1)
      Sptsum=lagvalue(Time - (tauPM+tauPS),1)
      Pt=lagvalue(Time - tauPM,2)
      Ptsum=lagvalue(Time - (tauPM+tauPS),2)
      # etaNt=etaNbar/(1+(Nss/the1script)^nu1)
      etaPt=etaPbar/(1+(Pt/the4script)^nu4)
      etaPtsum=etaPbar/(1+(Ptsum/the4script)^nu4)
      As=2*exp(-gamS*tauS)
      # Ant=Anmax*exp(-etaNt*tauNM)
      Apt=Apmax*exp(-etaPt*tauPM)
      Aptsum=Apmax*exp(-etaPtsum*tauPM)
      Kn=f0/(1+(Nss/the1)^s1)
      # Knt=f0/(1+(Nss/the1)^s1)
      Kp=Kpbar/(1+Kp*P^s4)
      Kpt=Kpbar/(1+Kp*Pt^s4)
      Kptsum=Kpbar/(1+Kp*Ptsum^s4)
      beta=k0/(1+(S/the2)^s2)
      betat=k0/(1+(St/the2)^s2)
      dS=-(beta+Kn+Kp+kR)*S + As*betat*St
      # dN=-gamN*Nss + Ant*Knt*Snt
      # dP=-gamP*P + Apt*(Kpt*Spt-exp(-gamP*tauPS)*Kptsum*Sptsum) 
      dP=-gamP*P + Apt*Kpt*Spt - exp(-gamP*tauPS)*Aptsum*Kptsum*Sptsum 
    }
    list(c(dS,dP))
  })
}

(SS<-zhugePars19[c("Sss","Pss")])/1e6
# Sss     Nss     Pss 
# 1.1   690.0 31071.0 

############################ ZHUGE19 SNP Figure 6 ##################
Kpi=zhugePars19["Kp"]
tauPMi=zhugePars19["tauPM"]

# C
zhugePars19["Kp"]=exp(-1.3)*Kpi
zhugePars19["tauPM"]=exp(0.6)*tauPMi

times <- seq(-(zhugePars19["tauPM"]+zhugePars19["tauPS"]),500,by=0.1)
yout <- dede(c(S=zhugePars19[["Sss"]],P=zhugePars19[["Pss"]]),
             times = times, func = zhuge19SP,	parms = zhugePars19)
plot(yout)

# D
zhugePars19["Kp"]=exp(2.5)*Kpi
zhugePars19["tauPM"]=exp(0.1)*tauPMi

times <- seq(-(zhugePars19["tauPM"]+zhugePars19["tauPS"]),500,by=0.1)
yout <- dede(c(S=zhugePars19[["Sss"]],P=zhugePars19[["Pss"]]),
             times = times, func = zhuge19SP,	parms = zhugePars19)
plot(yout)

# E
zhugePars19["Kp"]=exp(1.1)*Kpi
zhugePars19["tauPM"]=exp(0.2)*tauPMi

times <- seq(-(zhugePars19["tauPM"]+zhugePars19["tauPS"]),500,by=0.1)
yout <- dede(c(S=zhugePars19[["Sss"]],P=zhugePars19[["Pss"]]),
             times = times, func = zhuge19SP,	parms = zhugePars19)
plot(yout)


