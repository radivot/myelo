############ Section 3.5 of Craig et al BMB 2016: Modelling Exogenous Drug Administration
rm(list=ls())
library(tidyverse)
library(deSolve)
library(myelo)
attach(as.list(craigPars16))
search()

gcsfSQ<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    # if (Time <= 0) {
    #   dGs=0 #skin pool feed into G1
    # }	else {
    # }
      dGs=-ka750*Gs # G in skin, add F*Dg/Vd to Gs at each injection
    list(c(dGs))
  })
}

times <- seq(0,2,by=.1)
# (eventdat=data.frame(var="Gs",time=0,value=F750*750/Vd750,method="add"))
yout <- dede(c(Gs=F750*750e3/Vd750),times = times, func = gcsfSQ,	parms = craigPars16)
# yout <- dede(c(Gs=0),times = times, func = gcsfSQ,	parms = craigPars16,
             # events=list(data=eventdat),method="lsodar")
D=data.frame(yout)
head(D)
d=D%>%select(time,Gs)%>%mutate(y=F750*750e3*exp(-ka750*time)/Vd750) 
head(d)  # rate into G1 is then this times ka, i.e. Eq. 36


############# Section 4 of Craig et al BMB 2016: Parameter Estimation and Equation Constraints
fbeta=function(Q) fQ/(1+(Q/the2)^s2)
fkapN=function(G1) kapss + (kapss-kapmin)*(G1^s1-G1ss^s1)/(G1^s1+G1ss^s1)
phiNr=function(GBF) phiNrss  + (phiNrmax-phiNrss)*(GBF-GBFss)/(GBF-GBFss+bG)
GBFss=G2ss/(V*(Nrss+Nss))
# Q=Qss; G1=G1ss;Nr=Nrss;N=Nss
beta=fbeta(Qss) 
kapG1=fkapN(G1ss)

-(beta+kapG1+kapDel)*Qss + AQss*beta*Qss
(beta+kapG1+kapDel)*Qss 
AQss*beta*Qss
kapG1 # 0.0073325 = kapss=kapmin=0.0073325. By Eq. 45, this should equal
(AQss-1)*beta/3 # 0.007339529 (trust the2=0.0809 in beta or  AQss=1.5116?)
# think of (AQss-1) as surplus per cycle and beta as rate of cycling (third of profit to N)
Npss # 0.93
(kapG1/1e3)*Qss*(exp(etaNPss*tauNP)-1)/etaNPss  #Eq. 46  0.9296889 e9/Kg are dividing at SS
(kapG1/1e3)*Qss*exp(etaNPss*tauNP) # 1.547661e9/kg enter maturation each day (leave proliferation)

ANss  # 103735.4
exp(etaNPss*tauNP-gamNMss*aNM)  # 103780 #fewest sig digits in aNM=3.9
Nmss # 4.51
ANss*(kapG1/1e3)*Qss*(exp(gamNMss*aNM)-1)/gamNMss  #4.510228   Eq. 50
ANss*(kapG1/1e3)*Qss # 0.8370635 rate at which they leave maturation, 1e9/kg/day
                     # since ANss accounts for NP growth and NM killing

# Eq. 51
gamNss # 2.1875  fast killing once in blood
log(2)/(tauhalf/24) # 2.188886  redundant with 7.6 hr (.31 day) half life
(tauN=1/gamNss) # .4571429 days to bring down to exp(-1) i.e. 36% rather than 50% 

# Eq. 52
phiNr(GBFss)*Nrss  # 0.82264 rate out of reserve 
gamNss*Nss   # 0.8227187     equals rate out of N pool in blood
gamNr*Nrss   # 0.01438739    (both of which are not much less than rate in, 0.8370635 above) 
Nrss
gamNss*Nss/phiNr(GBFss)

1/(gamNr+phiNr(GBFss)) #2.700031  Eq. 53
tauNr  # 2.7, check

tauNr*phiNr(GBFss) # < 1 to ensure gamNr>0 Eq. 55

# Eq. 56
tauNr/tauN # this should be less than
Nrss/Nss  # this

Nrss*(exp(gamNMss*aNM)-1) #1.92
gamNMss*tauNr*Nmss  # 1.920313

Nrss/Nmss # should be less than
tauNr/aNM # this Eq. 58 to have apoptosis in maturation 

ANss #  103780
Nrss/((kapG1/1e3)*Qss*tauNr) # 103776.7  Eq 61

-gamNr*Nrss-phiNr(GBF)*Nrss + (ANss/1e3)*kapG1*Qss  # 1e3 maps units of Q to N
gamNr*Nrss
(ANss/1e3)*kapG1*Qss 

phiNr(GBFss)*Nrss-gamNss*Nss 
phiNr(GBFss)*Nrss
gamNss*Nss 

(dG1=Gprod - kren*G1ss - k12g*((Nrss+Nss)*V+G2ss)*G1ss^Pow+k21g*G2ss)
(dG2=k12g*((Nrss+Nss)*V+G2ss)*G1ss^Pow+k21g*G2ss -kint*G2ss - k21g*G2ss) # much bigger
(dG2=k12g*((Nrss+Nss)*V+G2ss)*G1ss^Pow -kint*G2ss - k21g*G2ss) # much bigger

k12g*((Nrss+Nss)*V+G2ss)*G1ss^Pow+k21g*G2ss -kint*G2ss - k21g*G2ss) # much bigger



detach(as.list(craigPars16))


