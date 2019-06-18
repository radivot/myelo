############################ craig16 Model ##################
library(tidyverse)
library(deSolve)
library(myelo)
craigPars16
# x0=c(Q=craigPars16[["Qss"]],Nr=craigPars16[["Nrss"]],N=craigPars16[["Nss"]],Gs=0,
#      G1=craigPars16[["G1ss"]],G2=craigPars16[["G2ss"]],
#      Tn=craigPars16[["tauNP"]]+craigPars16[["aNM"]],  #here time Tn is to reserves
#      An=craigPars16[["ANss"]],Cp=0,Cf=0,Cs1=0,Cs2=0,Aq=craigPars16[["AQss"]],)

x0=c(Q=craigPars16[["Qss"]],Nr=craigPars16[["Nrss"]],N=craigPars16[["Nss"]],
     G1=craigPars16[["G1ss"]],G2=craigPars16[["G2ss"]],
     Tn=craigPars16[["tauNP"]]+craigPars16[["aNM"]],  #here time Tn is to reserves
     An=craigPars16[["ANss"]])
# x0["An"]=103780 no better
times <- seq(-50,200,by=1)
yout <- dede(x0,times = times, func = craig16,	parms = craigPars16)

tail(yout,1)  # G1 and G2 are higher than in x0 and this makes Nr lower
x0[]=tail(yout,1)[-1]
yout <- dede(x0,times = times, func = craig16,	parms = craigPars16)


D=data.frame(yout)
head(D)
d=D%>%select(time:N)%>%gather(key="Cell",value="Counts",-time)%>%
  mutate(Cell=factor(Cell,levels=c("Q","Nr","N")))
tc=function(sz) theme_classic(base_size=sz)
gx=xlab("Days")
sbb=theme(strip.background=element_blank())
d%>%ggplot(aes(x=time,y=Counts))+facet_grid(Cell~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb
ggsave("~/Results/myelo/craig16ss.pdf",height=6,width=6.5)



# Section 4 of 2016 paper: Parameter Estimation and Equation Constraints
rm(list=ls())
library(tidyverse)
library(deSolve)
library(myelo)
attach(as.list(craigPars16))
search()
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


