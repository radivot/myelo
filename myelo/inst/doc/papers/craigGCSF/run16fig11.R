############################ craig16 Model ##################
rm(list=ls())
library(tidyverse)
library(deSolve)
library(myelo)
attach(as.list(craigPars16))
kapmin
kapmin=0.0052359  # value 1 in Table 2 matches the plot
fkapN=function(G1) kapss + (kapss-kapmin)*(G1^s1-G1ss^s1)/(G1^s1+G1ss^s1)
G1=seq(0,.1,0.001)
kap=fkapN(G1)
tc=function(sz) theme_classic(base_size=sz)
tibble(G1,kap)%>%ggplot(aes(x=G1,y=kap))+geom_point(size=3,data=tibble(G1=G1ss,kap=kapss))+
  geom_line(size=1)+tc(14)
ggsave("~/Results/myelo/craig16fig11a.png",height=4,width=4)
# conclude we should try overriding kapmin to 5.2e-3

fetaNP=function(G1) etaNPss  + (etaNPss -etaNPmin)*(bNP/G1ss)*(G1-G1ss)/(G1+bNP)
eta=fetaNP(G1)
tibble(G1,eta)%>%ggplot(aes(x=G1,y=eta))+geom_point(size=3,data=tibble(G1=G1ss,eta=etaNPss))+
  geom_line(size=1)+tc(14)
ggsave("~/Results/myelo/craig16fig11b.png",height=4,width=4)
# matches fines

# D=c(0.05,5) #mg/m2 in Gonzales-Sales
# D=c(2,4) #mg/m2 given bi-weekly in Gonzales-Sales
# BSA*D*1e6/(822*1e3) #ng/mL   # Vd=882L in 
# BSA*D*1e6/(32.7*1e3) #ng/mL   # V1=32.7L in Ruixo et al 2012 
# 43.7*24 # CL 1048.8 L/day
# 43.7*24/822 # CL 1.2 bodyfulls/day
# 43.7*24/32.7 # 32
# kelC #132
# 6.39*24 # 153 in 2012 pk paper is close
# EC50 # 0.7539  # guess it should be ~25 ng/mL based on Fig11c, assuming ng/mL are units there
# 32.7*0.754
# 25/EC50  # 33-fold gap
# EC50=30
# EC50=24.6
fetaNPchemo=function(G1,Cp) etaNPInf  + (fetaNP(G1)-etaNPInf)/(1+(Cp/EM50)^Sc)
# fetaNPchemo=function(Cp) etaNPInf  + (fetaNP(G1ss)-etaNPInf)/(1+(Cp/EC50)^Sc)
Cp=seq(0,300)
# eta=fetaNPchemo(G1ss,Cp)                   #no steady state Cp level so skip making the point
eta=fetaNPchemo(0,Cp) #starting value matches G1=0, not G1=G1ss
tibble(Cp,eta)%>%ggplot(aes(x=Cp,y=eta))+ # geom_point(size=3,data=tibble(G1=G1ss,eta=etaNPss))+
  geom_line(size=1)+tc(14)+ylim(0,1.7)
ggsave("~/Results/myelo/craig16fig11c.png",height=4,width=4)


bV
Vmax
bbarV*Vmax

Vn=function(G1) 1 + (Vmax-1)*(G1-G1ss)/(G1-G1ss+bV)
G1=seq(0,1.2,0.01)
v=Vn(G1)
vss=Vn(G1ss)
tibble(G1,v)%>%ggplot(aes(x=G1,y=v))+geom_point(size=3,data=tibble(G1=G1ss,v=vss))+
  geom_line(size=1)+tc(14)
ggsave("~/Results/myelo/craig16fig11d.png",height=4,width=4) # matches fines



GBFss=G2ss/(V*(Nrss+Nss))
GBF=seq(0,0.001,0.00001)
phiNr=function(GBF) phiNrss  + (phiNrmax-phiNrss)*(GBF-GBFss)/(GBF-GBFss+bG)
phi=phiNr(GBF)
phiss=phiNr(GBFss)
tibble(GBF,phi)%>%ggplot(aes(x=GBF,y=phi))+geom_point(size=3,data=tibble(GBF=GBFss,phi=phiss))+
  geom_line(size=1)+tc(14)
ggsave("~/Results/myelo/craig16fig11e.png",height=4,width=4) # matches fine also

# need to figure out Cp units and try kapmin at 0.0052

