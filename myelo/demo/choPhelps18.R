# Yu Kyoung Cho ...  Mitch A. Phelps
# CPT Pharmacometrics Syst. Pharmacol. (2018) 7, 748-758
library(tidyverse)
library(myelo)
library(mrgsolve)
(MTT=106/24)
(ktr=4/MTT)
(kcirc=log(2)/(7/24))
Circ0=5.69
(prol_0=kcirc*Circ0/ktr)
(ITT=14/24)
(kin=1/ITT)
# Melphlan PK from their 2017 paper CLINICAL PHARMACOLOGY & THERAPEUTICS | VOLUME 102 NUMBER 3 | SEPTEMBER 2017
CL=26.5 # L/hr
V1=17.7 # L
Q=28.3 # L/hr
V2=18.3 # L
(k10=CL*24/V1) # 35.9322
(k12=24*Q/V1)
(k21=k12*V1/V2)


code='
$PROB CPharmacokinetic-Pharmacodynamic Model of Neutropenia in Patients 
With Myeloma Receiving High- Dose Melphalan for Autologous Stem Cell Transplant 
Yu Kyoung Cho ...  Mitch A. Phelps
CPT Pharmacometrics Syst. Pharmacol. (2018) 7, 748-758
slope is in 1/(ug/mL), so mw conversion to uM not needed

$PARAM Circ0=5.69,ktr=0.9056604,kcirc=2.376505,  gam=0.221,slope=11.3
k12=38.37288,k21=37.11475,k10=35.9322,V1=17.7, kin=1.714286

$INIT C1=0,C2=0,Prol=14.93089,Trans1=14.93089,Trans2=14.93089
      Trans3=14.93089,Circ=5.69,Marg=0

$ODE 
double Cp=C1/V1;
double Edrug=slope*Cp;
double FB=1;

FB=Circ>=0.001 ? pow(Circ0/Circ,gam): pow(Circ0/0.001,gam);

dxdt_C1=-(k12+k10)*C1 + k21*C2;
dxdt_C2=k12*C1 - k21*C2;
dxdt_Prol = ktr*Prol*(1-Edrug)*FB - ktr*Prol;
dxdt_Trans1=ktr*Prol-ktr*Trans1;
dxdt_Trans2=ktr*Trans1-ktr*Trans2;
dxdt_Trans3=ktr*Trans2-ktr*Trans3;
dxdt_Circ=ktr*Trans3-kcirc*Circ+kin*Marg;
dxdt_Marg=-kin*Marg;
$CAPTURE FB
'
mod <- mread("choPhelps18", "~/tmp", code)
(e1=ev(time=5,amt=200*1.8,cmt=1)) 
out=mod%>%ev(e1)%>%mrgsim(start=0,end = 50, delta = 1)
out%>%plot(xlab="Days")

(e2=ev(time=15,amt=50,cmt=8,evid=8)) 
out=mod%>%ev(e1+e2)%>%mrgsim(start=0,end = 50, delta = 1)
out%>%plot(xlab="Days")


# tail(out)
# d=as.data.frame(out)
# D=d%>%select(time,Prol:Circ)%>%gather(key="Cell",value="Value",-time)%>%mutate(Cell=as_factor(Cell))
# gx=xlab("Days")
# sbb=theme(strip.background=element_blank())
# tc=function(sz) theme_classic(base_size=sz)
# D%>%ggplot(aes(x=time,y=Value,col=Cell))+geom_line(size=1)+gx+tc(14)+sbb
# ggsave("~/GitHubs/myelo/docs/waves.png",width=5, height=4)

