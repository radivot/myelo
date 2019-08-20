############################  Friberg02 ######################
library(tidyverse)
library(deSolve)
library(myelo)
library(mrgsolve)
(MTT=83.9/24)
(ktr=4/MTT)

code='
$PARAM Neu0=5.3,ktr=1.144219,gam=0.144,slope=14.5
k12=25.44,k21=36.24,k13=30.24,k31=2.016,k10=124.8,V1=7.4,mw=0.808
$INIT C1=0,C2=0,C3=0,Prol=5.3,Trans1=5.3,Trans2=5.3,Trans3=5.3,Neu=5.3 
$ODE 
double Cp=C1/V1/mw;
double Edrug=slope*Cp;
dxdt_C1=-(k12+k13+k10)*C1 + k21*C2 + k31*C3;
dxdt_C2=k12*C1 - k21*C2;
dxdt_C3=k13*C1 - k31*C3;
dxdt_Prol = ktr*Prol*(1-Edrug)*pow(Neu0/Neu,gam) - ktr*Prol;
dxdt_Trans1=ktr*Prol-ktr*Trans1;
dxdt_Trans2=ktr*Trans1-ktr*Trans2;
dxdt_Trans3=ktr*Trans2-ktr*Trans3;
dxdt_Neu=ktr*Trans3-ktr*Neu;
'
mod <- mread("kloft06", "~/tmp", code)
(e=ev(time=5,amt=180,cmt=1)) 
out=mod%>%ev(e)%>%mrgsim(start=0,end = 40, delta = .1)
out%>%plot(xlab="Days")
d=as.data.frame(out)
D=d%>%select(time,Prol:Neu)%>%gather(key="Cell",value="Value",-time)%>%mutate(Cell=as_factor(Cell))
gx=xlab("Days")
sbb=theme(strip.background=element_blank())
tc=function(sz) theme_classic(base_size=sz)
D%>%ggplot(aes(x=time,y=Value,col=Cell))+geom_line(size=1)+gx+tc(14)+sbb
# ggsave("~/GitHubs/myelo/docs/waves.png",width=5, height=4)


