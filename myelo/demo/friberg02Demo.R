############################  Friberg02 ######################
library(tidyverse)
library(deSolve)
library(myelo)
fribergPars02
# Figure 4 upper left panel (docetaxel), want nadir at 7 days and repeak at 20 day
(d=data.frame(x=c(0,6,9,13,16,20,22),y=c(5.4,2,0.6,1,1.7,7,8.6))) # mouse-balling it
plot(d)
times <- seq(-5,25,0.1)
200/0.808/7.4
(evnt=data.frame(var="C1",time=0,value=33.5,method="rep"))
x0=c(C1=0,C2=0,C3=0,Prol=5.05,Trans1=5.05,Trans2=5.05,Trans3=5.05,Circ=5.05)
yout=ode(x0,times=times,func=friberg02,events=list(data=evnt),parms=fribergPars02)
plot(yout) #use export in gui
D=as.data.frame(yout)
gx=xlab("Days")
sbb=theme(strip.background=element_blank())
tc=function(sz) theme_classic(base_size=sz)
head(D)
d=D%>%select(time,C1:Circ)%>%gather(key="Lab",value="Value",-time)
d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb#+cc
ggsave("~/Results/myelo/fri02.png",width=5, height=6)


tc=function(sz) theme_classic(base_size=sz)
gx=xlab("Days")
gy=ylab(quote(paste(10^3," Neutrophils/uL")))
D%>%ggplot(aes(x=time,y=Circ))+geom_line()+gx+gy+tc(14)+ylim(0,7)
ggsave("~/Results/myelo/friberg02.png",width=5,height=6)
D$Circ
min(D$Circ)


(f=file.path(system.file(paste("libs",Sys.getenv("R_ARCH"),sep=""), package = "myelo"),
             paste("myelo",.Platform$dynlib.ext,sep="")))
dyn.load(f)
yout=ode(x0,times=times,func="derivsFri02",
         dllname = "myelo",initfunc = "parmsFri02",
         ,events=list(data=evnt),parms=fribergPars02)
plot(yout) 



library(mrgsolve)
code='
$PARAM Circ0=5.05,ktr=1.0823,gam=0.161,slope=8.58,k12=25.44,k21=36.24,k13=30.24,k31=2.016,k10=124.8
$INIT C1=0,C2=0,C3=0,Prol=5.05,Trans1=5.05,Trans2=5.05,Trans3=5.05,Circ=5.05 
$ODE 
double Edrug=slope*C1;
dxdt_C1=-(k12+k13+k10)*C1 + k21*C2 + k31*C3;
dxdt_C2=k12*C1 - k21*C2;
dxdt_C3=k13*C1 - k31*C3;
dxdt_Prol = ktr*Prol*(1-Edrug)*pow(Circ0/Circ,gam) - ktr*Prol;
dxdt_Trans1=ktr*Prol-ktr*Trans1;
dxdt_Trans2=ktr*Trans1-ktr*Trans2;
dxdt_Trans3=ktr*Trans2-ktr*Trans3;
dxdt_Circ=ktr*Trans3-ktr*Circ;
'
mod <- mread("fri02", "~/tmp", code)
(e=ev(time=5,amt=33.5,cmt=1)) #adds (evid=1) 33.5 to cmt 1 at t=5 days
mod%>%ev(e)%>%mrgsim(end = 30, delta = 0.1)%>%plot(xlab="Days")


