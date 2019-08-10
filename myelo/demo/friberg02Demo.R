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
38.5/7.4
(evnt=data.frame(var="C1",time=0,value=180,method="rep"))
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
# ggsave("~/Results/myelo/fri02.png",width=5, height=6)

gy=ylab(quote(paste(10^3," Neutrophils/uL")))
D%>%ggplot(aes(x=time,y=Circ))+geom_line()+gx+gy+tc(14)+ylim(0,7)
# ggsave("~/Results/myelo/friberg02.png",width=5,height=6)
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
$PARAM Circ0=5.05,ktr=1.0823,gam=0.161,slope=8.58
k12=25.44,k21=36.24,k13=30.24,k31=2.016,k10=124.8,V1=7.4,mw=0.808
$INIT C1=0,C2=0,C3=0,Prol=5.05,Trans1=5.05,Trans2=5.05,Trans3=5.05,Circ=5.05 
$ODE 
double Cp=C1/V1/mw;
double Edrug=slope*Cp;
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
(e=ev(time=0,amt=180,cmt=1)) 
out=mod%>%ev(e)%>%mrgsim(start=0,end = 25, delta = 1)
out%>%plot(xlab="Days")
d=as.data.frame(out)
D=d%>%select(time,Prol:Circ)%>%gather(key="Cell",value="Value",-time)%>%mutate(Cell=as_factor(Cell))
D%>%ggplot(aes(x=time,y=Value,col=Cell))+geom_line(size=1)+gx+tc(14)+sbb
ggsave("~/GitHubs/myelo/docs/waves.png",width=5, height=4)


(e=ev_rx("50 q 21 x 6 then 25 q 14 x 6"))
mod%>%ev(e)%>%mrgsim(end = 300, delta = 0.1)%>%plot(xlab="Days")
(e=ev_rx("50 over 12 q 21 x 6 then 25 q 14 x 6"))
mod%>%ev(e)%>%mrgsim(end = 300, delta = 0.1)%>%plot(xlab="Days")

# cant have negative times, so forget about starting negative

sd=0.1
d$ANC=d$Circ+rnorm(dim(d)[1],sd=sd)
d%>%ggplot(aes(x=time,y=Circ))+ geom_line(size=.1)+
  geom_point(aes(x=time,y=ANC),size=1)+gx+tc(14)+sbb
ggsave("~/GitHubs/myelo/docs/noiseData",width=5, height=4)


dput(fribergPars02)
(pars=c(Circ0 = 5.05, ktr = 1.08229988726043, gam = 0.161, slope = 8.58))# recover 4 params
LF=function(pars) {
  evnt=ev(time=0,amt=180,cmt=1,Circ0=pars["Circ0"],ktr=pars["ktr"],gam=pars["gam"],slope=pars["slope"])
  out=mod%>%ev(evnt)%>%mrgsim(start=0,end = 25, delta = 1)
  out%>%plot(xlab="Days")
  as.data.frame(out)
}
D=LF(pars)%>%select(time,Circ)
dd=d%>%select(time,Circ=ANC)
library(FME)
LFcost <- function (pars) {
  out=LF(pars)%>%select(time,Circ)
  modCost(model = out, obs = dd,sd=sd)
}
(Fit <- modFit(f = PHcost, p = 1.5*pars))
pars
coef(Fit)
1.5*pars
summary(Fit)



