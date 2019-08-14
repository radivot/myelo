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
# ggsave("~/GitHubs/myelo/docs/waves.png",width=5, height=4)


(e=ev_rx("50 q 21 x 6 then 25 q 14 x 6"))
mod%>%ev(e)%>%mrgsim(end = 300, delta = 0.1)%>%plot(xlab="Days")
(e=ev_rx("50 over 12 q 21 x 6 then 25 q 14 x 6"))
mod%>%ev(e)%>%mrgsim(end = 300, delta = 0.1)%>%plot(xlab="Days")

# cant have negative times, so forget about starting negative


####### simulate data
(e=ev(time=0,amt=180,cmt=1)) 
END=50
DELTA=1
out=mod%>%ev(e)%>%mrgsim(start=0,end=END,delta=DELTA)
d=as.data.frame(out)[-1,]
d=d[!duplicated(d$time),]
sd=0.05
d$ANC=d$Circ+rnorm(dim(d)[1],sd=sd)
d%>%ggplot(aes(x=time,y=Circ))+ geom_line(size=.1)+
  geom_point(aes(x=time,y=ANC),size=1)+gx+tc(14)+sbb
# ggsave("~/GitHubs/myelo/docs/noiseData.png",width=5, height=4)


####### try FME with mrgsolve
library(FME)
dput(fribergPars02)
### start with two scale params for y and time axis
(pars=c(Circ0 = 5.05, ktr = 1.08229988726043))
LF=function(pars) {
  evnt=ev(time=0,amt=180,cmt=1,Circ0=pars["Circ0"],ktr=pars["ktr"],
          Prol_0=pars["Circ0"],Trans1_0=pars["Circ0"],Trans2_0=pars["Circ0"],
          Trans3_0=pars["Circ0"],Circ_0=pars["Circ0"] )
  as.data.frame(mod%>%ev(evnt)%>%mrgsim(start=0,end = END, delta = DELTA))[-1,]
}
D=LF(pars)%>%select(time,Circ)
dd=d%>%select(time,Circ=ANC)%>%mutate(sd=sd)
LFcost <- function (pars) {
  out=LF(pars)%>%select(time,Circ)
  modCost(model = out, obs = dd, err = "sd")
}
Pars=1.5*pars
Fit <- modFit(f = LFcost, p = Pars,method="Nelder-Mead")
Pars
coef(Fit)
pars
summary(Fit)

######### try adding in gam
dput(fribergPars02)
(pars=c(Circ0 = 5.05, ktr = 1.08229988726043,gam = 0.161))
LF3=function(pars) {
  evnt=ev(time=0,amt=180,cmt=1,Circ0=pars["Circ0"],ktr=pars["ktr"],gam = pars["gam"],
          Prol_0=pars["Circ0"],Trans1_0=pars["Circ0"],Trans2_0=pars["Circ0"],
          Trans3_0=pars["Circ0"],Circ_0=pars["Circ0"] )
  as.data.frame(mod%>%ev(evnt)%>%mrgsim(start=0,end = END, delta = DELTA))[-1,]
}
D=LF3(pars)%>%select(time,Circ)
dd=d%>%select(time,Circ=ANC)%>%mutate(sd=sd)
LF3cost <- function (pars) {
  out=LF3(pars)%>%select(time,Circ)
  modCost(model = out, obs = dd, err = "sd")
}
Pars=1.5*pars
LF3cost2 <- function(lpars)  LF3cost(exp(lpars))
Fit <- modFit(f = LF3cost2, p = log(Pars),method="Nelder-Mead") #crashes without log
exp(coef(Fit))
Pars
pars
summary(Fit)

######### try adding in slope
dput(fribergPars02)
(pars=c(Circ0 = 5.05, ktr = 1.08229988726043,gam = 0.161, slope = 8.58))
LF4=function(pars) {
  evnt=ev(time=0,amt=180,cmt=1,Circ0=pars["Circ0"],ktr=pars["ktr"],gam=pars["gam"],slope=pars["slope"],
          Prol_0=pars["Circ0"],Trans1_0=pars["Circ0"],Trans2_0=pars["Circ0"],
          Trans3_0=pars["Circ0"],Circ_0=pars["Circ0"] )
  as.data.frame(mod%>%ev(evnt)%>%mrgsim(start=0,end = END, delta = DELTA))[-1,]
}
D=LF4(pars)%>%select(time,Circ)
dd=d%>%select(time,Circ=ANC)%>%mutate(sd=sd)
LF4cost <- function (pars) {
  out=LF4(pars)%>%select(time,Circ)
  modCost(model = out, obs = dd, err = "sd")
}
parsIC=2*pars
LF4cost2 <- function(lpars)  LF4cost(exp(lpars))
Fit <- modFit(f = LF4cost2, p = log(parsIC),method="Nelder-Mead") #crashes without log
data.frame(parsIC,fit=exp(coef(Fit)),trueVals=pars)
(s=summary(Fit))
data.frame(point=exp(s$par[,1]),
           lowCI=exp(s$par[,1]-1.96*s$par[,2]),
           hiCI=exp(s$par[,1]+1.96*s$par[,2])  )



########### BBMLE approach using deSolve C
library(bbmle)
pars=fribergPars02
nLL<-function(Circ0,ktr,gam,slope) { # pass these globally: d,pars
  # attach(as.list(IC))
  Circ0=exp(Circ0)
  ktr=exp(ktr)
  gam=exp(gam)
  slope=exp(slope)
  x0=c(C1=0,C2=0,C3=0,Prol=Circ0,Trans1=Circ0,Trans2=Circ0,Trans3=Circ0,Circ=Circ0)
  pars["Circ0"]=Circ0 
  pars["ktr"]=ktr
  pars["gam"]=gam
  pars["slope"]=slope
  out=ode(x0,times=seq(0,END,DELTA),func="derivsFri02",
          dllname = "myelo",initfunc = "parmsFri02",
          ,events=list(data=evnt),parms=pars)
  y.pred=out[,"Circ"]
  sigma  <- sqrt(sum((d$ANC-y.pred)^2)/length(d$ANC))
  -sum(dnorm(d$ANC, mean=out[,"Circ"],sd=sigma,log=TRUE)) 
}

IC0=c(Circ0=5.05,ktr=1.083,gam=0.161,slope=8.58) 
IC=log(2*IC0)
(s=summary(M<-mle2(nLL,method="Nelder-Mead",
                start=as.list(IC),
                control = list(maxit=50000, parscale=IC) ) ) )
data.frame(IC=exp(IC),fit=exp(coef(M)),trueVals=IC0)
data.frame(point=exp(s@coef[,1]),
           lowCI=exp(s@coef[,1]-1.96*s@coef[,2]),
           hiCI=exp(s@coef[,1]+1.96*s@coef[,2])  )
######### END BBMLE with deSolve C

######### FME with deSolve C
pars=fribergPars02
LF4des=function(pars) {
  Circ0=pars[["Circ0"]]
  x0=c(C1=0,C2=0,C3=0,Prol=Circ0,Trans1=Circ0,Trans2=Circ0,Trans3=Circ0,Circ=Circ0)
  as.data.frame(ode(x0,times=seq(0,END,DELTA),func="derivsFri02",
          dllname = "myelo",initfunc = "parmsFri02",
          ,events=list(data=evnt),parms=pars))
}
D=LF4des(pars)%>%select(time,Circ)
dd=d%>%select(time,Circ=ANC)%>%mutate(sd=sd)
LF4desCost <- function (pars) {
  out=LF4des(pars)%>%select(time,Circ)
  modCost(model = out, obs = dd, err = "sd")
}
dput(pars)
LF4desCost2 <- function(lpars)  LF4desCost(c(exp(lpars), k12 = 25.44, k21 = 36.24, 
               k13 = 30.24, k31 = 2.016, k10 = 124.8, V1 = 7.4, mw = 0.808))
parsIC <- pars[1:4] * 1
parsIC <- pars[1:4] * 1.9 # crashes at 2, not sure why
Fit <- modFit(f = LF4desCost2, p = log(parsIC),method="Nelder-Mead") 
data.frame(parsIC,fit=exp(coef(Fit)),trueVals=pars[1:4])
(s=summary(Fit))
data.frame(point=exp(s$par[,1]),
           lowCI=exp(s$par[,1]-1.96*s$par[,2]),
           hiCI=exp(s$par[,1]+1.96*s$par[,2])  )
######### END FME with deSolve C

####### back to FME to display estimate distrubutions
ini   <- LF4(parsIC)
final <- LF4(exp(coef(Fit)))
head(d)
ini%>%ggplot(aes(x=time,y=Circ))+geom_line()+geom_point(data=d)+
  geom_line(data=final,col="red",size=1)+tc(14)+gx
# ggsave("~/Results/myelo/parsIC2xLF.png",width=5,height=5)
# ggsave("~/GitHubs/myelo/docs/parsIC2xLF.png",width=5, height=4)



var0 <- Fit$var_ms_unweighted
cov0 <- summary(Fit)$cov.scaled * 2.4^2/5
# MCMC <- modMCMC(f = LF4cost2, p = Fit$par, niter = 2000, jump = cov0,
#               var0 = var0, wvar0 = 0.1, updatecov = 50)
# save(MCMC,file="~/Results/myelo/mcmc.RData")
load("~/Results/myelo/mcmc.RData")
MCMC$pars <- exp(MCMC$pars)
summary(MCMC)
par(mar=c(4, 4, 3, 1) + .1)
plot(MCMC, Full = TRUE)
pairs(MCMC, nsample = 1000,cex.labels=1.4,cex=0.7)


######### data for sensitivity c(Circ0=5.05,ktr=1.083,gam=0.161,slope=8.58) 
DELTA=.1
(e=ev(time=0,amt=180,cmt=1)) 
Circ0=5.05
ic=c(Circ0=Circ0,Prol_0=Circ0,Trans1_0=Circ0,Trans2_0=Circ0,Trans3_0=Circ0,Circ_0=Circ0)
span=c(1/2, 1, 2)
idata=data.frame(ID=1:3,span%*%t(ic))
out=mod%>%idata_set(idata)%>%ev(e)%>%mrgsim(start=0,end=END,delta=DELTA)
fixout=function(out) {
  d=as.data.frame(out)[-1,]
  d$ID=as_factor(d$ID)
  d
}
lb=labs(col = "Circ0")
scd = scale_color_discrete(name = "Circ0", labels = idata$Circ0)
out%>%fixout%>%ggplot(aes(x=time,y=Circ,col=ID))+ geom_line(size=1)+gx+gy+tc(14)+lb+scd
# ggsave("~/GitHubs/myelo/docs/circ0Sens.png",width=5, height=4)

idata=data.frame(ID=1:3,ktr=1.083*span)
lb=labs(col = "ktr")
scd = scale_color_discrete(name = "ktr", labels = idata$ktr)
out=mod%>%idata_set(idata)%>%ev(e)%>%mrgsim(start=0,end=END,delta=DELTA)
out%>%fixout%>%ggplot(aes(x=time,y=Circ,col=ID))+ geom_line(size=1)+gx+gy+tc(14)+lb+scd
# ggsave("~/GitHubs/myelo/docs/ktrSens.png",width=5, height=4)

idata=data.frame(ID=1:3,gam=0.161*span)
lb=labs(col = "gam")
scd = scale_color_discrete(name = "gam", labels = idata$gam)
out=mod%>%idata_set(idata)%>%ev(e)%>%mrgsim(start=0,end=END,delta=DELTA)
out%>%fixout%>%ggplot(aes(x=time,y=Circ,col=ID))+ geom_line(size=1)+gx+gy+tc(14)+lb+scd
# ggsave("~/GitHubs/myelo/docs/gamSens.png",width=5, height=4)

idata=data.frame(ID=1:3,slope=8.58*span)
lb=labs(col = "slope")
scd = scale_color_discrete(name = "slope", labels = idata$slope)
out=mod%>%idata_set(idata)%>%ev(e)%>%mrgsim(start=0,end=END,delta=DELTA)
out%>%fixout%>%ggplot(aes(x=time,y=Circ,col=ID))+ geom_line(size=1)+gx+gy+tc(14)+lb+scd
# ggsave("~/GitHubs/myelo/docs/slopeSens.png",width=5, height=4)


######### steady state offset of control?
span=c(1/1000000)
idata=data.frame(ID=1:length(span),gam=0.161*span)
lb=labs(col = "gam")
scd = scale_color_discrete(name = "gam", labels = idata$gam)
D=mod%>%idata_set(idata)%>%ev(e)%>%mrgsim(start=0,end=100000000,delta=1000)%>%fixout
D%>%filter(time>99999999)
####### no!!! you can always go out to larger times to hit the setpoint Circ0 with smaller gam

