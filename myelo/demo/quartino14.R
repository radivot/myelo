library(tidyverse)
# library(myelo)
library(mrgsolve)
(MTT=133/24)
(ktr=5/MTT)
(kcirc=log(2)/(7/24))
(ke=.592*24)
(kanc=5.64*24)
GCSF0=24.3 #ng/ml
(kin=(ke+kanc*3.53)*GCSF0)  # 11956.3

(kcort=log(2)/(25.5/24))
 
Circ0=3.53
trans4_0=kcirc*Circ0/ktr

# $INIT C1=0,C2=0,C3=0,Prol=9.297876,Trans1=9.297876,Trans2=9.297876
# Trans3=9.297876,Trans4=9.297876,Circ=3.53
# //double FBprol=GCSF>0.001 ? pow(GCSF/GCSF0,beta) : pow(0.001/GCSF0,beta);
# //double FBtr=GCSF>0.001 ? pow(GCSF/GCSF0,gam) : pow(0.001/GCSF0,gam);


code='
$PROB Characterization of Endogenous G-CSF and the Inverse Correlation to 
Chemotherapy-Induced Neutropenia in Patients with Breast Cancer Using Population Modeling 
Angelica L. Quartino & Mats O. Karlsson & Henrik Lindman & Lena E. Friberg
Pharm Res (2014) 31:3390-3403
Only doing Docetaxel as in Fri02, so keeping PK from there.

$PARAM Circ0=3.53,ktr=0.9022556,kcirc=2.376505,gam=0.444,beta=0.234,slope=17.2
ke=14.208, kanc=135.36, GCSF0=24.3, kcort=0.6523738, kin=11956.3
k12=25.44,k21=36.24,k13=30.24,k31=2.016,k10=124.8,V1=7.4,mw=0.808

$INIT C1=0,C2=0,C3=0,Prol=9.562973,Trans1=9.562973,Trans2=9.562973
      Trans3=9.562973,Trans4=9.562973,Circ=3.630645
      GCSF=24.3,GCSFcort=0
      
$ODE 
double Cp=C1/V1/mw;
double Edrug=slope*Cp;
double FBprol=pow(GCSF/GCSF0,gam);
double FBtr=pow(GCSF/GCSF0,beta);
dxdt_C1=-(k12+k13+k10)*C1 + k21*C2 + k31*C3;
dxdt_C2=k12*C1 - k21*C2;
dxdt_C3=k13*C1 - k31*C3;
dxdt_Prol = ktr*Prol*(1-Edrug)*FBprol - ktr*Prol*FBtr;
dxdt_Trans1=ktr*Prol*FBtr-ktr*Trans1*FBtr;
dxdt_Trans2=ktr*Trans1*FBtr-ktr*Trans2*FBtr;
dxdt_Trans3=ktr*Trans2*FBtr-ktr*Trans3*FBtr;
dxdt_Trans4=ktr*Trans3*FBtr-ktr*Trans4*FBtr;
dxdt_Circ=ktr*Trans4*FBtr-kcirc*Circ;
dxdt_GCSF=kin-(ke+kanc*Circ*GCSF) + kcort*GCSFcort;
dxdt_GCSFcort=-kcort*GCSFcort;
$CAPTURE FBprol, FBtr 
'
mod <- mread("quart14", "~/tmp", code)
tstart=10
(e=ev(time=tstart,amt=80*1.8,cmt=1)) 
out=mod%>%ev(e)%>%mrgsim(start=0,end = 50, delta = 1)
out%>%plot(xlab="Days")
tail(out)
d=as.data.frame(out)
D=d%>%select(time,Prol:Circ)%>%gather(key="Cell",value="Value",-time)%>%mutate(Cell=as_factor(Cell))
gx=xlab("Days")
sbb=theme(strip.background=element_blank())
tc=function(sz) theme_classic(base_size=sz)
D%>%ggplot(aes(x=time,y=Value,col=Cell))+geom_line(size=1)+gx+tc(14)+sbb
# ggsave("~/GitHubs/myelo/docs/quart14waves.png",width=5, height=4)

# simulated data
tstart=10
(e=ev(time=tstart,amt=80*1.8,cmt=1)) 
END=50
DELTA=1
out=mod%>%ev(e)%>%mrgsim(start=0,end=END,delta=DELTA)
d=as.data.frame(out)%>%filter(!(time==tstart&C1==0))
sd=0.05
d$ANC=d$Circ+rnorm(dim(d)[1],sd=sd)
d%>%ggplot(aes(x=time,y=Circ))+ geom_line(size=.1)+
  geom_point(aes(x=time,y=ANC),size=1)+gx+tc(14)+sbb
ggsave("~/GitHubs/myelo/docs/quart14noiseData.png",width=5, height=4)

# Try fitting 4 params that could be fitted using Friberg 2002 model
library(FME)
# (pars=c(Circ0 = 5.05, ktr = 1.08229988726043,gam = 0.161, slope = 8.58)) # fri02 values
# (pars=c(    Circ0=3.53, ktr=0.9022556,         gam=0.444,   slope=17.2))  
# Q=function(pars) {
#   evnt=ev(time=tstart,amt=80*1.8,cmt=1,Circ0=pars["Circ0"],ktr=pars["ktr"],
#           gam=pars["gam"],
#           slope=pars["slope"],
#           Prol_0=pars["Circ0"],Trans1_0=pars["Circ0"],Trans2_0=pars["Circ0"],
#           Trans3_0=pars["Circ0"],Trans4_0=pars["Circ0"],Circ_0=pars["Circ0"] )
#   as.data.frame(mod%>%ev(evnt)%>%mrgsim(start=0,end = END, delta = DELTA))%>%filter(!(time==tstart&C1==0))
# }  # this Q would not move toward true values

# (pars=c(    Circ0=3.53, ktr=0.9022556,         gam=0.444,   beta=0.234))  
# Q=function(pars) {
#   evnt=ev(time=tstart,amt=80*1.8,cmt=1,Circ0=pars["Circ0"],ktr=pars["ktr"],
#           gam=pars["gam"],
#           beta=pars["beta"],
#           Prol_0=pars["Circ0"],Trans1_0=pars["Circ0"],Trans2_0=pars["Circ0"],
#           Trans3_0=pars["Circ0"],Trans4_0=pars["Circ0"],Circ_0=pars["Circ0"] )
#   as.data.frame(mod%>%ev(evnt)%>%mrgsim(start=0,end = END, delta = DELTA))%>%filter(!(time==tstart&C1==0))
# } #better but still not good


#drop to two parameters
(pars=c(    Circ0=3.53, ktr=0.9022556))  
Q=function(pars) {
  evnt=ev(time=tstart,amt=80*1.8,cmt=1,Circ0=pars["Circ0"],ktr=pars["ktr"],
          Prol_0=pars["Circ0"],Trans1_0=pars["Circ0"],Trans2_0=pars["Circ0"],
          Trans3_0=pars["Circ0"],Trans4_0=pars["Circ0"],Circ_0=pars["Circ0"] )
  as.data.frame(mod%>%ev(evnt)%>%mrgsim(start=0,end = END, delta = DELTA))%>%filter(!(time==tstart&C1==0))
} 


Q(pars)%>%select(time,Circ)%>%head(12)  # test that Q works
(dd=d%>%select(time,Circ=ANC)%>%mutate(sd=sd))%>%head(12)
Qcost <- function (pars) {
  out=Q(pars)%>%select(time,Circ)
  modCost(model = out, obs = dd, err = "sd")
}
parsIC=1.3*pars
Qcost2 <- function(lpars)  Qcost(exp(lpars))
Fit <- modFit(f = Qcost2, p = log(parsIC),method="Nelder-Mead")
Fit <- modFit(f = Qcost2, p = log(parsIC)) 
data.frame(parsIC,fit=exp(coef(Fit)),trueVals=pars)

(s=summary(Fit))
data.frame(point=exp(s$par[,1]),
           lowCI=exp(s$par[,1]-1.96*s$par[,2]),
           hiCI=exp(s$par[,1]+1.96*s$par[,2])  )


##### Perturbation Study
DELTA=.1
(e=ev(time=tstart,amt=80*1.8,cmt=1)) 
Circ0=3.53 
ktr=0.9022556
(kcirc=log(2)/(7/24))
(trans0=kcirc*Circ0/ktr)

ic=c(Circ0=Circ0,Prol_0=trans0,Trans1_0=trans0,Trans2_0=trans0,Trans3_0=trans0,Trans4_0=trans0,Circ_0=Circ0)
span=c(1/2, 1, 2)
idata=data.frame(ID=1:3,span%*%t(ic))
out=mod%>%idata_set(idata)%>%ev(e)%>%mrgsim(start=tstart,end=END,delta=DELTA)
fixout=function(out) {
  d$ID=as_factor(d$ID)
  d%>%filter(!(time==tstart&C1==0))
}
lb=labs(col = "Circ0")
scd = scale_color_discrete(name = "Circ0", labels = idata$Circ0)
out%>%fixout%>%ggplot(aes(x=time,y=Circ,col=ID))+ geom_line(size=1)+gx+gy+tc(14)+lb+scd
# ggsave("~/GitHubs/myelo/docs/circ0Sens.png",width=5, height=4)

idata=data.frame(ID=1:3,ktr=0.9022556*span)
lb=labs(col = "ktr")
scd = scale_color_discrete(name = "ktr", labels = idata$ktr)
out=mod%>%idata_set(idata)%>%ev(e)%>%mrgsim(start=0,end=END,delta=DELTA)
out%>%fixout%>%ggplot(aes(x=time,y=Circ,col=ID))+ geom_line(size=1)+gx+gy+tc(14)+lb+scd
# ggsave("~/GitHubs/myelo/docs/ktrSens.png",width=5, height=4)

idata=data.frame(ID=1:3,gam=0.444*span)   
lb=labs(col = "gam")
scd = scale_color_discrete(name = "gam", labels = idata$gam)
out=mod%>%idata_set(idata)%>%ev(e)%>%mrgsim(start=0,end=END,delta=DELTA)
out%>%fixout%>%ggplot(aes(x=time,y=Circ,col=ID))+ geom_line(size=1)+gx+gy+tc(14)+lb+scd
# ggsave("~/GitHubs/myelo/docs/gamSens.png",width=5, height=4)

idata=data.frame(ID=1:3,slope=17.2*span)
lb=labs(col = "slope")
scd = scale_color_discrete(name = "slope", labels = idata$slope)
out=mod%>%idata_set(idata)%>%ev(e)%>%mrgsim(start=0,end=END,delta=DELTA)
out%>%fixout%>%ggplot(aes(x=time,y=Circ,col=ID))+ geom_line(size=1)+gx+gy+tc(14)+lb+scd
# ggsave("~/GitHubs/myelo/docs/slopeSens.png",width=5, height=4)



