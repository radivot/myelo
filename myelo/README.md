# myelo
This R package houses models of myeloid hematopoiesis. To install it use:  
```
devtools::install_github("radivot/myelo",subdir="myelo")
```

<!--
# Pharmacokinetic and -dynamic modelling of G-CSF derivatives in humans

The model of   [Scholz et al  *Theoretical Biology and Medical Modelling* **9** 32 (2012)](https://www.ncbi.nlm.nih.gov/pubmed/22846180) 
-->


# Model of Chemotherapy-Induced Myelosuppression With Parameter Consistency Across Drugs

[Friberg et al  *J Clin Oncol* **20**  4713-4721 (2002)](https://www.ncbi.nlm.nih.gov/pubmed/12488418) 
provide a model of 5 neutrophil lineage cell state variables, one proliferating, 
3 transitioning through maturation stages, and one circulating in  blood. 
![](../docs/friberg02graph.png)


In R their model is 
```
friberg02<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dC1=-(k12+k13+k10)*C1 + k21*C2 + k31*C3
    dC2=k12*C1 - k21*C2
    dC3=k13*C1 - k31*C3
    Cp=C1/V1/mw  #mw in mg/umole => uM
    Edrug=slope*Cp
    dProl = ktr*Prol*(1-Edrug)*(Circ0/Circ)^gam - ktr*Prol
    dTrans1=ktr*Prol-ktr*Trans1
    dTrans2=ktr*Trans1-ktr*Trans2
    dTrans3=ktr*Trans2-ktr*Trans3
    dCirc=ktr*Trans3-ktr*Circ
    return(list(c(dC1,dC2,dC3,dProl,dTrans1,dTrans2,dTrans3,dCirc)))
  })
}

```
Using fits of this model to docetaxel response data, the parameter estimates, with time in days, are
```
fribergPars02=c(Circ0=5.05, ktr=24*4/88.7,gam=0.161,slope=8.58, 
k12 = 24*1.06,k21 = 24*1.51,k13 = 24*1.26, k31 = 24*0.084,k10 = 24*5.2,V1=7.4,mw=0.808)
```
Starting from a steady state at t=-5 days, the code below simulates adding
 a bolus  of  docetaxel (100 mg/m2 *1.8 m2 = 180 mg) at t=0. 
```
library(tidyverse)
library(deSolve)
library(myelo)
times <- c(-5:0,seq(0,1,0.01),1:25)
x0=c(C1=0,C2=0,C3=0,Prol=5.05,Trans1=5.05,Trans2=5.05,Trans3=5.05,Circ=5.05)
(evnt=data.frame(var="C1",time=0,value=180,method="rep"))
yout=ode(x0,times=times,func=friberg02,events=list(data=evnt),parms=fribergPars02)
D=as.data.frame(yout)
gx=xlab("Days")
sbb=theme(strip.background=element_blank())
tc=function(sz) theme_classic(base_size=sz)
d=D%>%select(time,C1:Circ)%>%gather(key="Lab",value="Value",-time)
d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb
ggsave("~/Results/myelo/fri02.png",width=5, height=6)
```
![](../docs/fri02.png)
Note that the initial number of circulating neutrophils (Circ),  5.05e3 cells/uL, is 
also the setpoint, Circ0, and from the equations, also the steady state of each  compartment. 

The same model in faster  deSolve C code, is
```
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
static double  parms[11];
#define Circ0  parms[0]
#define ktr    parms[1]
#define gam    parms[2]
#define slope  parms[3]
#define k12    parms[4]
#define k21    parms[5]
#define k13    parms[6]
#define k31    parms[7]
#define k10    parms[8]
#define V1     parms[9]
#define mw     parms[10]


void parmsFri02(void (* odeparms)(int *, double *))
{   int N=11;
    odeparms(&N, parms);
}


void derivsFri02(int *neq, double *t, double *y, double *ydot)
{
    double C1, C2,C3,Prol,Trans1,Trans2,Trans3,Circ, Edrug;
    double dC1, dC2,dC3,dProl,dTrans1,dTrans2,dTrans3,dCirc;
    C1=y[0];  C2=y[1];   C3=y[2];  Prol=y[3];   
    Trans1=y[4];  Trans2=y[5];   Trans3=y[6];  Circ=y[7];   
    
    dC1=-(k12+k13+k10)*C1 + k21*C2 + k31*C3;
    dC2=k12*C1 - k21*C2;
    dC3=k13*C1 - k31*C3;
    Edrug=slope*C1;
    dProl = ktr*Prol*(1-Edrug)*pow(Circ0/Circ,gam) - ktr*Prol;
    dTrans1=ktr*Prol-ktr*Trans1;
    dTrans2=ktr*Trans1-ktr*Trans2;
    dTrans3=ktr*Trans2-ktr*Trans3;
    dCirc=ktr*Trans3-ktr*Circ;
    
    ydot[0] = dC1 ; 
    ydot[1] = dC2;
    ydot[2] = dC3;
    ydot[3] = dProl;
    ydot[4] = dTrans1;
    ydot[5] = dTrans2;
    ydot[6] = dTrans3;
    ydot[7] = dCirc;
}

```

Running this C version of the model is done as follows
```
(f=file.path(system.file(paste("libs",Sys.getenv("R_ARCH"),sep=""), package = "myelo"),
paste("myelo",.Platform$dynlib.ext,sep="")))
dyn.load(f)
yout=ode(x0,times=times,func="derivsFri02",
       dllname = "myelo",initfunc = "parmsFri02",
       events=list(data=evnt),parms=fribergPars02)
plot(yout) 
```
![](../docs/fri02deSolve.png)


Using Metrum Research Group's mrgsolve, such C code is automatically generated and compiled using this neat R code.
```
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
```
which generates

![](../docs/mrgFri02.png)
![](../docs/waves.png)


mrgsolve also offers nice Rx short hand as follows: 
```
(e=ev_rx("50 q 21 x 6 then 25 q 14 x 6"))
mod%>%ev(e)%>%mrgsim(end = 300, delta = 0.1)%>%plot(xlab="Days")
```

![](../docs/Fribolus.png)

and the ability to easily change the dose to a continuous infusion, here over 12 of the 21 days
of a cycle
```
(e=ev_rx("50 over 12 q 21 x 6 then 25 q 14 x 6"))
mod%>%ev(e)%>%mrgsim(end = 300, delta = 0.1)%>%plot(xlab="Days")
```
![](../docs/infus12.png)

To see if we can recover 4 parameters from a simulation, we first simulate some data

```
(e=ev(time=0,amt=180,cmt=1)) 
out=mod%>%ev(e)%>%mrgsim(start=0,end = 50, delta = .1)
d=as.data.frame(out)
sd=0.1
d$ANC=d$Circ+rnorm(dim(d)[1],sd=sd)
d%>%ggplot(aes(x=time,y=Circ))+ geom_line(size=.1)+
  geom_point(aes(x=time,y=ANC),size=1)+gx+tc(14)+sbb
ggsave("~/GitHubs/myelo/docs/noiseData.png",width=5, height=4)
```
![](../docs/noiseData.png)


and then see if FME can retrieve the estimates 

```
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
```

The output of this below shows that while the 1st three parameters moved in the
correct direction, the fourth did not. Strong negative correlation between gamma and slope 
hints at the problem: stronger chemo can look weaker if ANC control is stronger.
```
> pars
 Circ0    ktr    gam  slope 
5.0500 1.0823 0.1610 8.5800 
> coef(Fit)
     Circ0        ktr        gam      slope 
 4.4233503  1.4563017  0.1647432 27.0110295 
> 1.5*pars
   Circ0      ktr      gam    slope 
 7.57500  1.62345  0.24150 12.87000 
> summary(Fit)

Parameters:
       Estimate Std. Error t value Pr(>|t|)    
Circ0  4.423350   0.138695   31.89   <2e-16 ***
ktr    1.456302   0.033024   44.10   <2e-16 ***
gam    0.164743   0.006138   26.84   <2e-16 ***
slope 27.011029   0.832308   32.45   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.572 on 248 degrees of freedom

Parameter correlation:
        Circ0     ktr     gam   slope
Circ0  1.0000 -0.6384  0.3664 -0.5419
ktr   -0.6384  1.0000  0.4584 -0.2954
gam    0.3664  0.4584  1.0000 -0.9570
slope -0.5419 -0.2954 -0.9570  1.0000
```




