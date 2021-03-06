#  Characterization of Endogenous G-CSF and the Inverse Correlation to Chemotherapy-Induced Neutropenia in Patients with Breast Cancer Using Population Modeling 


[Quartino et al  *Pharm Res* **31**  3390-3403 (2014)](https://www.ncbi.nlm.nih.gov/pubmed/24919931) 
provide a model of 6 neutrophil lineage cell state variables, one proliferating, 
4 transitioning through maturation stages, and one circulating in  blood. 
Their model also includes GCSF formation, including due to corticosteroids given to reduce emesis due to docetaxel, and removal both by the kidney (ke) and internalization by
circulating neutrophils (Circ). Their model has separate control parameters for the 
effects of GCSF on maturation rates (beta) and proliferation rates (gamma).

![](../../docs/quartino.png)


Using Metrum Research Group's mrgsolve, their model is
```
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
```

Some model parameter values, given in the paper in hours, were 
converted to rates in days as follows.    
```
(MTT=133/24)
(ktr=5/MTT)
(kcirc=log(2)/(7/24))
(ke=.592*24)
(kanc=5.64*24)
GCSF0=24.3 #ng/ml
(kin=(ke+kanc*3.53)*GCSF0)  # 11956.3
(kcort=log(2)/(25.5/24))
Circ0=3.53
(Trans4=kcirc*Circ0/ktr)
```
These values were manually pasted into the model definition string above; the last line 
generates initial (and steady state) numbers of precursor cells (Circ0 is 
the initial number of circulating neutrophils).

The folllowing  runs the model and plots results using mrgsolve's plot method.
```
library(tidyverse)
library(mrgsolve)
mod <- mread("quart14", "~/tmp", code)
(e=ev(time=10,amt=80*1.8,cmt=1)) 
out=mod%>%ev(e)%>%mrgsim(start=0,end = 40, delta = 1)
out%>%plot(xlab="Days")
```
![](../../docs/quart14mrg.png)
 

To use ggplot2 results are first converted to a long dataframe as follows
```
d=as.data.frame(out)
D=d%>%select(time,Prol:Circ)%>%gather(key="Cell",value="Value",-time)%>%mutate(Cell=as_factor(Cell))
gx=xlab("Days")
sbb=theme(strip.background=element_blank())
tc=function(sz) theme_classic(base_size=sz)
D%>%ggplot(aes(x=time,y=Value,col=Cell))+geom_line(size=1)+gx+tc(14)+sbb
```
![](../../docs/quart14waves.png)


## Simulated Data Fitting
We now see if we can recover parameters from simulated data
```
(e=ev(time=0,amt=80*1.8,cmt=1)) 
END=50
DELTA=1
out=mod%>%ev(e)%>%mrgsim(start=0,end=END,delta=DELTA)
d=as.data.frame(out)[-1,]
sd=0.05
d$ANC=d$Circ+rnorm(dim(d)[1],sd=sd)
d%>%ggplot(aes(x=time,y=Circ))+ geom_line(size=.1)+
  geom_point(aes(x=time,y=ANC),size=1)+gx+tc(14)+sbb
```
![](../docs/quart14noiseData.png)
