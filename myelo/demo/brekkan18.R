library(tidyverse)
library(myelo)
library(mrgsolve)
(ka=0.0114*24) # 0.2736 (down from old value in Quartino14)
# (kcort=log(2)/(25.5/24)) = 0.6523738  
(MTT=120/24)
(ktr=5/MTT)
(kcirc=log(2)/(7/24))
Circ0=2.7
(trans_0=kcirc*Circ0/ktr) # 6.416562  (bigger number => more reserves/inventory)
(R1=2.84*24)  #R1 in paper = 68.16 mg/day (rate of injection, ignore here by using bolus)
(CL=2.80e-3*24)  #0.0672 10^9 cells L/h  

#### check steady state without PG
# Goal 1: find Cpg natural level to get margination factor at 1 for natural steady state
EC50=9.24    # ng/mL
EC50scale=0.477
EmaxMat=102
EmaxProl=109
EmaxScale=0.0622

Cpg=.0994374 # trial and error to get 1 on marg factor below
(FBprol=EmaxProl*Cpg/(EC50+Cpg))   # 1.160528
(FBtr=EmaxMat*Cpg/(EC50+Cpg)  ) # 1.085998
(FBmarg=EmaxProl*EmaxScale*Cpg/(EC50*EmaxScale+Cpg)) # 1

## Goal 2, find kin to counter background GCSF degradation rates 
V1=1.81 
(GCSF0=V1*Cpg)
Circ0=2.7
(CL/V1)*Circ0  # 0.1002431
Vmax=1.1208 
Km=2.05   # ng/mL of G-CSF
Vmax*Cpg/(Km+Cpg)  #0.0518505
(kin=(CL/V1)*Circ0+Vmax*Cpg/(Km+Cpg))  # 0.1520936

# Goal 3, find initial values of internal cell number states, starting with last, and then prol
(Trans4_0=kcirc*Circ0/FBtr) # 5.908445 all trans are the same
(Prol=5.908445*FBtr/FBprol)  # 5.529004


code='
$PROB A Population Pharmacokinetic-Pharmacodynamic Model of Pegfilgrastim
Ari Brekkan, Luis Lopez-Lazaro, Gunnar Yngman, Elodie L. Plan, Chayan Acharya, 
Andrew C. Hooker, Suresh Kankanwadi, and Mats O. Karlsson
The AAPS Journal (2018) 20: 91
healthy subjects (n = 192) were treated with single 6 mg doses of PEG-GCSF
Like Quartino 2014, but 
1) no cell killing, so Prol state => constant pressure boundary condition.
2) remove power feedback laws on Prol state regrowth


$PARAM Circ0=2.7,ktr=1,kcirc=2.376505,Prol=5.529004, ka=0.2736, CL=0.0672
V1=1.8 // Liters
Vmax=1.1208 // 0.0467*24 mg/day out via neut uptake
Km=2.05   // ng/mL of G-CSF
kin=0.1520936  // mg of GCSF/day
EC50=9.24    // ng/mL
EC50scale=0.477
EmaxMat=102
EmaxProl=109
EmaxScale=0.0622

$INIT Trans1=5.908445,Trans2=5.908445,Trans3=5.908445,Trans4=5.908445
      Circ=2.7
      GCSF=0.1799817,GCSFdepot=0
      
$ODE 
double Cpg=GCSF/V1;
double FBprol=EmaxProl*Cpg/(EC50+Cpg);
double FBtr=EmaxMat*Cpg/(EC50+Cpg);
double FBmarg=EmaxProl*EmaxScale*Cpg/(EC50*EmaxScale+Cpg);
//V1=V1/FBmarg // V1 defines Cpg, which defines FBmarg, which defines V1 again?? ... not clear how to do this!!! 
dxdt_Trans1=ktr*Prol*FBprol-ktr*Trans1*FBtr;
dxdt_Trans2=ktr*Trans1*FBtr-ktr*Trans2*FBtr;
dxdt_Trans3=ktr*Trans2*FBtr-ktr*Trans3*FBtr;
dxdt_Trans4=ktr*Trans3*FBtr-ktr*Trans4*FBtr;
dxdt_Circ=ktr*Trans4*FBtr-kcirc*Circ;
dxdt_GCSF=kin + ka*GCSFdepot - ((CL/V1)*Circ+Vmax*Cpg/(Km+Cpg)) ;
dxdt_GCSFdepot=-ka*GCSFdepot;
$CAPTURE FBprol, FBtr, FBmarg, Cpg 
'
mod <- mread("brekkan18", "~/tmp", code)
# (e=ev(time=5,amt=0,cmt=7)) 
(e=ev(time=5,amt=6,cmt=7)) 
out=mod%>%ev(e)%>%mrgsim(start=0,end = 30, delta = .1)
out%>%plot(xlab="Days")
tail(out)
d=as.data.frame(out)
head(d)
D=d%>%select(time,Trans1:Circ)%>%gather(key="Cell",value="Value",-time)%>%mutate(Cell=as_factor(Cell))
gx=xlab("Days")
sbb=theme(strip.background=element_blank())
tc=function(sz) theme_classic(base_size=sz)
D%>%ggplot(aes(x=time,y=Value,col=Cell))+geom_line(size=1)+gx+tc(14)+sbb
# should peak at 30 at 50 hours = 2 days out = Day 7. Shy because I didn't implement demargination 


