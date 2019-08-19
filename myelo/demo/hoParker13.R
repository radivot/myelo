library(tidyverse)
library(myelo)
library(mrgsolve)
code='
$PROB A model of neutrophil dynamics in response to inflammatory and cancer chemotherapy challenges
Thang Ho, Gilles Clermont, and Robert S. Parker
Computers and Chemical Engineering 51 (2013) 187-196

$PARAM 
k1=0.0006
k2=14.5511
k3=2.36e-8
k4=8000
k5=0.0016
k6=0.25
k7=105
k8=107
k9=0.0858
k10=0.11843
k11=0.1184
k12=0.0253
k13=0.3603
k14=106
k15=8.48e-7
k16=24.9738
k17=50000
k18=0.0062
k19=0.0922
k20=0.002
k21=5e4
k22=0.0188
k23=3.2e4
k24=0.0085
k25=2.5e3
k26=0.001
k27=0.0066
k28=107
k29=0.004
k30=9e-11
k31=0.00073
k32=0.00159
k33=6e4
k34=8.4170
k35=104
k36=3.5e9
k37=77.5
k38=0.0015
k39=0.00045
k40=0.0091
k41=0.000035
Circ0=5.05,gam=0.161,slope=8.58
k12d=25.44,k21d=36.24,k13d=30.24,k31d=2.016,k10d=124.8,V1=7.4,mw=0.808

$INIT 
C1=0,C2=0,C3=0
S=5.05,Prol=5.05,Trans1=5.05,Trans2=5.05,Trans3=5.05,Nc=5.05,Np=5.05,Nd=5.05
I1=0,IL1=0,IL17=0,IL23=0,LPS=0,GCSF=24.3,GCSFT=0 
Ti=0,Ta=0

$ODE 
double Cp=C1/V1/mw;
double D=slope*Cp;
double ktr=k1*(1 + k2*(GCSF/(k4 + GCSF) )*k3/(k3+D) );
dxdt_C1=-(k12d+k13d+k10d)*C1 + k21d*C2 + k31d*C3;
dxdt_C2=k12d*C1 - k21d*C2;
dxdt_C3=k13d*C1 - k31d*C3;

dxdt_LPS=- k39*LPS;
dxdt_I1= k39*LPS + k40*I1;
dxdt_Prol = (k41 + k32*GCSF/(k33 + GCSF))*S - ktr*Prol - k37*Prol*D/(k38+D);
dxdt_Trans1=ktr*Prol-ktr*(1+1.5*I1/(k36+I1) + 1.5*IL1/(k36+IL1) )*Trans1;
dxdt_Trans2=ktr*Trans1-ktr*(1+1.5*I1/(k36+I1) + 1.5*IL1/(k36+IL1) )*Trans2;
dxdt_Trans3=ktr*Trans2-ktr*(1+1.5*I1/(k36+I1) + 1.5*IL1/(k36+IL1) )*Trans3;
dxdt_Nc=ktr*Trans3 - GCSFT - 1.5*ktr*(I1/(k36+I1) + 1.5*IL1/(k36+IL1))*(Trans1+Trans2+Trans3) +
       -ktr*Nc*(1-IL1/(k36+IL1)) + ktr*(1  +  10*IL1/(k36+IL1)  +  k34*GCSF/(k35+GCSF)  )*Np;
dxdt_Np=-ktr*(1  +  10*IL1/(k36+IL1)  +  k34*GCSF/(k35+GCSF)  )*Np + ktr*Nc*(1-IL1/(k36+IL1)) ;
dxdt_Nd=k5*Nc/(k8+Nc)*(1-k6*GCSF/(k7+GCSF))*Nc - k9*Nd;
dxdt_IL23=k10-k11*Nd/(k7+Nd) - k12*IL23;
dxdt_Ta= (k15 +  pow(IL23,2)/(pow(k16,2)+pow(IL23,2)))*Ti - Ta*k13*pow(k14,2)/(pow(k14,2)+pow(Ta,2));
dxdt_Ti= -(k15 +  pow(IL23,2)/(pow(k16,2)+pow(IL23,2)))*Ti + Ta*k13*pow(k14,2)/(pow(k14,2)+pow(Ta,2));
dxdt_IL17= k18*Ta - k19*IL17;
dxdt_S= k20*IL17*IL17/(k21+IL17) - k22*k23*S/(k23+S);
dxdt_GCSF= k24*S*pow(S,2)/(pow(k25,2)+pow(S,2))*Ti - (k26+k27*Nc/(k28+Nc))*GCSF +k29*GCSFT-k30*GCSF;
dxdt_GCSFT= -k29*GCSFT+k30*GCSF-k31*GCSFT;
'   # had to guess  tissue in dxdt_Nc was  GCSFT, but it doesn't make sense why this would be in linearly
mod <- mread("hoParker13", "~/tmp", code)
(e=ev(time=0,amt=180,cmt=1)) 
out=mod%>%ev(e)%>%mrgsim(start=0,end = 25, delta = 1)
out%>%plot(xlab="Days")
### Conclusion: model looks way overparameterized. No initial conditions => not really meant to be reproduced
