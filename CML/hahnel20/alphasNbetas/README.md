## Alpha and Betas in Hahnel et al. Cancer Research (2020).

Goal: Use BCR-ABL percentages above detection and before cessation to demo 
population level modeling of bi-exponential BCR-ABL 
decays often observed initially with TKI use. 

The model of interest is y=Ae<sup>αt</sup>+Be<sup>βt</sup>. When the first term dominates, the slope on a natural 
log scale is α and when the second term dominates the slope is β. On a log10 scale, these
values are log(10)=2.3 times smaller, as is the residual error of the fits. We will plot 
data and fits on a log10 scale, but we will perform fits on a natural log scale. Taking logs gives small y values more weight in estimates of β. 


First, as in Hahnel et al, maxLik is used to fit the model.
```
rm(list=ls()) 
library(tidyverse)
library(myelo)
head(d<-hahnelFigS2)
(d=d%>%filter(UL==0)%>%select(-UL)) # use only non-zero measurements to keep it simple
names(d)=c("TIME","DV","ID")  #change names to NONMEM style
d=d%>%mutate(LDV=log10(DV),LNDV=log(DV)) 
library(maxLik)
loglik <- function(param) {
  lA <- log(param[1])
  alpha <- param[2]
  lB <- log(param[3])
  beta <- param[4]
  mu <- log(exp(lA+alpha*t)+exp(lB+beta*t))
  sigma <- param[5]
 sum(dnorm(x, mu, sigma, log=TRUE))
} #indirect passing of x and t via global is a weakness of maxLik
(d21=d%>%group_by(ID)%>%summarize(Tf=max(TIME)))
d21$Tpred=lapply(d21$Tf,function(x) seq(0,x))
L=NULL
for (i in 1:21) {   # this weakness drove my use here of a for loop
  print(i)
  d1=d%>%filter(ID==i)
  print(d1)
  t=d1$TIME
  x=d1$LNDV
  N=length(x)
  L[[i]]=summary(maxLik(loglik, start=c(A=100,alpha=-1,B=0.1,beta=-0.05,sigma=0.5)))  
}
(L=lapply(L,coef))
(P=lapply(L,function(x) x[,1]))
d21$P=P
d21
simY <- function(param,tpred,id) {
  A <- param[1]
  alpha <- param[2]
  B <- param[3]
  beta <- param[4]
  data.frame(t=tpred,y=log10(A*exp(alpha*tpred)+B*exp(beta*tpred)),ID=id)
}
(D=mapply(simY,P,d21$Tpred,1:21,SIMPLIFY = F))
library(data.table)
(D=data.table(bind_rows(D)))
names(D)=c("TIME","LDV","ID")
mkDF=function(x) data.frame(TIME=100,LDV=2.5,a=round(x[["alpha"]],3),b=round(x[["beta"]],3))
(DTx=bind_rows(lapply(P,mkDF)))
DTx$ID=1:21
DTxb=DTx
DTxb$LDV=1
d%>%ggplot(aes(x=TIME,y=LDV))+facet_wrap(ID~.,ncol=4)+geom_point(size=1)+
   geom_line(data=D)+geom_text(aes(label=a),data=DTx,parse=T)+geom_text(aes(label=b),data=DTxb,parse=T)
setwd("~/githubs/myelo/CML")
ggsave("../docs/alphaNbetaDataNfits.png",width=7,height=8)
```

![](../../../docs/alphaNbetaDataNfits.png)
which shows α and β values ~2.3-fold higher than log10 scale α and β slopes in Fig S2. 
Fits can be improved by tuning initial estimates.

Using the same initial estimates and the same 
optimization method (Nelder-Mead), but using bbmle instead of maxLik, the code 
is cleaner, as variables are no longer passed to functions via globals.  
```
library(bbmle)
dn=d%>%group_by(ID)%>%nest()
fitBi=function(d1)  {
  coef(summary(mle2(LNDV~dnorm(mean=log((exp(lA+alpha*TIME) + exp(lB+beta*TIME))),sd=sigma),
               method="Nelder-Mead", 
               start=list(lA=log(100),alpha=-1,lB=log(0.1),beta=-0.05,sigma=0.5),data=d1,
               control = list(maxit=500))))  
}
dn=dn%>%mutate(M=map(data,fitBi))
dn$M[[20]]
getPars=function(x){
  x=x[,1] 
  x[c(1,3)]=exp(x[c(1,3)])
  names(x)[c(1,3)]=c("A","B")
  x
}
dn=dn%>%mutate(P=map(M,getPars))

(D=mapply(simY,dn$P,d21$Tpred,1:21,SIMPLIFY = F))
(D=data.table(bind_rows(D)))
names(D)=c("TIME","LDV","ID")
(DTx=bind_rows(lapply(dn$P,mkDF)))
DTx$ID=1:21
DTxb=DTx
DTxb$LDV=1
d%>%ggplot(aes(x=TIME,y=LDV))+facet_wrap(ID~.,ncol=4)+geom_point(size=1)+
  geom_line(data=D)+geom_text(aes(label=a),data=DTx,parse=T)+geom_text(aes(label=b),data=DTxb,parse=T)
ggsave("../docs/alphaNbetaDataNfitsBB.png",width=7,height=8)

```

![](../../../docs/alphaNbetaDataNfitsBB.png)
  

To bias estimates toward patient population averages, we can use the function nlmer in 
the mixed effects R package lme4. 

```
library(lme4)
(startvec <- c(lA = log(100), alpha = -1, lB=log(0.1), beta = -0.05))
nform <- ~ log(exp(lA+alpha*input) + exp(lB+beta*input))
nfun <- deriv(nform, namevec=c("lA", "alpha", "lB", "beta"),
              function.arg=c("input","lA", "alpha", "lB", "beta"))
nfun
(M<-nlmer(LNDV ~ nfun(TIME, lA, alpha, lB, beta) ~ (lA|ID) + (alpha|ID) +(lB|ID) + (beta|ID),
                        data=d,start = startvec,
      nAGQ = 0L,
      control = nlmerControl(tolPwrss = 1e-4)) )
(fe=fixef(M))
(RE=ranef(M)$ID)
apply(RE,2,mean) #check
(P=RE+t(t(rep(1,21)))%*%t(fe))
P=P%>%mutate(lA=exp(lA),lB=exp(lB))%>%rename(A=lA,B=lB)
P=lapply(1:21,function(i) unlist(P[i,]))
(D=mapply(simY,P,d21$Tpred,1:21,SIMPLIFY = F))
(D=data.table(bind_rows(D)))
names(D)=c("TIME","LDV","ID")
(DTx=bind_rows(lapply(P,mkDF)))
DTx$ID=1:21
DTxb=DTx
DTxb$LDV=1
d%>%ggplot(aes(x=TIME,y=LDV))+facet_wrap(ID~.,ncol=4)+geom_point(size=1)+
  geom_line(data=D)+geom_text(aes(label=a),data=DTx,parse=T)+geom_text(aes(label=b),data=DTxb,parse=T)
ggsave("../docs/alphaNbetaFitsLME4.png",width=7,height=8)

```

![](../../../docs/alphaNbetaFitsLME4.png)
We now have nice fits with parameter estimates that are more similar. 
Information was borrowed across patients in one joint fit of 
population parameter means and variances. 
Assuming normality, estimates for individuals were 
brought closer to the population means. 

A limitation of lme4 is that it requires closed form models. Models are often more simply described by 
ordinary differential equations (ODEs). An R package that allows this is nlmixr. In 
this pharmacokinetic data analysis package, states are masses 
(drug amounts),  rate constants are clearances divided by volumes, and measured
drug concentrations are central comparment drug masses divided by volumes.

Bi-exponential decays are captured by two-compartment models. We thus
view BCR-ABL as if injected into a central compartment from which it is cleared 
permanently or set aside temporarily in a peripheral compartment that protects it, until it returns to the
central compartment. Needed then is a mapping of our initial estimates above 
(A=100, alpha=-1, B=0.1, and beta=-0.05) to initial estmates of two-compartment PK model parameters.
To get the initial BCR-ABL value of 100.1, we set the dose to 1001 and 
our initial estimate of the central volume (Vc) to 10. We then set
the initial estimate of the elimination rate constant Ke to 1 (i.e. -alpha) 
and the initial estimate of the peripheral-to-central
rate constant Kpc to 0.05 (i.e. -beta). Setting the reverse,
Kcp, to a similarly small value of 0.05, the nlmixr code is

```
library(nlmixr)
dn=d # will add dose rows into dn 
dn$TIME[dn$TIME<=0.001] #need these to come after bolus injection for ODEs
dn$TIME[dn$TIME<=0.001]= 0.001 # some are slightly negative (near zero)
dn$EVID=0
dn$AMT=0
(dn=dn%>%group_by(ID)%>%nest())
dn$data[[1]]
# EVID = 10,000 × (1 If IV Infusion 0 If Bolus) +100 × (Compartment #)+1
(dtop=data.frame(TIME=0,DV=0,LDV=0,LNDV=0,EVID=101,AMT=1001)) #101 =>into compartment 1 
dn=dn%>%mutate(ndata=map(data,function(x) rbind(dtop,x)))
(dn=dn%>%select(-data)%>%unnest(cols = c(ndata)))

two.compartment.IV.model <- function(){
  ini({ # Where initial conditions/variables are specified
    lVc <- log(10)    #log Vc    
    lKe <- log(1)     #log Ke   
    lKpc<- log(0.05)  #log Kpc      
    lKcp<- log(0.05)  #log Kcp  
    prop.err <- 0.3 
    eta.Vc ~ 0.15   
    eta.Ke ~ 0.15
    eta.Kpc ~ 0.15  
    eta.Kcp  ~ 0.15
  })
  model({ # Where the model is specified
    Vc  <- exp(lVc + eta.Vc)
    Ke  <- exp(lKe + eta.Ke)
    Kpc <- exp(lKpc + eta.Kpc)
    Kcp <- exp(lKcp + eta.Kcp)
    # RxODE-style differential equations are supported
    d/dt(centr)  = Kpc*periph-Kcp*centr-Ke*centr;
    d/dt(periph) =-Kpc*periph+Kcp*centr;
    cp = centr / Vc;
    cp ~ prop(prop.err)
  })
}
(fit <- nlmixr(two.compartment.IV.model,dn,est="saem"))
pfit=nlmixrPred(fit,ipred=TRUE)
str(pfit)
dp=data.frame(ID=pfit$id,TIME=pfit$time,Y=pfit$ipred)%>%mutate(LDV=log10(Y))
d%>%ggplot(aes(x=TIME,y=LDV))+facet_wrap(ID~.,ncol=4)+geom_point(size=1)+geom_line(data=dp)
ggsave("../docs/nlmixrFit1.png",width=7,height=8)
```
![](../../../docs/nlmixrFit1.png)

Patient-specific model parameter estimates corresponding to the plots above can be retrieved as follows.

```
fe=coef(fit)$fixed[1:4]
re=coef(fit)$random[,2:5]
(P=re+t(t(rep(1,21)))%*%t(fe))
(P=exp(P))
names(P)=c("Vc","Ke","Kpc","Kcp")
P

           Vc        Ke         Kpc          Kcp
1  190.022803 0.9683176 0.023447195 1.024638e-02
2  410.873068 0.4121233 0.009171166 1.004519e-02
3  206.820753 0.8199347 0.089255357 2.691175e-02
4  281.549731 0.8425779 0.009688618 1.449773e-02
5  235.062980 0.7039342 0.026637963 2.904405e-03
6  165.676543 1.3410858 0.036949112 2.052292e-02
7  285.489040 0.9074194 0.025367178 3.498893e-02
8  764.512689 0.3948731 0.028553293 3.769384e-03
9  145.398173 1.0819822 0.071341269 4.845837e-03
10 170.026827 1.2494456 0.008808499 3.080718e-02
11 181.425513 0.7539363 0.015390202 4.596884e-03
12  72.967882 0.6335443 0.097858338 4.628165e-01
13   5.778486 0.5723832 0.005194983 9.861744e-05
14  28.600861 0.6537157 0.003690037 3.437080e-03
15   7.789124 0.2410669 0.006206287 3.068806e-05
16 206.344675 0.7796115 0.014283552 8.859635e-03
17 122.805873 0.8659598 0.023525108 1.587549e-02
18  24.659964 1.7006906 0.007634866 1.231971e-02
19  82.711098 0.4727110 0.023823067 3.658259e-03
20 134.782622 0.7838043 0.029108729 6.414310e-02
21 261.690853 0.8096706 0.014647168 6.283243e-03

```

Using the R package IQRtools (https://iqrtools.intiquan.com/) and IQdesktop
(https://iqdesktop.intiquan.com/) on DockerHub (https://hub.docker.com/r/intiquan/iqdesktop),
a single representation of the model can be applied
to three different backends, NONMEM, Monolix, and nlmixr. 
The code for this is as follows. 

```
library(IQRtools)
#setwd("~/PROJECTS/SHARE") # in IQ Desktop
#setwd("~/githubs/myelo/CML/hahnel20/alphasNbetas")  # in myelo folder
# data first needs to be put in a backend agnostic form
dn$EVID[dn$EVID==101]=1
dn$CMT=1
(d=dn%>%rename(USUBJID=ID,VALUE=DV)%>%mutate(TIMEUNIT="Months",
                                         VALUE=ifelse(EVID==1,1001,VALUE),
                                         ROUTE=ifelse(EVID==1,"IV",NA),
                                         UNIT=ifelse(EVID==1,"mg","prct"),
                                         NAME=ifelse(EVID==1,"Dose","BCRABL")))
(d=d%>%select(USUBJID,TIME,TIMEUNIT,VALUE,NAME,UNIT,ROUTE))
fQ="data/biExpQ.csv"
write.table(d,fQ,sep=",",row.names = FALSE)
(dq=IQRdataGENERAL(input=fQ)) 
exportNLME_IQRdataGENERAL(dq,filename = "data/biExport")
#now use on form of the data
data=data_IQRest(datafile= "data/biExport.csv")     
(dosing=dosing_IQRest(INPUT1=c(type="BOLUS")) )## Define dosing
# and one form of the model
modK=IQRmodel("mods/modK.txt")
cat(export_IQRmodel(modK))
modelSpecK=modelSpec_IQRest(POPvalues0=c(Ke=10,Vc=10,Kpc=0.05,Kcp=0.05), 
                             errorModel = list(OUTPUT1 = c(type = "rel", guess = 0.3) ))
estK <- IQRnlmeEst(model         = modK,
                    dosing        = dosing,
                    data          = data,
                    modelSpec     = modelSpecK)
                    
# finally, these last two function calls are backend specific                     
prNONK <- IQRnlmeProject(est = estK,
                         projectPath="projs/NONK",
                         tool="NONMEM",
                         algOpt.K1 = 50,
                         algOpt.K2 = 20,
                         comment = "BiExp Ke NONMEM version") 
run_IQRnlmeProject(prNONK) 

prMONK <- IQRnlmeProject(est = estK,
                         projectPath="projs/MONK",
                         tool="MONOLIX",
                         algOpt.K1 = 50,
                         algOpt.K2 = 20,
                         comment = "BiExp Ke  MONOLIX version")
run_IQRnlmeProject(prMONK) 

projMIXK <- IQRnlmeProject(est = estK,
                           projectPath="projs/MIXK",
                           tool="NLMIXR",
                           algOpt.K1 = 50,
                           algOpt.K2 = 20,
                           comment = "BiExp Ke  NLMIXR version")
run_IQRnlmeProject(projMIXK) 

png(filename="outs/NONK.png",width = 680, height = 680)
plotINDIV_IQRnlmeProject("projs/NONK",outputNr = 1,filename =NULL ,plotLog = TRUE,nindiv = 25)
dev.off()
png(filename="outs/MONK.png",width = 680, height = 680)
plotINDIV_IQRnlmeProject("projs/MONK",outputNr = 1,filename =NULL ,plotLog = TRUE,nindiv = 25)
dev.off()
png(filename="outs/MIXK.png",width = 680, height = 680)
plotINDIV_IQRnlmeProject("projs/MIXK",outputNr = 1,filename =NULL ,plotLog = TRUE,nindiv = 25)
dev.off()

###### now we update the model to its CL form
modCL=IQRmodel("mods/modCL.txt")
cat(export_IQRmodel(modCL))
modelSpecCL=modelSpec_IQRest(POPvalues0=c(CL=10,Vc=10,Q1=0.5,Vp1=10), 
                             IIVdistribution = c(CL="L",Vc="L",Q1="L",Vp1 ="L"),
                             errorModel = list(OUTPUT1 = c(type = "rel", guess = 0.3) ))
estCL <- IQRnlmeEst(model         = modCL,
                    dosing        = dosing,
                    data          = data,
                    modelSpec     = modelSpecCL)

prNONCL <- IQRnlmeProject(est = estCL,
                          projectPath="projs/NONCL",
                          tool="NONMEM",
                          algOpt.K1 = 50,
                          algOpt.K2 = 20,
                          comment = "BiExp CL NONMEM version") 
run_IQRnlmeProject(prNONCL) 

prMONCL <- IQRnlmeProject(est = estCL,
                          projectPath="projs/MONCL",
                          tool="MONOLIX",
                          algOpt.K1 = 50,
                          algOpt.K2 = 20,
                          comment = "BiExp CL  MONOLIX version")
run_IQRnlmeProject(prMONCL) 

projMIXR <- IQRnlmeProject(est = estCL,
                              projectPath="projs/MIXCL",
                              tool="NLMIXR",
                              algOpt.K1 = 50,
                              algOpt.K2 = 20,
                              comment = "BiExp CL  NLMIXR version")
run_IQRnlmeProject(projMIXR) 

png(filename="outs/NONCL.png",width = 680, height = 680)
plotINDIV_IQRnlmeProject("projs/NONCL",outputNr = 1,filename =NULL ,plotLog = TRUE,nindiv = 25)
dev.off()
png(filename="outs/MONCL.png",width = 680, height = 680)
plotINDIV_IQRnlmeProject("projs/MONCL",outputNr = 1,filename =NULL ,plotLog = TRUE,nindiv = 25)
dev.off()
png(filename="outs/MIXCL.png",width = 680, height = 680)
plotINDIV_IQRnlmeProject("projs/MIXCL",outputNr = 1,filename =NULL ,plotLog = TRUE,nindiv = 25)
dev.off()
```

Using the K parameterization of the two compartment model in mods/modK.txt
```
********** MODEL NAME
model_2cpt_linear_iv_Ke

********** MODEL NOTES

PK model for simulation of drug concentration in central compartment
with following characteristics:
Compartments:  2
Elimination :  linear
Absorption  :  none (infusion/iv)

Unit convention
Dose: mg
Concentration: ug/mL
Time: hours

The annotation of the parameter units is consistent with the given unit convention.
Units of the inputs (dose) and outputs (concentration) in the dataset for parameter estimation need to match the unit convention.

********** MODEL STATES

d/dt(Ac) 	=  - Kcp*Ac + Kpc*Ap1 - Ke*Ac + INPUT1
d/dt(Ap1) =    Kcp*Ac - Kpc*Ap1

Ac(0)    	= 0
Ap1(0)    	= 0

********** MODEL PARAMETERS

Vc 			  = 10 	# Central volume (L)
Ke 			  = 1 	# -alpha
Kpc			  = 0.05 	# - beta
Kcp			  = 0.05 	# 

********** MODEL VARIABLES

% Calculation of concentration in central compartment
Cc 			  = Ac/Vc

% Defining an output (only needed when interfacing with NLME
% parameter estimation tools such as NONMEM and MONOLIX)
OUTPUT1  	  = Cc  # Compound concentration (ug/mL)
********** MODEL REACTIONS
********** MODEL FUNCTIONS
********** MODEL EVENTS
```

yielded

for NONMEM

![](../../../docs/NONK.png)

for Monolix 

![](../../../docs/MONK.png)

and for nlmixr

![](../../../docs/MIXK.png)


Using  the CL-Q-Vp parameterization in mods/modCL.txt 

```
********** MODEL NAME
model_2cpt_linear_iv_CL

********** MODEL NOTES

PK model for simulation of drug concentration in central compartment
with following characteristics:
Compartments:  2
Elimination :  linear
Absorption  :  none (infusion/iv)

Unit convention
Dose: mg
Concentration: ug/mL
Time: hours

The annotation of the parameter units is consistent with the given unit convention.
Units of the inputs (dose) and outputs (concentration) in the dataset for parameter estimation need to match the unit convention.

********** MODEL STATES

d/dt(Ac) 	=  - Q1/Vc*Ac + Q1/Vp1*Ap1 - CL/Vc*Ac  + INPUT1
d/dt(Ap1) 	=  + Q1/Vc*Ac - Q1/Vp1*Ap1

Ac(0)    	= 0
Ap1(0)    = 0

********** MODEL PARAMETERS
CL 			  = 10 	# Clearance (L/hour)
Vc 			  = 10 	# Central volume (L)
Q1			  = 0.5 	# Intercompartmental clearance (L/hour)
Vp1			  = 10 	# Peripheral volume (L)
INPUT1 = 0
Tlag1 = 0

********** MODEL VARIABLES

% Calculation of concentration in central compartment
Cc 			  = Ac/Vc

% Defining an output (only needed when interfacing with NLME
% parameter estimation tools such as NONMEM and MONOLIX)
OUTPUT1  	  = Cc  # Compound concentration (ug/mL)

********** MODEL REACTIONS
********** MODEL FUNCTIONS
********** MODEL EVENTS
```


yielded 

for NONMEM

![](../../../docs/NONCL.png)

for Monolix

![](../../../docs/MONCL.png)


and for nlmixr 

![](../../../docs/MIXCL.png)



To compare results to those of Michor et al 2005, time is converted to days 
and limited to the first 365. 

```
head(d)
d1y=d%>%mutate(TIME=30.5*TIME)%>%filter(TIME<365)
(startvec <- c(lA = log(100), alpha = -0.05, lB=log(0.1), beta = -0.005))
nform <- ~ log(exp(lA+alpha*input) + exp(lB+beta*input))
nfun <- deriv(nform, namevec=c("lA", "alpha", "lB", "beta"),
              function.arg=c("input","lA", "alpha", "lB", "beta"))
nfun
(M<-nlmer(LNDV ~ nfun(TIME, lA, alpha, lB, beta) ~ (lA|ID) + (alpha|ID) +(lB|ID) + (beta|ID),
                        data=d1y,start = startvec,
      nAGQ = 0L,
      control = nlmerControl(tolPwrss = 1e-4)) )
(fe=fixef(M))
(RE=ranef(M)$ID)
apply(RE,2,mean) #check
(P=RE+t(t(rep(1,21)))%*%t(fe))
P=P%>%mutate(lA=exp(lA),lB=exp(lB))%>%rename(A=lA,B=lB)
P=lapply(1:21,function(i) unlist(P[i,]))
(d21y=d1y%>%group_by(ID)%>%summarize(Tf=max(TIME)))
d21y$Tpred=lapply(d21y$Tf,function(x) seq(0,x))

(D=mapply(simY,P,d21y$Tpred,1:21,SIMPLIFY = F))
(D=data.table(bind_rows(D)))
names(D)=c("TIME","LDV","ID")
mkDF=function(x) data.frame(TIME=250,LDV=2.5,a=round(x[["alpha"]],3),b=round(x[["beta"]],3))
(DTx=bind_rows(lapply(P,mkDF)))
DTx$ID=1:21
DTxb=DTx
DTxb$LDV=1
d1y%>%ggplot(aes(x=TIME,y=LDV))+facet_wrap(ID~.,ncol=4)+geom_point(size=1)+
  geom_line(data=D)+geom_text(aes(label=a),data=DTx,parse=T)+geom_text(aes(label=b),data=DTxb,parse=T)
ggsave("../docs/alphaNbetaFitsLME4d1y.png",width=7,height=8)

```

![](../../../docs/alphaNbetaFitsLME4d1y.png)
This shows that the rate constants all shrank to the population mean. 


Looking at the random effects
```
> RE
   lA alpha         lB beta
1   0     0 -1.1721089    0
2   0     0  0.3703322    0
3   0     0  0.2168800    0
4   0     0 -1.1838943    0
5   0     0 -0.1784630    0
6   0     0 -1.1395010    0
7   0     0 -0.8421898    0
8   0     0  0.5092889    0
9   0     0 -0.1917855    0
10  0     0 -1.5505433    0
11  0     0 -0.9056834    0
12  0     0  1.8923407    0
13  0     0  2.2577791    0
14  0     0  0.3605456    0
15  0     0  3.3665251    0
16  0     0 -1.4363420    0
17  0     0 -0.1196345    0
18  0     0 -1.5758449    0
19  0     0  2.0487526    0
20  0     0  0.6604085    0
21  0     0 -1.3868622    0
```
only the amplitude parameter B varies across patients.  
