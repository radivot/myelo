# myelo
This R package houses models of myeloid hematopoiesis. To install it use:  
```
devtools::install_github("radivot/myelo",subdir="myelo")
```

<!--
# Pharmacokinetic and -dynamic modelling of G-CSF derivatives in humans

The model of   [Scholz et al  *Theoretical Biology and Medical Modelling* **9** 32 (2012)](https://www.ncbi.nlm.nih.gov/pubmed/22846180) 
-->

[//]: # comment 2

# Model of Chemotherapy-Induced Myelosuppression With Parameter Consistency Across Drugs

[Friberg et al  *J Clin Oncol* **20**  4713-4721 (2002)](https://www.ncbi.nlm.nih.gov/pubmed/12488418) 
provide a model of 5 neutrophil lineage cell state variables, one proliferating, 
3 transitioning through maturation stages, and one circulating in  blood. In R their model is 
```
friberg02<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dC1=-(k12+k13+k10)*C1 + k21*C2 + k31*C3
    dC2=k12*C1 - k21*C2
    dC3=k13*C1 - k31*C3
    Edrug=slope*C1
    dProl = ktr*Prol*(1-Edrug)*(Circ0/Circ)^gam - ktr*Prol
    dTrans1=ktr*Prol-ktr*Trans1
    dTrans2=ktr*Trans1-ktr*Trans2
    dTrans3=ktr*Trans2-ktr*Trans3
    dCirc=ktr*Trans3-ktr*Circ
    return(list(c(dC1,dC2,dC3,dProl,dTrans1,dTrans2,dTrans3,dCirc)))
  })
}
```
Using fits of this model to docetaxel response data, the parameter estimates in days are
```
fribergPars02=c(Circ0=5.05, ktr=24*4/88.7,gam=0.161,slope=8.58, 
		k12 = 24*1.06,k21 = 24*1.51,k13 = 24*1.26, k31 = 24*0.084,k10 = 24*5.2)
```
and initial state  is
```
x0=c(C1=0,C2=0,C3=0,Prol=5.05,Trans1=5.05,Trans2=5.05,Trans3=5.05,Circ=5.05)
```
in  uM and 1e3 cells/uL (of blood). 
The initial number of circulating neutrophils,  5.05,  is 
also the setpoint(i.e. Circ0) and thus steady state. 
At steady state this is also the number of cells in each of the 
three transition/maturation compartments and in the
proliferating compartment. A bolus  of docetaxel making C1 33.5 uM at t=0 yields a 
circulating neutrophil response similar to the curve PREDe 
in the upper left panel of Figure 4: at doses of 100 mg/m2 *2 m2 = 200 mg, 
the number of umoles is 247.524 = 200/0.808 (since 808 gm/mole=0.808 mg/umole),
so a volume of 7.4L => C1(0)=247.524/7.4 = 33.5 uM.
The code for this is
```
library(myelo)
times <- c(-5:25)
(evnt=data.frame(var="C1",time=0,value=33.5,method="rep"))
yout=ode(x0,times=times,func=friberg02,events=list(data=evnt),parms=fribergPars02)
plot(yout)
```
![](../docs/fri02.png)

```
library(tidyverse)
D=as.data.frame(yout)
tc=function(sz) theme_classic(base_size=sz)
gx=xlab("Days")
gy=ylab(quote(paste(10^3," Neutrophils/uL")))
D%>%ggplot(aes(x=time,y=Circ))+geom_line()+gx+gy+tc(14)+ylim(0,7)
ggsave("~/Results/myelo/friberg02.png",width=5,height=6)
```
![](../docs/friberg02.png)

