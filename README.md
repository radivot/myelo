# myelo
This R package houses models of myeloid hematopoiesis. To install it use:  
```
devtools::install_github("radivot/myelo",subdir="myelo")
```
The help pages can then be reached via:

help(pack="myelo")


## A mathematical model for chronic myelogenous leukemia (CML) and T cell interaction (Moore and Li 2004)

This model captures CML cells (C) interacting with naive T cells (Tn) and effector T cells (Te). 
The differential equations of this model are:
![](docs/moore04model.png)


The following code reproduces Figures 6-8 in [Moore and Li 2004](https://www.ncbi.nlm.nih.gov/pubmed/15038986) 

![](docs/fig6to8.png)


```
library(myelo)  
library(deSolve)
library(tidyverse)
#Table 1
(p1=c(sn=0.073,dn=0.04,de=0.06,dc=0.2,  kn=.001,eta=100,alfn=0.41,alfe=0.2, Cmax=3e5, rc=.03, ge=.005, gc=.005)) 
#Table 3
(p6=c(sn=0.37, dn=0.23,de=0.30,dc=0.024,kn=.062,eta=720,alfn=0.14,alfe=0.98,Cmax=23e4,rc=.0057,ge=.057,gc=.0034)) 
(p7=c(sn=0.29, dn=0.35,de=0.40,dc=0.012,kn=.066,eta=140,alfn=0.39,alfe=0.65,Cmax=16e4,rc=.011, ge=.079, gc=.058)) 
(p8=c(sn=0.071,dn=0.05,de=0.12,dc=0.68, kn=.063,eta=43, alfn=0.56,alfe=0.53,Cmax=19e4,rc=.23,  ge=.0077,gc=.047)) 
(P=list(Table1=p1,Fig.6=p6,Fig.7=p7,Fig.8=p8))
(ic=c(Tn=1510,Te=20,C=1e4)) # units are cells/uL 
t= seq(0,750,1)             #  and days
######## see http://eriqande.github.io/2015/01/22/solving-differential-equations-in-R.html
L=lapply(P, function(x) { 
  params=x
  ode(y=ic,times=t,func=moore04,parms=params)
})
L
D=do.call(rbind, lapply(names(L), function(x) data.frame(L[[x]],x,stringsAsFactors=F)))
###########
ltb=theme(legend.margin=margin(0,0,0,0),legend.title=element_blank())
ltp=theme(legend.position="top",legend.direction="horizontal")
tc=function(sz) theme_classic(base_size=sz);
gy=ylab("CML Cells/uL")
gx=xlab("Days")
D%>%ggplot(aes(x=time,y=C,color=x))+geom_line(size=1)+gx+gy+tc(14)+ltb+ltp
ggsave("~/Results/CML/fig6to8.png",height=3,width=3.5)
```



The code below further zooms in on 0 to 100 of the plot in Figure 8

```
cc=coord_cartesian(xlim=c(0,100))
library(scales)
D%>%filter(x=="Fig.8")%>%
  ggplot(aes(x=time,y=C))+geom_line(size=1,color=hue_pal()(4)[3])+gx+gy+tc(14)+ltb+ltp+cc+
  annotate("text", x =50, y =10000, label = "Fig. 8")
ggsave("~/Results/CML/fig8.png",height=3,width=3)
```
![](docs/fig8.png)


We can now examine time courses of naive T cells
```
gy=ylab("Naive T Cells/uL")
D%>%ggplot(aes(x=time,y=Tn,color=x))+geom_line(size=1)+gx+gy+tc(14)+ltb+ltp+cc
ggsave("~/Results/CML/naive.png",height=3,width=3.5)
```
![](docs/naive.png)

and effector T cells
```
gy=ylab("Effector T Cells/uL")
D%>%ggplot(aes(x=time,y=Te,color=x))+geom_line(size=1)+gx+gy+tc(14)+ltb+ltp+cc
ggsave("~/Results/CML/effector.png",height=3,width=3)
```
![](docs/effector.png)

It seems the initial effecter T cell population is far higher than it should be. The following code zooms in on the first 3 time points
```
D%>%group_by(x)%>%nest()%>%mutate(top=map(data,function(x) x[1:3,]))%>%unnest(top)
#  A tibble: 12 x 6
#    x       time    Tn      Te      C     T
#    <chr>  <dbl> <dbl>   <dbl>  <dbl> <dbl>
#  1 Table1     0 1510  20      10000  1530 
#  2 Table1     1 1449.  0.0130  9062. 1449.
#  3 Table1     2 1391.  0.0137  8252. 1391.
#  4 Fig.6      0 1510  20      10000  1530 
#  5 Fig.6      1 1133.  0.0162  9937. 1133.
#  6 Fig.6      2  850.  0.0122  9877.  850.
#  7 Fig.7      0 1510  20      10000  1530 
#  8 Fig.7      1  997.  0.0316 10148.  997.
#  9 Fig.7      2  659.  0.0205 10319.  659.
# 10 Fig.8      0 1510  20      10000  1530 
# 11 Fig.8      1 1349.  0.645   9606. 1350.
# 12 Fig.8      2 1205.  0.588   9415. 1206.
```

A help page is available for the function called by ode() of the R package deSolve, i.e. moore04(). Its definition is
```
moore04<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)),{
    dTn = sn - dn*Tn - kn*Tn*C/(C+eta) 
    dTe = alfn*kn*Tn*C/(C+eta) + alfe*Te*C/(C+eta) - de*Te - ge*C*Te
    dC  = rc*C*log(Cmax/C) - dc*C - gc*C*Te 
    list( c(dTn,dTe,dC), c(T=Tn+Te) )
  }) 
}
```
## A computational model to understand mouse iron physiology and disease [(Parmar and Mendes 2019, Plos Comp Biol)](https://www.ncbi.nlm.nih.gov/pubmed/30608934) 

This code first runs the model from -30 to 0 days to show the wildtype steady state. 
It then emulates hemachromatosis formation by setting
the hepcidin synthesis rate constant (ksHepci) to 0 at time 0 before running out to 365 days to yield values in Table 1.
The  set of plots below this code are from the plot method provided by deSolve: in the last line
the first argument to plot() is a matrix of class deSolve, so the plot method plot.deSolve() gets invoked. 

```
# use COPASI to save parmar sbml as XPPAUT *.ODE and do rest by hand + find/replace
library(myelo)   # load definition of function parmar19
library(deSolve)
library(rodeoExt) #provides rbind for deSolve class matrices

ic=c(NTBI=5.2e-11,Tf=1.50934e-08,FeDuo=3.85469e-07,FeBM=4.45824e-07,FeRest=7.91257e-06,FeLiver=1.51048e-06,
     Fe1Tf=1.05653e-08,FeSplee=1.52679e-07,EPO=3.15471e-15,Hepcidi=2.9912e-11,FeRBC=1.817e-05,Fe2Tf=2.46513e-08) 
(parameters=parmarPars19)

out=ode(y = ic, times = seq(-30,0,1), func = parmar19, parms = parameters) 
(N=length(ic))
n=dim(out)[1]
X0=out[n,2:(N+1)]
names(X0)<-names(ic)
parameters["ksHepci"]=0
out2   <- ode(y=X0, times=seq(1, 365, by = 1), func = parmar19, parms = parameters)
out=rbind(out,out2)
head(out)
graphics.off()
quartz(width=8,height=7)
vars=c("NTBI_c","Tf_c","FeDuo_c","FeBM_c","FeRBC_c","FeLiver_c","FeSplee_c","FeRest_c","Hepcidi_c")
plot(out,which=vars,xlab="Days",ylab="Concentration (M)") #screen capture to png
```
![](docs/parmar19base.png)

A newer way to plot things is to use ggplot2. For this we first need to convert the data into a long format using gather().
Using ggplot2 then yields nicer y-axes and also has the convenience of writing the file to a png automatically.

```
library(tidyverse)
ltb=theme(legend.margin=margin(0,0,0,0),legend.title=element_blank())
ltp=theme(legend.position="top",legend.direction="horizontal")
tc=function(sz) theme_classic(base_size=sz)
sbb=theme(strip.background=element_blank())
gy=ylab("Concentration in M")
gx=xlab("Days")
D=data.frame(out)
head(D)
D=D%>%select(time,NTBI_c:Hepcidi_c)
d=D%>%gather(key=variable,value=Concentration,-time)
d%>%ggplot(aes(x=time,y=Concentration))+facet_wrap(~variable,scales="free",nrow=3)+geom_line(size=1)+gx+gy+tc(14)+ltb+ltp+sbb
ggsave("~/Results/myelo/parmar19.png",height=6,width=6.5)
```

![](docs/parmar19.png)

The transferrin concentration (Tf_c) is seen to dip to negative values.  COPASI does not do this.

![](docs/copasiTfC.png)

The reason for this difference is unclear. 


Fortunately, there is now a COPASI R Connector R package called CoRC, see https://github.com/jpahle/CoRC. CoRC allows us to run  the iron model through COPASI directly via model files created by COPASI, which was used to develop the model. The code (below) is straightforward. Plots are made using ggplot2 as above. Note that the runTC() (run Time Course)  argument save_result_in_memory must be overriden to TRUE.  

```
library(CoRC)
path="~/ccf/jarek/grants/msb/iron/parmar19sup/cps/"
(m0=loadModel(paste0(path,"IronMousePV3.cps")))
r0=runTC( model=m0,duration=30,save_result_in_memory = TRUE)
D0=r0$result%>%select(Time:FeRBC)%>%select(-Fe1Tf,-EPO)
D0$Time=D0$Time-30
(m1=loadModel(paste0(path,"IronMousePV3_Hemochromatosis.cps")))
r1=runTC(model=m1,save_result_in_memory = TRUE) 
D1=r1$result%>%select(Time:FeRBC)%>%select(-Fe1Tf,-EPO)
D=bind_rows(D0,D1)
d=D%>%gather(key=variable,value=Concentration,-Time)
d%>%ggplot(aes(x=Time,y=Concentration))+facet_wrap(~variable,scales="free",nrow=3)+geom_line(size=1)+gx+gy+tc(14)+ltb+ltp+sbb
ggsave("~/Results/myelo/parmar19CoRC.png",height=6,width=6.5)
```

![](docs/parmar19CoRC.png)
The high frequency glitches at early times of the perturbation of zapping the hepicidin synthesis rate instantly to zero have now dissappeared. 
