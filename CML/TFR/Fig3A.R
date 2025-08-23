# Fig3A.R
setwd("~/cml/TFR")
library(patchwork)
library(tidyverse)
library(deSolve) 
source("setup.R")
mo=0;dur=6;base=2e-4; height=2e-3  #MR4 to MR3, L=2R/(R+1) = ~2R
Lh=height; Rh=Lh/(2-Lh)
(x0=c(L=base,X1=1,X2=0,X3b=0,H=0))
times <- seq(-2,12*6,by=0.1) #run from -2 months to 6 years
y0 <- ode(x0,times = times, func = markov,	parms = pars,method="lsodar")
d0=data.frame(y0) # baseline = continue dosing to hold L at 1e-4
(eventdat=data.frame(var="L",time=c(mo,mo+dur),value=c(height, base),method="replace"))
y1 <- ode(x0,times = times, func = markov,	parms = pars,events=list(data=eventdat),method="lsodar")
(d1=data.frame(y1)) # TFR attempt failure causes a 10-fold increase in L for 6 months
d0|>select(time:P)|>gather(key="Var",value="Value",-time)|>  
  mutate(Var=factor(Var,levels=c("L","X1","X2","X3b","H","P")))|>
  ggplot(aes(x=time,y=Value))+facet_grid(Var~.,scales = "free")+gl+xlab("Months")+tc(14)+sbb
d1|>select(time:P)|>gather(key="Var",value="Value",-time)|>
  mutate(Var=factor(Var,levels=c("L","X1","X2","X3b","H","P")))|>
  ggplot(aes(x=time,y=Value))+facet_grid(Var~.,scales = "free")+gl+xlab("Months")+tc(14)+sbb
(D0=d0%>%select(time,L,P)%>%mutate(Duration="0 Mos"))
(D1=d1%>%select(time,L,P)%>%mutate(Duration="6 Mos"))
D=bind_rows(D0,D1)
D=D%>%mutate(R=L/(2-L))
(p1=D%>%ggplot(aes(x=time/12,y=100*R,linetype=Duration))+geom_step()+sy+gxTFR+gyIS+cc1+tc(12)+top) 
print(DR0<-D%>%filter(time==mo+60,Duration=="0 Mos"))  # risk 5-years without stopping
print(DR1<-D%>%filter(time==mo+60,Duration=="6 Mos"))  # risk 5-years after stopping
print(progRisk <- DR1$P-DR0$P)  # 0.031% risk
(tit=paste0("TFR attempt adds ",round(100*progRisk,3),"% risk in 5y"))
(p2=D%>%ggplot(aes(x=time/12,y=P,linetype=Duration))+gl+ggtitle(tit)+tc(12)+ 
    theme(legend.position="none", plot.title=element_text(size=10,hjust=0.7,vjust=-5))+gxTFR+gyP )
p1/p2 + plot_layout(design = layout)
ggsave("outs/Fig3A.pdf",height=4.5,width=4)
