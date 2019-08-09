############################  Friberg02 ######################
library(tidyverse)
library(deSolve)
library(myelo)
fribergPars02
# Figure 4 upper left panel (docetaxel), want nadir at 7 days and repeak at 20 day
(d=data.frame(x=c(0,6,9,13,16,20,22),y=c(5.4,2,0.6,1,1.7,7,8.6))) # mouse-balling it
plot(d)
times <- c(-5:25)
200/0.808/7.4
(evnt=data.frame(var="C1",time=0,value=33.5,method="rep"))
x0=c(C1=0,C2=0,C3=0,Prol=5.05,Trans1=5.05,Trans2=5.05,Trans3=5.05,Circ=5.05)
yout=ode(x0,times=times,func=friberg02,events=list(data=evnt),parms=fribergPars02)
plot(yout) #use export in gui


D=as.data.frame(yout)
tc=function(sz) theme_classic(base_size=sz)
gx=xlab("Days")
gy=ylab(quote(paste(10^3," Neutrophils/uL")))
D%>%ggplot(aes(x=time,y=Circ))+geom_line()+gx+gy+tc(14)+ylim(0,7)
ggsave("~/Results/myelo/friberg02.png",width=5,height=6)

D$Circ
min(D$Circ)
