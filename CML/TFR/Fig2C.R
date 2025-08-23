## Fig2C.R
graphics.off();rm(list=ls())
setwd("~/cml/TFR")
library(tidyverse)
library(pracma)
source("setup.R")
del=1e-3
(tf=c(0,10^seq(log10(0.01),log10(200),del))) 
D100=he(tf)                # L=1
D10=he(tf,k1=0.035121676) # L=0.1
D1=he(tf, k1=0.0035121676 ) #L=0.01
Dp1=he(tf,k1=0.00035121676) #L=0.001
D100=D100%>%mutate(cg=cumtrapz(tf,h2),Sg=exp(-cg),yg=1-Sg)
D10=D10%>%mutate(cg=cumtrapz(tf,h2),Sg=exp(-cg),yg=1-Sg)
D1=D1%>%mutate(cg=cumtrapz(tf,h2),Sg=exp(-cg),yg=1-Sg)
Dp1=Dp1%>%mutate(cg=cumtrapz(tf,h2),Sg=exp(-cg),yg=1-Sg)
(t0=round(D100[[k<-which.min(abs(D100$yg-0.05)),"t"]],2))
(t1=round(D10[[k<-which.min(abs(D10$yg-0.05)),"t"]],2))
(t2=round(D1[[which.min(abs(D1$yg-0.05)),"t"]],2))
(t3=round(Dp1[[which.min(abs(Dp1$yg-0.05)),"t"]],2))
D100%>%ggplot(aes(x=tf,y=yg))+gl+gx+ gyP+
  geom_line(aes(x=tf,y=yg),col="black",data=D10)+ # add L=0.1 line
  geom_line(aes(x=tf,y=yg),col="black",data=D1)+  # add L=0.01 line 
  geom_line(aes(x=tf,y=yg),col="black",data=Dp1)+ # add L=0.001 line
  geom_hline(yintercept=0.05,col="gray")+
  geom_vline(xintercept=t0,col="gray")+
  geom_vline(xintercept=t1,col="gray")+
  geom_vline(xintercept=t2,col="gray")+
  geom_vline(xintercept=t3,col="gray")+
  scale_x_log10(breaks=c(t0,t1,t2,t3))+
  scale_y_continuous(breaks=c(0,0.025,0.05))+
  annotate("text", x = 0.3, y = 0.053, label ="Load = 100%")+ #, size = 12)+
  annotate("text", x = t1, y = 0.053, label ="10%        ",col="black")+ 
  annotate("text", x = t2, y = 0.053, label ="1%       ",col="black")+ 
  annotate("text", x = t3, y = 0.053, label ="0.1%         ",col="black")+ 
  tc(14)+coord_cartesian(xlim=c(.1,200),ylim=c(0,0.07)) #+
ggsave("outs/Fig2C.pdf",width=4.5,height=3.5)
