## Fig2C.R
graphics.off();rm(list=ls())
setwd("~/cml/TFR")
library(tidyverse)
library(pracma)
source("setup.R")

del=1e-6  #good as it comes for this 
(tf=c(seq(0,0.69,0.01),10^seq(log10(0.7),log10(200),del)))

(x=10^0.05) # 1.12 => x^20 = 10 (20 samples per decade)
(Ls=0.001*x^(0:60)) #all same ratio apart
(L=lapply(Ls,function(x) he(tf,k1=x*0.35121676)))
(L=lapply(L,function(x) x%>%mutate(cg=cumtrapz(t,h2),Sg=exp(-cg),yg=1-Sg)))
(L=lapply(L,function(x) x[which.min(abs(x$yg-0.05)),]))
(D=bind_rows(L))
D$Ls=Ls

gx=xlab("Load (%)")
gy=ylab("Years to 5% Probability of Death by CML")
(yL=-148*0.05)
myt=theme(axis.title.y=element_text(size=12))

D%>%ggplot(aes(x=Ls*100,y=t))+gl+gx+ gy+tc(14)+
  scale_x_log10(breaks=c(0.1,1,10,100))+
  geom_segment(aes(x=0.1,y=yL,xend=0.1,yend=148),col="gray")+
  geom_segment(aes(x=0.1,y=148,xend=1,yend=148),col="gray")+
  annotate("text", x = 1.8, y = 148, label ="148 y",col="black",size=5)+ 
  annotate("text", x = 17, y = 16, label ="16 y",col="black",size=5)+ 
  annotate("text", x = 72, y = 8, label ="3 y",col="black",size=5)+ 
  annotate("text", x = 140, y = 1, label ="1 y",col="black",size=5)+ 
  geom_segment(aes(x=1,y=yL,xend=1,yend=16),col="gray")+
  geom_segment(aes(x=1,y=16,xend=10,yend=16),col="gray")+
  geom_segment(aes(x=10,y=yL,xend=10,yend=3),col="gray")+
  geom_segment(aes(x=10,y=3,xend=50,yend=3),col="gray")+
  geom_segment(aes(x=100,y=yL,xend=100,yend=1),col="gray")+
  scale_y_continuous(breaks=c(0,16,25,50,75,100,148),expand = c(0,4))+myt
ggsave("outs/Fig2C.pdf",width=4.5,height=3.5)
