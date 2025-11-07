## Fig1A2A.R
graphics.off();rm(list=ls())#clear plots and environment 
setwd("~/cml/TFR")
library(tidyverse);library(SEERaBomb);library(ggsci)#load packages          
library(bbmle); library(broom)
source("setup.R")
if (0) {
  load("~/data/seer25/CML25.RData")
  d=d%>%select(yrdx,agedx,sex,surv,status)
  #see HMD section of https://github.com/radivot/SEERaBomb to make next file
  load("~/data/mrt/mrtUSA.RData")#loads US mortality data
  mrt[[1]]
  d=d%>%filter(yrdx<1986) #2732
  d=d%>%filter(agedx<90) # 2666
  sum(d$status) #2635 have died
  sum(d$status==0) # so 31 are still alive
  (D=msd(d,mrt,brkst=c(0,0.5,1,2,3,4,5,6,8,10,15,20,25)))
  D=D%>%rename(Group="sex")%>%select(Group,int,everything())
  (D=foldD(D)%>%mutate(Group="Pooled")%>%select(int,Group,everything()))
  D[,-(1:2)]=round(D[,-(1:2)],4)
  # datapasta::tribble_paste(D)
} # this code chunk makes this tibble, switch 0 to 1 (true) is you want to recreate it
D=tibble::tribble(   
  ~int,   ~Group,  ~O,      ~E,       ~PY,      ~t,   ~EAR,     ~LL,    ~UL,     ~RR,    ~rrL,    ~rrU,
  "(0,0.5]", "Pooled", 574, 37.1116, 1143.7083,  0.2147, 0.4694,  0.4284, 0.5105, 15.4669, 14.2273, 16.7856,
  "(0.5,1]", "Pooled", 320, 26.2885,   959.875,  0.7304,  0.306,  0.2695, 0.3425, 12.1726, 10.8753, 13.5821,
  "(1,2]", "Pooled", 435, 40.8485,  1557.625,  1.4414,  0.253,  0.2268, 0.2793, 10.6491,  9.6718, 11.6984,
  "(2,3]", "Pooled", 326, 29.2867, 1167.1667,  2.4381, 0.2542,  0.2239, 0.2845, 11.1313,  9.9557, 12.4077,
  "(3,4]", "Pooled", 255, 21.5045,     870.5,  3.4323, 0.2682,  0.2323, 0.3042,  11.858, 10.4471, 13.4063,
  "(4,5]", "Pooled", 207, 15.3585,   641.125,   4.427, 0.2989,  0.2549, 0.3429, 13.4779, 11.7043, 15.4443,
  "(5,6]", "Pooled", 139, 10.2999,  466.0417,  5.4282, 0.2762,  0.2266, 0.3257, 13.4953, 11.3452, 15.9345,
  "(6,8]", "Pooled", 169, 12.9088,  628.4167,  6.7772, 0.2484,  0.2078, 0.2889, 13.0919, 11.1924, 15.2213,
  "(8,10]", "Pooled",  74,  7.2761,  388.9167,  8.8276, 0.1716,  0.1282, 0.2149, 10.1703,  7.9859, 12.7679,
  "(10,15]", "Pooled",  81, 10.3884,  556.2917, 11.7211, 0.1269,  0.0952, 0.1586,  7.7971,  6.1921,  9.6911,
  "(15,20]", "Pooled",  24,  4.7688,    320.75, 17.0512,   0.06,    0.03, 0.0899,  5.0328,  3.2246,  7.4884,
  "(20,25]", "Pooled",   9,  4.3181,   256.125, 22.3613, 0.0183, -0.0047, 0.0412,  2.0843,  0.9531,  3.9566,
  "(25,100]", "Pooled",  22,  8.6815,  513.7083, 30.7376, 0.0259,   0.008, 0.0438,  2.5341,  1.5881,  3.8367
)

nLL=function(lk1,lk2,lk3,lk4,lp1,lp2) {  # nLL = negative Log Likelihood
  k1=exp(lk1);k2=exp(lk2);k3=exp(lk3);k4=exp(lk4);p1=exp(lp1);p2=exp(lp2)
  t=D$t
  x3A=p1*exp(-k3*t)  
  x1=p2*exp(-k1*t)  
  x2=k1*p2*(exp(-k1*t)-exp(-k2*t))/(k2-k1)
  x3B=k1*k2*p2*(exp(-k1*t)/((k2-k1)*(k3-k1))+exp(-k2*t)/((k1-k2)*(k3-k2))+exp(-k3*t)/((k1-k3)*(k2-k3)))
  x4=(1-p1-p2)*exp(-k4*t)  
  denom=(x1+x2+x3A+x3B+x4)
  h1=k3*x3A/denom
  h2=k3*x3B/denom
  h3=k4*x4/denom
  he=h1+h2+h3
  he[he<=0]=1e-6
  -sum(dpois(D$O, lambda=he*D$PY+D$E, log=TRUE))
}
fit=1:6
IC=c(lk1=log(0.3), lk2=log(1.2), lk3=log(1.4), lk4=log(0.03),lp1=log(0.33),lp2=log(0.62))
(s1=summary(m1<-mle2(nLL,
                     method="Nelder-Mead",
                     fixed=as.list(IC)[-fit],
                     start=as.list(IC), 
                     control = list(maxit=50000, parscale=IC[fit]) ) ))
(M=coef(s1))
(CI=cbind(exp(M[,1]),exp(M[,1]-1.96*M[,2]),exp(M[,1]+1.96*M[,2])))
(nms=str_replace(rownames(CI),"l",""))
rownames(CI)<-nms
CI
#### april 2025 values
# k1 0.35121676 0.29636703 0.41621772
# k2 1.11192339 0.69989370 1.76651628
# k3 1.61921238 1.21894520 2.15091600
# k4 0.02503601 0.01380432 0.04540621
# p1 0.35909041 0.30104431 0.42832872
# p2 0.58603095 0.52034433 0.66000963
p1=0.35909041
p2=0.58603095
1-p1-p2 # 0.0549 
(pars=as.list(CI[,1])) # save full res for simulations
(CIk=apply(CI,2,round,c(3,3,3,4,3,3))) # round to make text string for MS
paste0(nms," = ",CIk[,1],"(",CIk[,2],", ",CIk[,3],")",collapse=", ") 
#k1 = 0.351(0.296, 0.416), k2 = 1.112(0.7, 1.767), k3 = 1.619(1.219, 2.151), 
# k4 = 0.025(0.0138, 0.0454), p1 = 0.359(0.301, 0.428), p2 = 0.586(0.52, 0.66)
(CIt=apply(1/CI,2,round,c(3,3,3,1,2,2)))
(nmsT=str_replace(nms,"k","Tau"))
paste0(nmsT," = ",CIt[,1]," (",CIt[,3],", ",CIt[,2],")",collapse=", ")
# "Tau1 = 2.847 (2.403, 3.374), Tau2 = 0.899 (0.566, 1.429), Tau3 = 0.618
# (0.465, 0.82), Tau4 = 39.9 (22, 72.4), p1 = 2.78 (2.33, 3.32), p2 = 1.71 (1.52, 1.92)"
# redefine function for fine grain simulation and showing all 3 components
he=function(t,k1=0.35121676,k2=1.11192339,k3=1.61921238,k4=0.02503601,p1=0.35909041,p2=0.58603095) {
  x3A=p1*exp(-k3*t)  
  x1=p2*exp(-k1*t)  
  x2=k1*p2*(exp(-k1*t)-exp(-k2*t))/(k2-k1)
  x3B=k1*k2*p2*(exp(-k1*t)/((k2-k1)*(k3-k1))+exp(-k2*t)/((k1-k2)*(k3-k2))+exp(-k3*t)/((k1-k3)*(k2-k3)))
  x4=(1-p1-p2)*exp(-k4*t)  
  denom=(x1+x2+x3A+x3B+x4)
  h1=k3*x3A/denom
  h2=k3*x3B/denom
  h3=k4*x4/denom
  he=h1+h2+h3
  tibble(t,h1,h2,h3,he,x3A,x1,x2,x3B,x4)
}
del=0.02
tf=seq(0,30,del) #fine grid of times to 30 years
tpars=c(list(t=tf),pars)
(D2=do.call(he,tpars))
D2%>%ggplot(aes(x=t,y=h1))+gx+gyE+
  geom_line(linewidth=.5,col="gray50",linetype="dashed")+
  geom_line(aes(y=h2),linewidth=.5,col="gray50")+
  geom_line(aes(y=h3),linewidth=.5,col="gray50",linetype="dotted")+
  geom_line(aes(y=he),linewidth=.75,col="black")+
  geom_point(aes(x=t,y=EAR),data=D,color="black")+ 
  geom_errorbar(aes(x=t,y=EAR,ymin=LL,ymax=UL),width=.2,data=D,color="black")+
  tc(14)
ggsave("outs/Fig1A.pdf",width=4.5,height=3.5)

tpars1=tpars
tpars1$k1=0.01*tpars$k1  # L=0.01
tpars1$p1=0.0
tpars1$p2=1.0  # all in CP
(D01=do.call(he,tpars1)) #for L=0.01
tpars1$k1=tpars$k1  #now bring L back up to 1
(D2=do.call(he,tpars1)) #for L=1, i.e. 2 as in 10^2 %
library(pracma)
(D01=D01%>%mutate(c01=cumtrapz(t,h2),S01=exp(-c01),y01=1-S01))
D2=D2%>%mutate(c2=cumtrapz(t,h2),S2=exp(-c2),y2=1-S2) 
D2%>%ggplot(aes(x=t,y=y2))+geom_line()+gx+ gyP+
  geom_line(aes(x=t,y=y01),col="black",data=D01)+
  geom_hline(yintercept=0.05,col="gray")+
  geom_vline(xintercept=1,col="gray")+
  geom_vline(xintercept=16,col="gray")+
  scale_x_continuous(breaks=c(0,1,5,10,16,20,30))+
  scale_y_continuous(breaks=c(0,0.05,0.25,0.50,0.75,1.0))+
  annotate("text",6.0,y=0.25, label ="100% Load")+ 
  annotate("text",x=25,y=0.13,label="1% Load",col="black")+ 
  tc(14)
ggsave("outs/Fig2A.pdf",width=4.5,height=3.5)
