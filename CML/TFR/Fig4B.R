# Fig4B.R  
setwd("~/cml/TFR")
library(patchwork)
library(tidyverse)
library(deSolve) 
source("setup.R")
markovHahnel<-function(Time, State, Pars) {
  with(as.list(c(Time, State, Pars)),{
    L=Y/Ky
    dUtki   = 0
    dUpz   = 0
    dX1  = -k1*L*X1
    dX2  = k1*L*X1 - k2*X2
    dX3b = k2*X2-k3*X3b # H is the cumulative hazard of progression, 
    dH   = k3*X3b/(X1+X2+X3b) # conditional on it not having happened yet 
    dX = -pxy*X + pyx*Y 
    dY =  pxy*X - pyx*Y + py*Y*(1-Y/Ky)   -  m*Y*Z - TKI*Utki*Y
    dZ =  rz    -   az*Z                   + pz*Upz*Y*Z/(Kz^2+Y^2) 
    list(c(dUtki,dUpz,dX1,dX2,dX3b,dH,dX,dY,dZ),c(P=1-exp(-H),R=L/(2-L)))
  })                                # P=1-RS where RS=exp(-H)
}
mpars=c(k1=0.35121676,k2=1.11192339,k3=1.61921238)/12 # updated to April 2025 values 
## bring in the smoothing model parameters fitted in Fig4A.R
hpars=c(pyx = 2.891981305438, pxy = 0.117493735003685, TKI = 13.3427401229456,
        pz = 4373.5, Kz = 632.5, py = 1.658, Ky = 1e+06, az = 2, rz = 200,
        m = 1e-04, X0 = 688039.114941631, Y0 = 27953.2531178926, Z0 = 108.482082597056,
        CessT = 69.9846394984326) # hahnel model parameters from Fig4A.R
# make new Markov-Hahnel (MH)  params and ICs
(mhPars=c(mpars,hpars[1:10]))
(mhIC0=c(Utki=1,Upz=1,X1=1,X2=0,X3b=0,H=0,X = 771.745489104474, Y = 6.26564295379477, Z = 103.729321625289))
(mhIC1=c(Utki=0,Upz=1,X1=1,X2=0,X3b=0,H=0,X = 771.745489104474, Y = 6.26564295379477, Z = 103.729321625289))

times <- seq(-12*0,12*10,by=1) #run from 0 months to 10 more years
# no TFR attempt, skip immunosuppression event since not relevant 
y0 <- ode(mhIC0,times = times, func = markovHahnel,	parms = mhPars,method="lsodar")
(d0=data.frame(y0)) 
d0|>select(time:R)|>gather(key="Var",value="Value",-time)|>
  mutate(Var=factor(Var,levels=c("Utki","Upz","X1","X2","X3b","H","P","R")))|>
  ggplot(aes(x=time,y=Value))+facet_grid(Var~.,scales = "free")+geom_line(size=1)+xlab("Months")+tc(14)+sbb

# now do TFR from time of attempt
(evPz=data.frame(var="Upz",time=c(25,35),value=c(0, 1),method="replace"))
# last measurement at 112 is under 0.1%, so assume TKI restarted at 115, i.e 45 after 70
# (evTKI=data.frame(var="Utki",time=c(45,57),value=c(1,0),method="replace")) #try TFR2 1y later
(evTKI=data.frame(var="Utki",time=c(45),value=c(1),method="replace")) #keep dosing after TFR1 failure
(eventdat=rbind(evPz,evTKI))
y1 <- ode(mhIC1,times = times, func = markovHahnel,	parms = mhPars,events=list(data=eventdat),method="lsodar")
(d1=data.frame(y1)) 
d1|>select(time:R)|>gather(key="Var",value="Value",-time)|>
  mutate(Var=factor(Var,levels=c("Utki","Upz","X1","X2","X3b","H","P","R")))|>
  ggplot(aes(x=time,y=Value))+facet_grid(Var~.,scales = "free")+geom_line(size=1)+xlab("Months")+tc(14)+sbb

(D0=d0%>%select(time,R,P)%>%mutate(TFR="No"))
(D1=d1%>%select(time,R,P)%>%mutate(TFR="Yes"))
D=bind_rows(D0,D1)
bks=c(0.001,0.01,0.1)
sy=scale_y_log10(breaks=bks,labels=bks)
cc1=coord_cartesian(ylim=c(1e-4,1))
(p1=D%>%ggplot(aes(x=time/12,y=100*R,linetype=TFR))+gl+sy+gxTFR+gyIS+cc1+tc(12)+top) 

print(DR0<-D%>%filter(time==120,TFR=="No"))  # risk 10-years without stopping
print(DR1<-D%>%filter(time==120,TFR=="Yes"))  # risk 10-years after stopping
print(progRisk <- DR1$P-DR0$P) # 0.0001555865
(tit=paste0("TFR attempt adds ",round(100*progRisk,3),"% risk in 10y"))
(p2=D%>%ggplot(aes(x=time/12,y=P,linetype=TFR)) + geom_line()+ ggtitle(tit)+ tc(12)+ 
    theme(legend.position="none", plot.title=element_text(size=10,hjust=0.7,vjust=-4))+gxTFR+gyP )
p1/p2 + plot_layout(design = layout)
ggsave("outs/Fig4B.pdf",height=4.5,width=4)
