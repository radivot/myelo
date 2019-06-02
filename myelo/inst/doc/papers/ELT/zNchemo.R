############################ ZHUGE12 N only CHEMO only  ##################
# This code checks out fine with Figure 2b and code without events in gcsfDemo.R
library(tidyverse)
library(deSolve)
library(myelo)  
zhugePars["Nss"]=639805537  #SS found in readme.md
Tf=200
mkEventsC=function(zhugePars,Tf) {  #chemo events only
  (Tc1=seq(0,Tf,zhugePars["T"]))
  N1=length(Tc1)
  (Tc2=seq(1,Tf,zhugePars["T"]))
  N2=length(Tc2)
  (eventdat1 <- data.frame(var = rep("eta",N1),
                           time = Tc1 ,
                           value = rep(zhugePars["etaMinNP"],N1),
                           method = rep("rep",N1)))
  (eventdat2 <- data.frame(var = rep("eta",N2),
                           time = Tc2 ,
                           value = rep(zhugePars["etaNP"],N2),
                           method = rep("rep",N2)))
  bind_rows(eventdat1,eventdat2)%>%arrange(time)
}

eventDF=mkEventsC(zhugePars,Tf)
head(eventDF)
sapply(eventDF,class)


times= seq(-zhugePars[["tauN"]],Tf,by=0.01)

# Chemo acts only on proliferation,  so add one state to integrate eta 
zhuge12Nchemo<-function(Time, State, Pars) {  # model with stem cells Q treated as constant
  with(as.list(c(State, Pars)), {
    deta=0
    dEta=eta
    if (Time < 0) {
      An=exp(etaNP*tauNP-gam0*tauNM)  # no gcsf or chemo perturbations for negative times
      dN=-gamN*N + An*f0/(1+(Nss/the1)^s1)*Qss
    }
    else{
      delEta=lagvalue(Time - tauNM,3)-lagvalue(Time - tauN,3)
      An=exp(delEta - gam0*tauNM)
      dN=-gamN*N + An*f0/(1+(lagvalue(Time - tauN,1)/the1)^s1)*Qss
    }
    list(c(dN,deta,dEta))
  })
}



yout=dede(c(N=zhugePars[["Nss"]],eta=zhugePars[["etaNP"]],Eta=0),
          times=times,func=zhuge12Nchemo,
          parms=zhugePars,events=list(data=mkEventsC(zhugePars,Tf)),method="lsodar")
D=data.frame(yout)
tc=function(sz) theme_classic(base_size=sz)
gy=ylab("Neutrophil Counts")
gx=xlab("Days")
D%>%ggplot(aes(x=time,y=N))+geom_line(size=1)+gx+gy+tc(14) 

zhugePars["T"]=18
yout=dede(c(N=zhugePars[["Nss"]],eta=zhugePars[["etaNP"]],Eta=0),
          times=times,func=zhuge12Nchemo,
          parms=zhugePars,events=list(data=mkEventsC(zhugePars,Tf)),method="lsodar")
D18=data.frame(yout)
D18%>%ggplot(aes(x=time,y=N))+geom_line(size=1)+gx+gy+tc(14) 

zhugePars["T"]=23
yout=dede(c(N=zhugePars[["Nss"]],eta=zhugePars[["etaNP"]],Eta=0),
          times=times,func=zhuge12Nchemo,
          parms=zhugePars,events=list(data=mkEventsC(zhugePars,Tf)),method="lsodar")
D23=data.frame(yout)
D23%>%ggplot(aes(x=time,y=N))+geom_line(size=1)+gx+gy+tc(14) 
sy=scale_y_log10()
D23%>%ggplot(aes(x=time,y=N))+geom_line(size=1)+gx+gy+tc(14)+sy
D23$Tc="23 Days"
D18$Tc="18 Days"
D=bind_rows(D18,D23)%>%mutate(N=N/1e8)
ltp=theme(legend.position="top")
D%>%ggplot(aes(x=time,y=N,col=Tc))+geom_line(size=1)+gx+gy+tc(14)+sy+ltp
cc=coord_cartesian(ylim=c(1e-2,1e4))
gh=geom_hline(yintercept=0.63)
D%>%ggplot(aes(x=time,y=N,col=Tc))+geom_line(size=1)+gx+gy+tc(14)+sy+ltp+cc+gh
ggsave("~/Results/myelo/zhugeNchemoEventsfig2B.png",height=6,width=6.5)
