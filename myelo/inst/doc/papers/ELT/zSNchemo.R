############################ ZHUGE12 SN + CHEMO only  ##################
# This code checks out fine with Figure 2b and code without events in gcsfDemo.R
library(tidyverse)
library(deSolve)
library(myelo)  
zhugePars["Nss"]=639805537  #SS found in readme.md
Tf=200
mkEventsC=function(zhugePars,Tf) {  #chemo events only
  # (Tc1=seq(zhugePars["T"],Tf,zhugePars["T"]))
  (Tc1=seq(0,Tf,zhugePars["T"]))
  Nc1=length(Tc1)
  # (Tc2=seq(zhugePars["T"]+1,Tf,zhugePars["T"]))
  (Tc2=seq(1,Tf,zhugePars["T"]))
  Nc2=length(Tc2)
  (etacB <- data.frame(var = rep("eta",Nc1),
                       time = Tc1 ,
                       value = rep(zhugePars["etaMinNP"],Nc1),
                       method = rep("rep",Nc1)))
  (etacE <- data.frame(var = rep("eta",Nc2),
                       time = Tc2 ,
                       value = rep(zhugePars["etaNP"],Nc2),
                       method = rep("rep",Nc2)))
  (gamsB <- data.frame(var = rep("gams",Nc1),
                       time = Tc1 ,
                       value = rep(zhugePars["gamMaxS"],Nc1),
                       method = rep("rep",Nc1)))
  (gamsE <- data.frame(var = rep("gams",Nc2),
                       time = Tc2 ,
                       value = rep(zhugePars["gamS"],Nc2),
                       method = rep("rep",Nc2)))
  # bind_rows(eventdat1,eventdat2)%>%arrange(time)
  rbind(etacB,etacE,gamsB,gamsE)%>%arrange(time) #no warnings
}

eventDF=mkEventsC(zhugePars,Tf)
head(eventDF)
sapply(eventDF,class)


times= seq(-zhugePars[["tauN"]],Tf,by=0.01)

# Chemo acts only on proliferation,  so add one state to integrate eta 
zhuge12SNchemo<-function(Time, State, Pars) {  # model with stem cells Q treated as constant
  with(as.list(c(State, Pars)), {
    deta=0
    dgams=0
    if (Time < 0) {
      Aq=2*exp(-gamS*tauS)
      dQ=-(k0/(1+(Q/the2)^s2) + f0/(1+(N/the1)^s1))*Q + Aq*k0/(1+(Qss/the2)^s2)*Qss
      An=exp(etaNP*tauNP-gam0*tauNM)  # no gcsf or chemo perturbations for negative times
      dN=-gamN*N + An*f0/(1+(Nss/the1)^s1)*Qss
      dEta=0
      dGams=0
    }
    else{
      dEta=eta
      dGams=gams
      delEta=lagvalue(Time - tauNM,4)-lagvalue(Time - tauN,4)
      An=exp(delEta - gam0*tauNM)
      dN=-gamN*N + An*f0/(1+(lagvalue(Time - tauN,1)/the1)^s1)*Qss
      delGamS=Gams -lagvalue(Time - tauS,5)
      Aq=2*exp(-delGamS)
      Qts=lagvalue(Time - tauS,1)
      Qtn=lagvalue(Time - tauN,1)
      Ntn=lagvalue(Time - tauN,2)
      dQ=-(k0/(1+(Q/the2)^s2) + f0/(1+(N/the1)^s1))*Q + Aq*k0/(1+(Qts/the2)^s2)*Qts
      dN=-gamN*N + An*f0/(1+(Ntn/the1)^s1)*Qtn
    }
    # list(c(dN,deta,dEta))
    list(c(dQ,dN,deta,dEta,dgams,dGams))
  })
}


yout=dede(c(Q=zhugePars[["Qss"]],N=zhugePars[["Nss"]],eta=zhugePars[["etaNP"]],Eta=0,
            gams=zhugePars[["gamS"]],Gams=0),
          times=times,func=zhuge12SNchemo,
          parms=zhugePars)
tail(yout,4)
head(yout,4)

yout=dede(c(Q=zhugePars[["Qss"]],N=zhugePars[["Nss"]],eta=zhugePars[["etaNP"]],Eta=0,
                                                     gams=zhugePars[["gamS"]],Gams=0),
          times=times,func=zhuge12SNchemo,
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
# cc=coord_cartesian(ylim=c(1e6,1e12)) #fig 2B has many logs padding (decieves to be in SS)
cc=coord_cartesian(ylim=c(1e-2,1e4))
gh=geom_hline(yintercept=0.63)
D%>%ggplot(aes(x=time,y=N,col=Tc))+geom_line(size=1)+gx+gy+tc(14)+sy+ltp+cc+gh
ggsave("~/Results/myelo/zhugeNchemoGeventsfig2B.png",height=6,width=6.5)
