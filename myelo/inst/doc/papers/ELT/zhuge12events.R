############################ ZHUGE12 CHEMO+GCSF  ##################

# Holding the chemo treatment cycle T fixed at 23 days (at resonance), now change the time to 
# GCSF in days after chemo last began

library(tidyverse)
library(deSolve)
library(myelo)  

zhugePars["T"]=23
zhugePars["T1"]=10  # later do 1 day
zhugePars["Nss"]=639805537  #SS found in readme.md
Tf=200

mkEventsCG=function(zhugePars,Tf) {
  (Tc1=seq(0,Tf,zhugePars["T"]))
  (Tg1=seq(zhugePars["T1"],Tf,zhugePars["T"]))
  Nc1=length(Tc1) 
  Ng1=length(Tg1)
  (Tc2=seq(1,Tf,zhugePars["T"]))
  (Tg2=seq(zhugePars["T1"]+1,Tf,zhugePars["T"]))
  Nc2=length(Tc2) 
  Ng2=length(Tg2)
  (etaB <- data.frame(var = rep("eta",Nc1),
                           time = Tc1 ,
                           value = rep(zhugePars["etaMinNP"],Nc1),
                           method = rep("rep",Nc1)))
  (etaE <- data.frame(var = rep("eta",Nc2),
                           time = Tc2 ,
                           value = rep(zhugePars["etaNP"],Nc2),
                           method = rep("rep",Nc2)))
  (gamB <- data.frame(var = rep("gam",Ng1),
                      time = Tg1 ,
                      value = rep(zhugePars["gamMin0"],Ng1),
                      method = rep("rep",Ng1)))
  (gamE <- data.frame(var = rep("gam",Ng2),
                      time = Tg2 ,
                      value = rep(zhugePars["gam0"],Ng2),
                      method = rep("rep",Ng2)))
  (tauB <- data.frame(var = rep("tau",Ng1),
                      time = Tg1 ,
                      value = rep(zhugePars["tauNMgcsf"],Ng1),
                      method = rep("rep",Ng1)))
  (tauE <- data.frame(var = rep("tau",Ng2),
                      time = Tg2 ,
                      value = rep(zhugePars["tauNM"],Ng2),
                      method = rep("rep",Ng2)))
  # bind_rows(etaB,etaE,gamB,gamE,tauB,tauE)%>%arrange(time)
  rbind(etaB,etaE,gamB,gamE,tauB,tauE)%>%arrange(time) #no warnings
}

eventDF=mkEventsCG(zhugePars,Tf)
sapply(eventDF,class)
head(eventDF,10)


zhuge12NchemoG<-function(Time, State, Pars) {  # model with stem cells Q treated as constant
  with(as.list(c(State, Pars)), {
    deta=0    #this is etaNP as a state
    dEta=eta
    dgam=0    #this is gam0 as a state
    dGam=gam
    dtau=0  # this is tauNM as a state. Go from 6 days of maturation to 2 while on Gcsf
    tauN=tanNP+tau #total tauN changes with maturation time changes
    if (Time < 0) {
      An=exp(etaNP*tauNP-gam0*tauNM)  # no gcsf or chemo perturbations for negative times
      dN=-gamN*N + An*f0/(1+(Nss/the1)^s1)*Qss
    }
    else{# quotients and remaninders: 21.4%/%10  21.4%%10
      # if (Time%%T < 1) dEta=etaMinNP # in chemo period
      # if ( (Time%%T -T1 > 0)&(Time%%T -T1 < 1) ) { # in G-CSF exposure period
      #   dEta=etaMaxNP
      #   dGam=gamMin0
      # }
      delEta=lagvalue(Time - tau,3)-lagvalue(Time - tauN,3)
      delGam=Gam -lagvalue(Time - tau,5)
      An=exp(delEta - delGam)
      dN=-gamN*N + An*f0/(1+(lagvalue(Time - tauN)[1]/the1)^s1)*Qss
    }
    list(c(dN,deta,dEta,dgam,dGam,dtau))
  })
}

yout=dede(c(N=zhugePars[["Nss"]],eta=zhugePars[["etaNP"]], Eta=0,
                                 gam=zhugePars[["gam0"]], Gam=0,
            ),
          times=times,func=zhuge12NchemoG,
          parms=zhugePars,events=list(data=mkEventsCG(zhugePars,Tf)),method="lsodar")



zhugePars["T1"]=1
zhugePars["T"]=23
yout1 <- dede(c(N=zhugePars[["Nss"]],Eta=0,Gam=0), times = times, func = ,	parms = zhugePars)
plot(yout1)

zhugePars["T1"]=10
yout10 <- dede(c(N=zhugePars[["Nss"]],Eta=0,Gam=0), times = times, func = zhuge12NgcsfChemo,	parms = zhugePars)
plot(yout10)

# the following block makes Figure 3B 
myplot(yout10[,1],yout10[,2],yout1[,2])




##### Now add a stem cell state variable to the model
library(myelo)
zhugePars
times <- seq(-zhugePars[["tauN"]],300,by=0.1)
zhugePars["T"]=18
zhugePars["T1"]=NA
yout18 <- dede(c(Q=zhugePars[["Qss"]],N=zhugePars[["Nss"]],Eta=0,Gam=0,GamS=0), 
               times = times, func = zhuge12,	parms = zhugePars)
plot(yout18)

zhugePars["T"]=21
yout21 <- dede(c(Q=zhugePars[["Qss"]],N=zhugePars[["Nss"]],Eta=0,Gam=0,GamS=0),
               times = times, func = zhuge12,	parms = zhugePars)
plot(yout21)
# the following attempts to make Figure 4B, but T=21 destabilize ~100 days later 
myplot(yout21[,1],yout21[,3],yout18[,3])

# the zhuge12 help page renders Figure 6B fairly closely 
times <- seq(-zhugePars[["tauN"]],200,by=0.1)
zhugePars["T"]=21  #chemo coming in at natural freqency
zhugePars["T1"]=4  # G-CSF a one day dose 4 days after chemo
yout4 <- dede(c(Q=zhugePars[["Qss"]],N=zhugePars[["Nss"]],Eta=0,Gam=0,GamS=0),
              times = times, func = zhuge12,	parms = zhugePars)
yout4=data.frame(yout4)%>%mutate(T1="4 Days",N=N/1e8)%>%select(time,N,T1)
head(yout4)
zhugePars["T1"]=14 ##G-CSF a one day dose 14 days after chemo
yout14 <- dede(c(Q=zhugePars[["Qss"]],N=zhugePars[["Nss"]],Eta=0,Gam=0,GamS=0),
               times = times, func = zhuge12,	parms = zhugePars)
yout14=data.frame(yout14)%>%mutate(T1="14 Days",N=N/1e8)%>%select(time,N,T1)
D=bind_rows(yout4,yout14)
# D$T1=factor(D$T1)
gy=ylab(expression(paste("Neutrophils in ",10^8," per kg")))
sy=scale_y_log10()
ltp=theme(legend.position="top")
D%>%ggplot(aes(x=time,y=N,col=T1))+geom_line(size=1)+gx+gy+sy+tc(14)+ltp 
ggsave("~/Results/myelo/zhugeResonance.png",height=6,width=6.5)

