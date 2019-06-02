############################ ZHUGE12 CHEMO+GCSF  ##################
# Holding the chemo treatment cycle T fixed at 23 days (at resonance), change the time to 
# GCSF in days after chemo last began. Compare T1 = 1 and 10 days

library(tidyverse)
library(deSolve)
library(myelo)  

zhugePars["T"]=23
zhugePars["T1"]=10  # later do 1 day
zhugePars["Nss"]=639805537  #SS found in readme.md
Tf=200

mkEventsCG=function(zhugePars,Tf) {
  (Tc1=seq(0,Tf,zhugePars["T"]))
  (Tg1=seq(zhugePars["T1"],Tf,zhugePars["T"])) #1 day after start of chemo
  # (Tg1=seq(zhugePars["T1"]+1,Tf,zhugePars["T"])) #1 day after end of chemo
  Nc1=length(Tc1) 
  Ng1=length(Tg1)
  (Tc2=seq(1,Tf,zhugePars["T"]))
  (Tg2=seq(zhugePars["T1"]+1,Tf,zhugePars["T"])) #1 day after start of chemo
  # (Tg2=seq(zhugePars["T1"]+2,Tf,zhugePars["T"]))  #1 day after end of chemo
  Nc2=length(Tc2) 
  Ng2=length(Tg2)
  (etacB <- data.frame(var = rep("eta",Nc1),
                           time = Tc1 ,
                           value = rep(zhugePars["etaMinNP"],Nc1),
                           method = rep("rep",Nc1)))
  (etacE <- data.frame(var = rep("eta",Nc2),
                           time = Tc2 ,
                           value = rep(zhugePars["etaNP"],Nc2),
                           method = rep("rep",Nc2)))
  (etagB <- data.frame(var = rep("eta",Ng1),
                      time = Tg1 ,
                      value = rep(zhugePars["etaMaxNP"],Ng1),
                      method = rep("rep",Ng1)))
  (etagE <- data.frame(var = rep("eta",Ng2),
                      time = Tg2 ,
                      value = rep(zhugePars["etaNP"],Ng2),
                      method = rep("rep",Ng2)))
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
  rbind(etacB,etacE,etagB,etagE,gamB,gamE,tauB,tauE)%>%arrange(time) #no warnings
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
    tauN=tauNP+tau #total tauN changes with maturation time changes
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

times= seq(-zhugePars[["tauN"]],Tf,by=0.01)
yout=dede(c(N=zhugePars[["Nss"]],eta=zhugePars[["etaNP"]], Eta=0,
                                 gam=zhugePars[["gam0"]], Gam=0, 
                                 tau=zhugePars[["tauNM"]] ),
          times=times,func=zhuge12NchemoG,
          parms=zhugePars,events=list(data=mkEventsCG(zhugePars,Tf)),method="lsodar")

D10=data.frame(yout)

zhugePars["T1"]=1
yout=dede(c(N=zhugePars[["Nss"]],eta=zhugePars[["etaNP"]], Eta=0,
            gam=zhugePars[["gam0"]], Gam=0, 
            tau=zhugePars[["tauNM"]] ),
          times=times,func=zhuge12NchemoG,
          parms=zhugePars,events=list(data=mkEventsCG(zhugePars,Tf)),method="lsodar")
D1=data.frame(yout)

D1$Tg="1 Day"
D10$Tg="10 Days"
D=bind_rows(D1,D10)%>%mutate(N=N/1e8)

tc=function(sz) theme_classic(base_size=sz)
gy=ylab("Neutrophil Counts")
gx=xlab("Days")
sy=scale_y_log10()
ltp=theme(legend.position="top")
# cc=coord_cartesian(ylim=c(1e6,1e12))
cc=coord_cartesian(ylim=c(1e-2,1e4))
gh=geom_hline(yintercept=0.63)
D%>%ggplot(aes(x=time,y=N,col=Tg))+geom_line(size=1)+gx+gy+tc(14)+sy+ltp+cc+gh
ggsave("~/Results/myelo/zhugeNchemoGeventsfig3B.png",height=6,width=6.5)
# ggsave("~/Results/myelo/zhugeNchemoGeventsfig3Bplus1.png",height=6,width=6.5)
# two prongs of peaks with T1=10 match paper, as do troughs (as without events [in pdf])
# T1=1 does not match paper, neither counting 1 from beginning or end (plus1) of chemo interval
