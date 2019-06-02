############################ ZHUGE12  ##################
library(tidyverse)
library(deSolve)
library(myelo)  
zhugePars #matches Table 1
#        Qss       gamS    gamMinS    gamMaxS       tauS         k0       the2         s2        Nss 
# 1.1000e+06 7.0000e-02 3.0000e-02 2.0000e-01 2.8000e+00 8.0000e+00 3.0000e+05 4.0000e+00 6.3000e+08 
#       gamN      tauNP      tauNM  tauNMgcsf       tauN      etaNP   etaMinNP   etaMaxNP       gam0 
# 2.4000e+00 5.0000e+00 6.0000e+00 2.0000e+00 1.1000e+01 2.5420e+00 2.0420e+00 3.0552e+00 2.7000e-01 
#    gamMin0         f0       the1         s1       kdel          T         T1 
# 1.2000e-01 4.0000e-01 3.6000e+07 1.0000e+00 1.0000e-02 2.1000e+01 4.0000e+00 
#HSC constant
# attach(as.list(zhugePars))
# An=exp(etaNP*tauNP-gam0*tauNM)
# -gamN*Nss + An*f0/(1+(Nss/the1)^s1)*Qss
zhuge12N<-function(Time, State, Pars) {  
  with(as.list(c(State, Pars)), {
    An=exp(etaNP*tauNP-gam0*tauNM)
    if (Time < 0) 
      dN=-gamN*N + An*f0/(1+(Nss/the1)^s1)*Qss
    else
      dN=-gamN*N + An*f0/(1+(lagvalue(Time - tauN)/the1)^s1)*Qss
    list(c(dN))
  })
}

# times= -zhugePars[["tauN"]]:100
times= seq(-zhugePars[["tauN"]],200,by=0.1)

yout <- dede(c(N=zhugePars[["Nss"]]), times = times, func = zhuge12N,	parms = zhugePars)
D=data.frame(yout)
tc=function(sz) theme_classic(base_size=sz)
gy=ylab("Neutrophil Counts")
gx=xlab("Days")
D%>%ggplot(aes(x=time,y=N))+geom_line(size=1)+gx+gy+tc(14) #### take to steady state

tail(D,1)["N"]  # Nss=639805537 (run out to 5000 days)
zhugePars["Nss"]=639805537
# 
# derivs <- function(t, var, parms)   list(dvar = rep(0, 2))
# yini <- c(v1 = 1, v2 = 2)
# times <- seq(0, 10, by = 0.1)
# (eventdat <- data.frame(var = c("v1", "v2", "v2", "v1"),
#                        time = c(1, 1, 5, 9) ,
#                        value = c(1, 2, 3, 4),
#                        method = c("add", "mult", "rep", "add")))
# out <- vode(func = derivs, y = yini, times = times, parms = NULL, events = list(data = eventdat))
# plot(out)
(eventdat <- data.frame(var = c("N"),
                       time = c(25) ,
                       # value = c(10),
                       value = c(2e9),
                       method = c("rep")))
                       # method = c("mult")))
times= seq(-zhugePars[["tauN"]],100,by=0.01)
yout=dede(c(N=zhugePars[["Nss"]]),times=times,func=zhuge12N,
          parms=zhugePars,events=list(data=eventdat),method="lsodar")
D=data.frame(yout)
D%>%ggplot(aes(x=time,y=N))+geom_line(size=1)+gx+gy+tc(14) #### take to steady state




# now add just chemo to this model to recreate Fig 2B
# Chemo acts only on proliferation,  so add one state to integrate eta 
zhuge12Nchemo<-function(Time, State, Pars) {  # model with stem cells Q treated as constant
  with(as.list(c(State, Pars)), {
    dEta=etaNP
    if (Time < 0) {
      An=exp(etaNP*tauNP-gam0*tauNM)  # no gcsf or chemo perturbations for negative times
      dN=-gamN*N + An*f0/(1+(Nss/the1)^s1)*Qss
    }
    else{# quotients and remaninders: 21.4%/%10  21.4%%10
      if (Time%%T < 1) dEta=etaMinNP  # in chemo period
      delEta=lagvalue(Time - tauNM)[2]-lagvalue(Time - tauN)[2]
      An=exp(delEta - gam0*tauNM)
      dN=-gamN*N + An*f0/(1+(lagvalue(Time - tauN)[1]/the1)^s1)*Qss
    }
    list(c(dN,dEta))
  })
}
times <- seq(-zhugePars[["tauN"]],200,by=0.1)
zhugePars["T"]=18
yout18 <- dede(c(N=zhugePars[["Nss"]],Eta=0), times = times, func = zhuge12Nchemo,	parms = zhugePars)
plot(yout18)
zhugePars["T"]=23
yout23 <- dede(c(N=zhugePars[["Nss"]],Eta=0), times = times, func = zhuge12Nchemo,	parms = zhugePars)
plot(yout23)

# the following block makes Figure 2B 
myplot=function(times, y1,y2) {
  graphics.off()
  quartz(width=6,height=6)
  par(mfrow=c(1,1),mar=c(4.7,6,3.3,1),oma=c(0,0,0,0),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
  plot(times,y1/1e8,type="l",lty=2,log="y",yaxt="n",ylim=c(1e-2,1e4),ylab="",xlab="days")
  lines(times,y2/1e8,lty=1)
  abline(h=0.63)
  axis(side=2,las=1, at=c(1e-2,1,1e2,1e4),labels=expression(10^-2,10^0,10^2,10^4))
  mtext(expression(paste("Neutrophils in ",10^8," per kg")),side=2,line=3.5,cex=2)
}
myplot(yout23[,1],yout23[,2],yout18[,2])




# Holding the chemo treatment cycle T fixed at 23 days (at resonance), now change the time to 
# GCSF in days after chemo last began
zhuge12NgcsfChemo<-function(Time, State, Pars) {  # model with stem cells Q treated as constant
  with(as.list(c(State, Pars)), {
    dEta=etaNP
    dGam=gam0
    if (Time < 0) {
      An=exp(etaNP*tauNP-gam0*tauNM)  # no gcsf or chemo perturbations for negative times
      dN=-gamN*N + An*f0/(1+(Nss/the1)^s1)*Qss
    }
    else{# quotients and remaninders: 21.4%/%10  21.4%%10
      if (Time%%T < 1) dEta=etaMinNP # in chemo period
      if ( (Time%%T -T1 > 0)&(Time%%T -T1 < 1) ) { # in G-CSF exposure period
        dEta=etaMaxNP
        dGam=gamMin0
      }
      delEta=lagvalue(Time - tauNM)[2]-lagvalue(Time - tauN)[2]
      delGam=Gam -lagvalue(Time - tauNM)[3]
      An=exp(delEta - delGam)
      dN=-gamN*N + An*f0/(1+(lagvalue(Time - tauN)[1]/the1)^s1)*Qss
    }
    list(c(dN,dEta,dGam))
  })
}

zhugePars["T1"]=1
zhugePars["T"]=23
yout1 <- dede(c(N=zhugePars[["Nss"]],Eta=0,Gam=0), times = times, func = zhuge12NgcsfChemo,	parms = zhugePars)
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

