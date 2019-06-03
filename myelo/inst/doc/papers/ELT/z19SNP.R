############################ ZHUGE12 SN + CHEMO only  ##################
# This code checks out fine with Figure 2b and code without events in gcsfDemo.R
library(tidyverse)
library(deSolve)
library(myelo)  
zhugePars19

zhuge19<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    if (Time < 0) {
      # etaN=etaNbar/(1+(N/the1script)^nu1)
      # etaP=etaPbar/(1+(P/the4script)^nu4)
      # An=Anmax*exp(-etaN*tauNM)
      # Ap=Apmax*exp(-etaP*tauPM)
      # beta=k0/(1+(S/the2)^s2)
      # Kn=f0/(1+(N/the1)^s1)
      # Kp=Kpbar/(1+Kp*P^s4)
      # dS=-(beta+Kn+Kp+kR)*S + Aq*beta*Sss
      # dN=-gamN*N + An*Kn*Sss
      # dP=-gamP*P + Ap*(1-exp(-gamP*tauPS))*Kp*Sss 
      # dP=-gamP*P + Ap*(Kp*Sss-exp(-gamP*tauPS)*Kp*Sss)
      dS=0
      dN=0
      dP=0
      # dP=Ap*(Kp*Sss-exp(-gamP*tauPS)*Kp*Sss)
    }	else {
      St=lagvalue(Time - tauS,1)
      Snt=lagvalue(Time - tauNM,1)
      Spt=lagvalue(Time - tauPM,1)
      Sptsum=lagvalue(Time - (tauPM+tauPS),1)
      Nt=lagvalue(Time - tauNM,2)
      Pt=lagvalue(Time - tauPM,3)
      Ptsum=lagvalue(Time - (tauPM+tauPS),3)
      etaNt=etaNbar/(1+(Nt/the1script)^nu1)
      etaPt=etaPbar/(1+(Pt/the4script)^nu4)
      etaPtsum=etaPbar/(1+(Ptsum/the4script)^nu4)
      As=2*exp(-gamS*tauS)
      Ant=Anmax*exp(-etaNt*tauNM)
      Apt=Apmax*exp(-etaPt*tauPM)
      Aptsum=Apmax*exp(-etaPtsum*tauPM)
      Kn=f0/(1+(N/the1)^s1)
      Knt=f0/(1+(Nt/the1)^s1)
      Kp=Kpbar/(1+Kp*P^s4)
      Kpt=Kpbar/(1+Kp*Pt^s4)
      Kptsum=Kpbar/(1+Kp*Ptsum^s4)
      beta=k0/(1+(S/the2)^s2)
      betat=k0/(1+(St/the2)^s2)
      dS=-(beta+Kn+Kp+kR)*S + As*betat*St
      dN=-gamN*N + Ant*Knt*Snt
      # dP=-gamP*P + Apt*(Kpt*Spt-exp(-gamP*tauPS)*Kptsum*Sptsum) 
      dP=-gamP*P + Apt*Kpt*Spt - exp(-gamP*tauPS)*Aptsum*Kptsum*Sptsum 
    }
    list(c(dS,dN,dP))
  })
}

(SS<-zhugePars19[c("Sss","Nss","Pss")])/1e6
times <- seq(-(zhugePars19["tauPM"]+zhugePars19["tauPS"]),5000,by=0.1)
# Sss     Nss     Pss 
# 1.1   690.0 31071.0 

yout <- dede(c(S=zhugePars19[["Sss"]],N=zhugePars19[["Nss"]],P=zhugePars19[["Pss"]]),
		times = times, func = zhuge19,	parms = zhugePars19)
plot(yout,select="S")
plot(yout,select="N")
plot(yout,select="P")
x0=as.vector(tail(yout,1))[2:4]
nms=colnames(yout)
names(x0)=nms[2:4]
x0/1e6
# S            N            P 
# 1.288574   454.968701 83814.645360 
newSS=x0
zhugePars19[["Sss"]]=x0[1]
zhugePars19[["Nss"]]=x0[2]
zhugePars19[["Pss"]]=x0[3]
times <- seq(-(zhugePars19["tauPM"]+zhugePars19["tauPS"]),500,by=0.1)
yout <- dede(c(S=zhugePars19[["Sss"]],N=zhugePars19[["Nss"]],P=zhugePars19[["Pss"]]),
             times = times, func = zhuge19,	parms = zhugePars19)
plot(yout)
plot(yout,select="S")
plot(yout,select="N")
plot(yout,select="P")
# steady states all look good

# kick N as before
(eventdat <- data.frame(var = c("N"),
                        time = c(25) ,
                        value = c(2e9),
                        method = c("rep")))

times <- seq(-(zhugePars19["tauPM"]+zhugePars19["tauPS"]),100,by=0.1)

yout <- dede(c(S=zhugePars19[["Sss"]],N=zhugePars19[["Nss"]],P=zhugePars19[["Pss"]]),
             times = times, func = zhuge19,	 parms = zhugePars19,
             events=list(data=eventdat),method="lsodar")
D=data.frame(yout)
tc=function(sz) theme_classic(base_size=sz)
gy=ylab("Neutrophil Counts")
gx=xlab("Days")
D%>%ggplot(aes(x=time,y=N))+geom_line(size=1)+gx+gy+tc(14) 
ggsave("~/Results/myelo/z12snpJumpN.png",height=6,width=6.5)
# 4 cycles in 25 days => now 6 day cycles instead of 21 day cycles

# 0.63e8*70/5# 882/uL as severly low (500 was my guess, close enough)
zhugePars19["Kpbar"] # 0.372 normal
zhugePars19["tauPM"] # 7 days is normal



zhugePars19["Kpbar"]=0.0744 # pattern 1a
zhugePars19["tauPM"]=13      # pattern 1a
times <- seq(-(zhugePars19["tauPM"]+zhugePars19["tauPS"]),300,by=0.1)
yout <- dede(c(S=zhugePars19[["Sss"]],N=zhugePars19[["Nss"]],P=zhugePars19[["Pss"]]),
             times = times, func = zhuge19,	parms = zhugePars19)
D=data.frame(yout)
d=D%>%gather(key="Cell",value="Counts",-time)%>%mutate(Cell=factor(Cell,levels=c("S","N","P")))
d%>%ggplot(aes(x=time,y=Counts))+facet_grid(Cell~.,scales = "free")+geom_line(size=1)+gx+tc(14) 
ggsave("~/Results/myelo/patt1a.png",height=6,width=6.5)


#try again with old SS off target to creat an impulse response

zhugePars19[c("Sss","Nss","Pss")]=SS
times <- seq(-(zhugePars19["tauPM"]+zhugePars19["tauPS"]),300,by=0.1)
yout <- dede(c(S=zhugePars19[["Sss"]],N=zhugePars19[["Nss"]],P=zhugePars19[["Pss"]]),
             times = times, func = zhuge19,	parms = zhugePars19)
D=data.frame(yout)
d=D%>%gather(key="Cell",value="Counts",-time)%>%mutate(Cell=factor(Cell,levels=c("S","N","P")))
d%>%ggplot(aes(x=time,y=Counts))+facet_grid(Cell~.,scales = "free")+geom_line(size=1)+gx+tc(14) 
###### ringing of transients is not sustained!!!!
ggsave("~/Results/myelo/patt1a.png",height=6,width=6.5)




zhugePars19[c("Sss","Nss","Pss")]=newSS
zhugePars19["Kpbar"]=2      # pattern 1b
zhugePars19["tauPM"]=14      # pattern 1b
times <- seq(-(zhugePars19["tauPM"]+zhugePars19["tauPS"]),300,by=0.1)
yout <- dede(c(S=zhugePars19[["Sss"]],N=zhugePars19[["Nss"]],P=zhugePars19[["Pss"]]),
             times = times, func = zhuge19,	parms = zhugePars19)
D=data.frame(yout)
d=D%>%gather(key="Cell",value="Counts",-time)%>%mutate(Cell=factor(Cell,levels=c("S","N","P")))
d%>%ggplot(aes(x=time,y=Counts))+facet_grid(Cell~.,scales = "free")+geom_line(size=1)+gx+tc(14) 
ggsave("~/Results/myelo/patt1b.png",height=6,width=6.5)


zhugePars19[c("Sss","Nss","Pss")]=SS
yout <- dede(c(S=zhugePars19[["Sss"]],N=zhugePars19[["Nss"]],P=zhugePars19[["Pss"]]),
             times = times, func = zhuge19,	parms = zhugePars19)
D=data.frame(yout)
d=D%>%gather(key="Cell",value="Counts",-time)%>%mutate(Cell=factor(Cell,levels=c("S","N","P")))
d%>%ggplot(aes(x=time,y=Counts))+facet_grid(Cell~.,scales = "free")+geom_line(size=1)+gx+tc(14) 

zhugePars19["Kpbar"]=9.5442 # pattern 2
zhugePars19["tauPM"]=11.583  # pattern 2
yout <- dede(c(S=zhugePars19[["Sss"]],N=zhugePars19[["Nss"]],P=zhugePars19[["Pss"]]),
             times = times, func = zhuge19,	parms = zhugePars19)
D=data.frame(yout)
d=D%>%gather(key="Cell",value="Counts",-time)%>%mutate(Cell=factor(Cell,levels=c("S","N","P")))
d%>%ggplot(aes(x=time,y=Counts))+facet_grid(Cell~.,scales = "free")+geom_line(size=1)+gx+tc(14) 
ggsave("~/Results/myelo/patt2.png",height=6,width=6.5)
tail(D)
