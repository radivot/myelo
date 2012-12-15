# Zhuge12
library(myelo)
zhugePars
# baseline 1-D model without any perturbations 
zhuge12N<-function(Time, State, Pars) {  # model with stem cells Q treated as constant
	with(as.list(c(State, Pars)), {
				An=exp(etaNP*tauNP-gam0*tauNM)
				if (Time < 0) 
					dN=-gamN*N + An*f0/(1+(Nss/the1)^s1)*Qss
				else
					dN=-gamN*N + An*f0/(1+(lagvalue(Time - tauN)/the1)^s1)*Qss
				list(c(dN))
			})
}
times <- c(-zhugePars[["tauN"]]:600)
yout <- dede(c(N=zhugePars[["Nss"]]), times = times, func = zhuge12N,	parms = zhugePars)
plot(yout)

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
	windows(width=6,height=6)
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




############################  BROOKS12 ######################
#### validate steady state solution with G=X=C=0 using a simplified version of the model
brooks12ss<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				Aq=2*exp(-gamS*tauS)
				An=exp(etaNP*tauNP-gam0*tauNMmax)
				if (Time < 0) {
					dQ=-(k0/(1+(Q/the2)^s2) + f0/(1+(N/the1)^s1)+kdel)*Q + Aq*k0/(1+(Qss/the2)^s2)*Qss
					dN=-gamN*N + An*f0/(1+(Nss/the1)^s1)*Qss
				}	else {
					Qts=lagvalue(Time - tauS)[1]
					Qtn=lagvalue(Time - tauN)[1]
					Ntn=lagvalue(Time - tauN)[2]
					dQ=-(k0/(1+(Q/the2)^s2) + f0/(1+(N/the1)^s1)+kdel)*Q + Aq*k0/(1+(Qts/the2)^s2)*Qts
					dN=-gamN*N + An*f0/(1+(Ntn/the1)^s1)*Qtn
				}
				list(c(dQ,dN))
			})
}

library(myelo)
brooksPars
times <- seq(-20,5000,by=1)
yout <- dede(c(Q=brooksPars[["Qss"]],N=brooksPars[["Nss"]]),
		times = times, func = brooks12ss,	parms = brooksPars)
tail(yout)
#### Yes, published steady state matches the model

# Now confirm for the full blown implementation of the model, i.e. brooks12
times <- seq(-20,200,by=.5)
brooksPars["hS"]=0  
brooksPars["hNP"]=0 

yout <- dede(c(Q=brooksPars[["Qss"]],N=brooksPars[["Nss"]],EtaNP=0,Gam0=0,GamS=0,G=0,X=0,C=0),
		times = times, func = brooks12,	parms = brooksPars)
tail(yout[,1:3])   # validates steady state 

# The following code chunk is also on the brook12 help page
brooksPars["hS"]=0.0702 
brooksPars["hNP"]=0.4275  
times <- c(-20:-1,seq(-.1,1,by=0.02),2:14,seq(14.01,15,by=0.02),16:100)
yout <- dede(c(Q=brooksPars[["Qss"]],N=brooksPars[["Nss"]],EtaNP=0,Gam0=0,GamS=0,G=0,X=0,C=0),
		times = times, func = brooks12,	parms = brooksPars)
# this should make Figure 3A. It has some similarity in the envelope time constant
par(mfrow=c(1,1),mar=c(4.7,6,3.3,1),oma=c(0,0,0,0),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
plot(times,yout[,3]/1e8,type="l",xlab="days",ylab="Neutrophils")

plot(yout) # This shows what a lot of other variables are doing


# Now take an initial stab at Figure 5
brooksPars["cycles"]=3  # three cycles
times <- c(-20:-1,seq(-.1,1,by=0.02),2:14,seq(14.01,15,by=0.02),16:300)
system.time({
brooksPars["kdel"]=0.0134 # default value used in Figure 5A
			youtA <- dede(c(Q=brooksPars[["Qss"]],N=brooksPars[["Nss"]],EtaNP=0,Gam0=0,GamS=0,G=0,X=0,C=0),
					times = times, func = brooks12,	parms = brooksPars)
			brooksPars["kdel"]=0.145 # value used in Figure 5B
			youtB <- dede(c(Q=brooksPars[["Qss"]],N=brooksPars[["Nss"]],EtaNP=0,Gam0=0,GamS=0,G=0,X=0,C=0),
					times = times, func = brooks12,	parms = brooksPars)
		}
		)
#Comment: We saved ~0.5 seconds, 10.8 -> 10.3, by changing lagvalue(Time - tauNv)[3] to lagvalue(Time - tauNv,3)
#in brooks12, and saved ~0.8 seconds by making only 3 lagvalue calls by using this block upfront.
#lagM<-lagvalue(Time - tauNMv)
#lagN<-lagvalue(Time - tauNv)
#lagS<-lagvalue(Time - tauS)  
# similar gains are obtained with fewer sample times via e.g. by=0.02 -> by=0.10
		
windows(width=8,height=4)
par(mfrow=c(1,2),mar=c(4.7,6,3.3,1),oma=c(0,0,0,0),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
plot(times,youtA[,3]/1e8,type="l",xlab="days",ylab="Neutrophils")
plot(times,youtB[,3]/1e8,type="l",xlab="days",ylab="Neutrophils") # ringing heavier but damped

# The following forcing function version of brooks12 gives 
# the user more work setting up interpolation functions & is slower, so this was 
# not promoted to be the brooks12() default.   

brooks12f<-function(Time, State, Pars,I0c,I0g) {
	with(as.list(c(State, Pars)), {
				gamSChemo=gamS+hS*C
				etaNPChemo=etaNP-hNP*C
				gamSv=gamMinS+(gamSChemo-gamMinS)*bS/(bS+G)
				gam0v=gamMin0+(gam0-gamMin0)*bn/(bn+G)
				etaNPv=etaNPChemo+(etaMaxNP-etaNPChemo)*G/(cn+G)
				tauNMv=tauNMmax/(1+(Vmax-1)*G/(bv+G))
				tauNv=tauNMv+tauNP
				if (Time < 0) {
					delEtaNP=etaNP*tauNP
					delGam0=gam0*tauNMmax
					delGamS=gamS*tauS
					An=exp(delEtaNP - delGam0)
					Aq=2*exp(-delGamS)
					dQ=-(k0/(1+(Q/the2)^s2) + f0/(1+(N/the1)^s1)+kdel)*Q + Aq*k0/(1+(Qss/the2)^s2)*Qss
					dN=-gamN*N + An*f0/(1+(Nss/the1)^s1)*Qss
					dC=0;	dX=0;	dG=0; trt=FALSE
				}	else {
					lagM<-lagvalue(Time - tauNMv)
					lagN<-lagvalue(Time - tauNv)
					lagS<-lagvalue(Time - tauS)
					delEtaNP=lagM[3]-lagN[3]
					delGam0=Gam0 -lagM[4]
					delGamS=GamS -lagS[5]
					An=exp(delEtaNP - delGam0)
					Aq=2*exp(-delGamS)
					Qts=lagS[1]
					Qtn=lagN[1]
					Ntn=lagN[2]
					dQ=-(k0/(1+(Q/the2)^s2) + f0/(1+(N/the1)^s1)+kdel)*Q + Aq*k0/(1+(Qts/the2)^s2)*Qts
					dN=-gamN*N + An*f0/(1+(Ntn/the1)^s1)*Qtn
					dC= I0c(Time)/phi - del*C
					dX= I0g(Time) + kT*VB*G -kB*X  
					dG= Gprod - kT*G + kB*X/VB -gamG*G - sig*N*G^2/(kG+G^2)
				}
				dEtaNP=etaNPv
				dGam0=gam0v
				dGamS=gamSv
				list(c(dQ,dN,dEtaNP,dGam0,dGamS,dG,dX,dC),
						c(delEtaNP=delEtaNP,delGam0=delGam0,delGamS=delGamS,
						etaNPv=etaNPv,gam0v=gam0v,gamSv=gamSv,tauNMv=tauNMv,tauNv=tauNv,An=An,Aq=Aq))
			})
}


ftimes=seq(-10,300,0.1)
gtab=ctab= as.data.frame(list(times = ftimes,
				import = rep(0, length(ftimes))))
ctab$import[(ctab$times>= 0)&(ctab$times%/%brooksPars["T"]<brooksPars["cycles"])&
				(ctab$times%%brooksPars["T"] < brooksPars["Delc"])] <- brooksPars["Dc"]/brooksPars["Delc"]
# the following isn't needed here since G-CSF=0 but would be needed in general
gtab$import[(gtab$times>= 0)&(gtab$times%/%brooksPars["T"]<brooksPars["cycles"])&
				(gtab$times%%brooksPars["T"] > brooksPars["T1"])&
				(gtab$times%%brooksPars["T"] < brooksPars["T1"]+brooksPars["Delg"])] <- brooksPars["Dg"]/brooksPars["Delg"]

cimp <- approxfun(ctab$times, ctab$import, rule = 2)
gimp <- approxfun(gtab$times, gtab$import, rule = 2)
windows(width=10,height=4)
par(mfrow=c(1,3),mar=c(4.7,6,3.3,1),oma=c(0,0,0,0),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
plot(ftimes,cimp(ftimes),type="l")
plot(ftimes,cimp(ftimes),type="l",xlim=c(-2,23))
plot(ftimes,gimp(ftimes),type="l",xlim=c(-2,23))

system.time({
			brooksPars["kdel"]=0.0134 # default value used in Figure 5A
			youtA <- dede(c(Q=brooksPars[["Qss"]],N=brooksPars[["Nss"]],EtaNP=0,Gam0=0,GamS=0,G=0,X=0,C=0),
					times = times, func = brooks12f,	parms = brooksPars,I0c=cimp,I0g=gimp)
			brooksPars["kdel"]=0.145 # value used in Figure 5B
			youtB <- dede(c(Q=brooksPars[["Qss"]],N=brooksPars[["Nss"]],EtaNP=0,Gam0=0,GamS=0,G=0,X=0,C=0),
					times = times, func = brooks12f,	parms = brooksPars,I0c=cimp,I0g=gimp)
		}
)
# slower at 11.5 secs. The same plots are produced 
windows(width=8,height=4)
par(mfrow=c(1,2),mar=c(4.7,6,3.3,1),oma=c(0,0,0,0),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
plot(times,youtA[,3]/1e8,type="l",xlab="days",ylab="Neutrophils")
plot(times,youtB[,3]/1e8,type="l",xlab="days",ylab="Neutrophils") # ringing heavier but damped




########## !!!!!!!!!!! SCHOLZ ET AL BELOW NEEDS FIXING SO IGNORE FOR NOW ########
# The following code pertains to the G-CSF model of M Scholz TBMM 2012
library(myelo)
# first we show that the Z function form they use is more complicated than need be
Z<-function(x,Amin,Amax,Anor,b){
	if((Anor<Amax)&(Anor>Amin))
		res=Amax-(Amax-Amin)*exp(-log((Amax-Amin)/(Amax-Anor))*x^b) else res=Anor
	res }  # since this 
#equivalent form shows better that x=1=>Z=Anor,  x=inf=>Z=Amax, and  x=0=>Z=Amin
Z2<-function(x,Amin,Amax,Anor,b) res = Amax - (Amax-Amin)*((Amax-Anor)/(Amax-Amin))^(x^b)
(x=10^c(-12:5)) #length(scholzPars) # 128 parameters
with(as.list(scholzPars),{ #prove equal by plotting both: y as points, y2 as lines
			y=Z(x,AminPGBF,AmaxPGBF,AnorPGBF,AbPGBF)
			y2=Z2(x,AminPGBF,AmaxPGBF,AnorPGBF,AbPGBF)
			plot(x,y,log="xy",ylab="Z",xlab="G-CSF",ylim=c(.01,1e3))
			lines(x,y2,col="green")})
# now run a simulation to show that RHS code works at least syntactically
X0=rep(1,24)
names(X0)<-c("g","g1","g2","g3","g4","S","CG","PGB","G4a","G4b","G4c","G4d","G4e",
		"G5a","G5b","G5c","G5d","G5e","G6a","G6b","G6c","G6d","G6e","GRA")
# WARNING: the following is needed because the model is not completely specified! 
APGBout=APGBin=ACGout=ACGin=GSS=G6nor=GRAnor=1
out   <- ode(X0, c(0,10^(1:4)), scholz12, scholzPars)
tail(out) # it hit steady state
(X0=out[dim(out)[1],2:25])
(eventdat <- data.frame(var = "g", time = 0, value = 2, method = "mult"))
out   <- ode(X0, -20:600, scholz12, scholzPars,events = list(data = eventdat))
plot(out)
