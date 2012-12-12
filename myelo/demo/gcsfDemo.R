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

# FYI, the zhuge12 help page example renders Figure 6B fairly closely 



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
