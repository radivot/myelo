rm(list=ls(all=TRUE))
graphics.off()
setwd("/Users/radivot/case/grants/sachs/R/epo")
library(XLConnect)
wb <- loadWorkbook("EpoData.xls")
(d<- readWorksheet(wb, sheet = "Epo Endocytosis Triplicates"))
tpts=unique(d$time)
sig=c(summary(lm(y1~as.factor(time),data=d))$coef[2,2],
        summary(lm(y2~as.factor(time),data=d))$coef[2,2],
        summary(lm(y3~as.factor(time),data=d))$coef[2,2])
Sig=rep(1,24)%*%t(sig)
epo=d[,1:4]
names(epo)<-c("time","epo.e","epo.m","epo.i")
epo
save(epo,file="epo.rda")

library(deSolve)
dyn.load("epo.dll")
#test DLL
log10(c(2100,164,516))

thetahat=c(Bmax = 2.809384, Epo0 = 3.338145,Kd= 2.583673, kde  = -1.871801,
	kdi=-2.738496, ke= -1.195154, kex=-2.557234,kon=-4.066954,kt=-1.701394,scale= 2.191042)    
p0=10^thetahat
np=as.numeric(p0) # strip names
yout <- ode(y = c(Epo=np[2],EpoR=np[1],EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), times = tpts, 
		func = "derivsc", parms = p0, 
		dllname = "epo",initfunc = "parmsc",
		nout = 3, outnames = c("y1", "y2","y3")) 
yout  # checks out OK, so code in epo.c must be OK

fopt<-function(p,d)
{	p0=10^p
	np=as.numeric(p0) # strip names
	yout <- ode(y = c(Epo=np[2],EpoR=np[1],EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), times = tpts, 
				          func = "derivsc", parms = p0, dllname = "epo",initfunc = "parmsc",
				          nout = 3, outnames = c("y1", "y2","y3")) 
				  E=cbind(rep(yout[,8],each=3),rep(yout[,9],each=3),rep(yout[,10],each=3))
				  res=(E-d[,2:4])/Sig
				  sum(res^2)
}

strt=c(Bmax=3,Epo0=3,Kd=3,kde=-3,kdi=-3,ke=-3,kex=-3,kon=-3,kt=-3,scale=3)
strt=c(Bmax=2,Epo0=3,Kd=2,kde=-2,kdi=-3,ke=-1,kex=-2,kon=-4,kt=-2,scale=2)
strt=thetahat
#ssolm=optim(strt,fopt,method="L-BFGS-B",control=list(maxit=400),d=d,hessian=TRUE)
ssolm=optim(strt,fopt,control=list(maxit=400),d=d,hessian=TRUE)
(sig=sqrt(diag(solve(ssolm$hessian))))
ssolm$par
cbind(strt,thetahat,point=ssolm$par,lower=ssolm$par-1.96*sig,upper=ssolm$par+1.96*sig)

# mle2 is an optim wrapper in bbmle that has some conveniences 
require(bbmle) # note that dataframes used globally below are already centered in age
nLL<-function(Bmax,Epo0,Kd,kde,kdi,ke,kex,kon,kt,scale,x)  	
			{	p0=10^c(Bmax,Epo0,Kd,kde,kdi,ke,kex,kon,kt,scale)
				np=as.numeric(p0) # strip names
				yout <- ode(y = c(Epo=np[2],EpoR=np[1],EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), times = tpts, 
						func = "derivsc", parms = p0, dllname = "epo",initfunc = "parmsc",
						nout = 3, outnames = c("y1", "y2","y3")) 
				E=cbind(rep(yout[,8],each=3),rep(yout[,9],each=3),rep(yout[,10],each=3))
				res=(E-x[,2:4])/Sig
#				res[1:3,2]=0 # take out of sse since ic = zero => no headway
				sum(res^2)+(np[2]-2100)^2/44100  # 44100=(210^2)
				}

#fit0m <- mle2(nLL,start=list(Bmax=3,Epo0=3,Kd=3,kde=-3,kdi=-3,ke=-3,kex=-3,kon=-3,kt=-3,scale=3),
fit0m <- mle2(nLL,start=as.list(strt),
		method="L-BFGS-B",
		data=list(x=d))
summary(fit0m)
#cbind(coef(fit0m),confint(fit0m))
#system.time(p2<-profile(fit0m,which=2,alpha=0.05,maxsteps=20))
#lowLim=c(Bmax=2.7,Epo0=3,Kd=1.5,kde=-2,kdi=-3.2,ke=-1.23,kex=-2.8,kon=-4.3,kt=-1.85,scale=2.1)
#upLim=c(Bmax=3,Epo0=3.5,Kd=3,kde=-1.8,kdi=-2.4,ke=-1.14,kex=-2.2,kon=-3.95,kt=-1.65,scale=2.35)
#system.time(p2<-profile(fit0m,alpha=0.05,maxsteps=20,prof.lower=lowLim,prof.upper=upLim	)) # takes 936s = ~15 minutes
system.time(p2<-profile(fit0m,alpha=0.05,maxsteps=20)) # 920 secs
plot(p2)
windows()
plot(p2,absVal=FALSE)
confint(p2)


