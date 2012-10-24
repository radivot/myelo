library(myelo)
############ Step 1: Show how much faster C coded RHS are relative to R coded RHS  ###############3
dyn.load(file.path(system.file(paste("libs",Sys.getenv("R_ARCH"),sep=""), package = "myelo"),"myelo.dll"))

thetahat=c(Bmax = 2.821, Epo0 = 3.322,Kd= 2.583, kde  = -1.884,
		kdi=-2.730, ke= -1.177, kex=-2.447,kon=-4.091,kt=-1.758,scale= 2.210)  # from Table 1 in Chao 2010  

(p0=10^thetahat) # in thetahat optimization, logs constrain p0 estimates to be positive. () echos result
tpts=c(0, 5, 20, 60, 120, 180, 240, 300) # the time points at which we have data, i.e. tpts=unique(epo$time)

system.time(for (i in 1:50) yout <-ode(y = c(Epo=2178,EpoR=645,EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), 
				times = tpts, 
				func = raue10, parms = p0))
with(epo,matplot(time,epo[,-1]))
matlines(yout[,1],yout[,-(1:7)])

system.time(for (i in 1:50) yout <-ode(y = c(Epo=2178,EpoR=645,EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), 
					times = tpts, 
					func = "derivsc", parms = p0, 
					dllname = "myelo",initfunc = "parmsc",
					nout = 3, outnames = c("y1", "y2","y3")) )
with(epo,matplot(time,epo[,-1]))
matlines(yout[,1],yout[,-(1:7)])
# so the C code is ~40 times faster, and it yields the same result. 

##############  Step 2: show how the optimum parameters used above were obtained.
#### From here on out we only use the C coded ODEs

# we first note that there were some prior point estimates 
(xcel=c(Epo0=2100,Kd=164,Bmax=516))  #in original Excel data file of Raue et al 2010. 
p0[names(xcel)] # that differ somewhat from the optimum values in thetahat found here
#Epo0        Kd      Bmax 
#2178.4370  383.4184  644.7391 


# in the nonlinear least squares fit objective functions defined below 
np=as.numeric(p0) # strips names before sending values in as the IC 
c(Epo=p0[2],EpoR=p0[1]) #to avoid this sort of thing
c(Epo=np[2],EpoR=np[1]) # and instead have this sort of thing

#Also, there three model outputs to fit to data, so residuals need to be 
# weighted by the reciprocols of the measurement error, obtained here using
# saturated one-way anova models of the time courses. 
sig=c(summary(lm(epo.e~as.factor(time),data=epo))$coef[2,2],
		summary(lm(epo.m~as.factor(time),data=epo))$coef[2,2],
		summary(lm(epo.i~as.factor(time),data=epo))$coef[2,2])
(Sig=rep(1,24)%*%t(sig))  # make the sd's match the dims of the data


# first let's see how this is done with optim
fopt<-function(p,d,W)
{	p0=10^p
	np=as.numeric(p0) # strip names
	yout <- ode(y = c(Epo=np[2],EpoR=np[1],EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), times = tpts, 
				          func = "derivsc", parms = p0, dllname = "myelo",initfunc = "parmsc",
				          nout = 3, outnames = c("y1", "y2","y3")) 
				  E=cbind(rep(yout[,8],each=3),rep(yout[,9],each=3),rep(yout[,10],each=3))
				  res=(E-d[,2:4])/Sig
				  sum(res^2)+(np[2]-2100)^2/44100  # only the Epo0 prior is used. Is a weight of 1 datapoint fair?
}

strt=c(Bmax=3,Epo0=3,Kd=3,kde=-3,kdi=-3,ke=-3,kex=-3,kon=-3,kt=-3,scale=3)
ssolm=optim(strt,fopt,control=list(maxit=400),d=epo,hessian=TRUE)
(sig=sqrt(diag(solve(ssolm$hessian))))
ssolm$par
cbind(strt,thetahat,point=ssolm$par,lower=ssolm$par-1.96*sig,upper=ssolm$par+1.96*sig)

myplot=function(p) {
		p0=10^p
		np=as.numeric(p0) # strip names
		yout <-ode(y = c(Epo=np[2],EpoR=np[1],EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), 
			times = tpts, func = "derivsc", parms = p0, 
			dllname = "myelo",initfunc = "parmsc",
			nout = 3, outnames = c("y1", "y2","y3")) 
	with(epo,matplot(time,epo[,-1]))
	matlines(yout[,1],yout[,-(1:7)])
}
myplot(ssolm$par) # maybe stuck in a local minimum. 

# mle2 is an optim wrapper in bbmle that has some conveniences 
require(bbmle) # note that dataframes used globally below are already centered in age
nLL<-function(Bmax,Epo0,Kd,kde,kdi,ke,kex,kon,kt,scale,x)  	
			{	p0=10^c(Bmax,Epo0,Kd,kde,kdi,ke,kex,kon,kt,scale)
				np=as.numeric(p0) # strip names for optimized IC parameters. Bmax is both a param and EpoR0
				yout <- ode(y = c(Epo=np[2],EpoR=np[1],EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), times = tpts, 
						func = "derivsc", parms = p0, dllname = "myelo",initfunc = "parmsc",
						nout = 3, outnames = c("y1", "y2","y3")) 
				E=cbind(rep(yout[,8],each=3),rep(yout[,9],each=3),rep(yout[,10],each=3))
				res=(E-x[,2:4])/Sig
#				res[1:3,2]=0 # take out of sse since ic = zero => no headway
				sum(res^2)+2*(np[2]-2100)^2/44100  # 44100=210^2  soft constraint = prior on Epo IC 
				}


fit0 <- mle2(nLL,start=as.list(strt),
		method="Nelder-Mead",  # same as optim default, so should be the same,
		data=list(x=epo))
(s=summary(fit0))   # but even with more trouble with CI, 
myplot(coef(s)[,1])  # the fit looks better


fit0 <- mle2(nLL,start=as.list(strt),
		method="L-BFGS-B", # see if this helps 
		data=list(x=epo))
(s=summary(fit0))   # still trouble with CIs 
myplot(coef(s)[,1])  # but fit is now even better


fit0 <- mle2(nLL,start=as.list(thetahat), # now start at Table 1 point estimate
		method="Nelder-Mead", # L-BFGS-B and SANNN yield similarly tight std devs 
		data=list(x=epo))
(s=summary(fit0))   # trouble with CIs is now gone, 
myplot(coef(s)[,1])  # and fit still looks great, but the std devs seem too tight given 
# all of the previous CI trouble. 



# this gives an error, and may not work with RHS hard coded in C
#fit0 <- mle2(nLL,start=as.list(thetahat[-3]), 
#		fix=list(Kd=2.583), 
#		method="Nelder-Mead", 
#		data=list(x=epo))
#(s=summary(fit0))   
#myplot(coef(s)[,1]) 

lowLim=c(Bmax=2.7,Epo0=3.22,Kd=1.4,kde=-1.95,kdi=-3.1,ke=-1.205,kex=-2.8,kon=-4.25,kt=-1.84,scale=2.13)
upLim=c(Bmax=2.95,Epo0=3.42,Kd=3,kde=-1.82,kdi=-2.5,ke=-1.145,kex=-2.2,kon=-3.96,kt=-1.66,scale=2.33)
fit0 <- mle2(nLL,start=as.list(thetahat), 
		optimizer="constrOptim", ui=rep(1,10),ci=lowLim,
		method="Nelder-Mead", 
		data=list(x=epo))
(s=summary(fit0))   
myplot(coef(s)[,1]) 


# with finite CI, we can now use the profile likelihood function on the output fit0
# set the limits of the scan to those seen in Figure 3 of Chaos 2010 
system.time(p2<-profile(fit0,alpha=0.05,maxsteps=20,prof.lower=lowLim,prof.upper=upLim	)) 
# returms object of class mle2 because it finds a better fit. 
fit0 <- mle2(nLL,start=as.list(coef(p2)), # start at output p2 of profile at its termination
		method="L-BFGS-B", 
		data=list(x=epo))
(s=summary(fit0))   # CIs still tight 
myplot(coef(s)[,1])  # and fit still looks good
# but this continues to find a better fit
system.time(p2<-profile(fit0,alpha=0.05,maxsteps=20,prof.lower=lowLim,prof.upper=upLim	)) 
# and it doesn't matter if fit0



# this takes ~900 seconds = 15 minutes = too long => use multicores
plot(p2)
windows()
plot(p2,absVal=FALSE)
confint(p2) # CI still seem too tight, and indeed, they are relative to Table 1
#        2.5 %    97.5 %        95%CI in Table1
#Bmax   2.742880  2.846018   +2.710  +2.932
#Epo0   3.277322  3.362898   +3.227  +3.400
#Kd     2.288379  2.809257   +1.641  +2.993
#kde   -1.902429 -1.847011   -1.941  -1.829
#kdi   -2.862994 -2.654893   -3.083  -2.535
#ke    -1.216857 -1.176136   -1.203  -1.150
#kex   -2.761644 -2.419782   -2.764  -2.225
#kon   -4.105063 -3.994924   -4.208  -3.973
#kt    -1.733934 -1.671456   -1.828  -1.683
#scale  2.166216  2.251927   +2.133  +2.305
# Wald CI = 2*sd from summary(fit0) are consistent with our profiles

system.time({
myprof=function(k) profile(fit0,which=k,alpha=0.05,prof.lower=lowLim,prof.upper=upLim,maxsteps=20)
L=lapply(1:10,myprof)
		})
# this crashed with "hit maximum number of steps (Kd)"
# after 880 secs with max steps=20 and after 1400 with max steps =100
