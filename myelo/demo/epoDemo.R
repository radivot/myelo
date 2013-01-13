library(myelo)
thetahat=c(Bmax = 2.821, Epo0 = 3.322,Kd= 2.583, kde  = -1.884,
		kdi=-2.730, ke= -1.177, kex=-2.447,kon=-4.091,kt=-1.758,scale= 2.210)  # from Table 1 in Chao 2010  
# NOTE that in the Science 2010 supplement these values are given
#log10(c(kt=0.03294, kon=0.10496e-3, koff=0.01721, ke=0.07483, kex=0.00994, kdi=0.003179, kde=0.01640, Epo=2030.19,
#				konS=2.294e-6, koffS=0.006799, kexS=0.0110, S=999.293))
#    kt       kon      koff        ke       kex       kdi       kde       Epo      konS     koffS      kexS         S 
#-1.482276 -3.978976 -1.764219 -1.125924 -2.002614 -2.497709 -1.785156  3.307537 -5.639407 -2.167555 -1.958607  2.999693 
(p0=10^thetahat) # in thetahat optimization, logs constrain p0 estimates to be positive. () echos result
tpts=c(0, 5, 20, 60, 120, 180, 240, 300) # the time points at which we have data, i.e. tpts=unique(epo$time)
system.time(for (i in 1:50) yout <-ode(y = c(Epo=2178,EpoR=645,EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), 
					times = tpts, 
					func = raue10, parms = p0))  # do it 50 times to make it take long enough to measure 
with(epo,matplot(time,epo[,-1]))  # time is the first column, so all but this are plotted on the y
matlines(yout[,1],yout[,-(1:7)])  # now plot the fits to epo.e, m and i, i.e. the auxillary variables y1, y2, y3

(f=file.path(system.file(paste("libs",Sys.getenv("R_ARCH"),sep=""), package = "myelo"),
					paste("myelo",.Platform$dynlib.ext,sep="")))
dyn.load(f)
system.time(for (i in 1:50) yout <-ode(y = c(Epo=2178,EpoR=645,EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), 
					times = tpts, 
					func = "derivsEpo", parms = p0, 
					dllname = "myelo",initfunc = "parmsEpo",
					nout = 3, outnames = c("y1", "y2","y3")) )
with(epo,matplot(time,epo[,-1]))
matlines(yout[,1],yout[,-(1:7)])
# so the C code is ~40 times faster, and it yields the same result. 
# so from here on out we only use the C coded ODEs

########## The rest of this demo concerns fitting the model to the data #########################
# for clarification, in the nonlinear least squares objective functions defined below 
np=as.numeric(p0) # strips names before sending values in as ODE ICs 
c(Epo=p0[2],EpoR=p0[1]) #to avoid this sort of thing
c(Epo=np[2],EpoR=np[1]) # and instead have this sort of thing
#Also, there are three model outputs to fit to data, so residuals need to be 
# weighted by the reciprocols of the measurement error, obtained here using
# saturated one-way anova models of the time courses. 
sig=c(summary(lm(epo.e~as.factor(time),data=epo))$coef[2,2],
		summary(lm(epo.m~as.factor(time),data=epo))$coef[2,2],
		summary(lm(epo.i~as.factor(time),data=epo))$coef[2,2])
(Sig=rep(1,24)%*%t(sig))  # make the sd's match the dims of the data


# first let's see how this is done using optim directly
fopt<-function(p,d,ts,IW)   # p = parameters, d = data; Note that tpts and Sig are passed globally
{	p0=10^p
	np=as.numeric(p0) # strip names
	yout <- ode(y = c(Epo=np[2],EpoR=np[1],EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), times = ts, 
			func = "derivsEpo", parms = p0, dllname = "myelo",initfunc = "parmsEpo",
			nout = 3, outnames = c("y1", "y2","y3")) 
	E=cbind(rep(yout[,8],each=3),rep(yout[,9],each=3),rep(yout[,10],each=3))
	res=(E-d[,2:4])/IW
	sum(res^2)+(np[2]-2100)^2/44100  # = Epo0 prior of 2100 with weight of 1 datapoint
}

strt=c(Bmax=3,Epo0=3,Kd=3,kde=-3,kdi=-3,ke=-3,kex=-3,kon=-3,kt=-3,scale=3)
ssolm=optim(strt,fopt,control=list(maxit=1e4),d=epo,ts=tpts,IW=Sig,hessian=TRUE)
(sig=sqrt(diag(solve(ssolm$hessian))))
ssolm$par
cbind(strt,thetahat,point=ssolm$par,lower=ssolm$par-1.96*sig,upper=ssolm$par+1.96*sig)
# so the parameter vector moved from its initial value toward Table 1 values,
# but it did not hit them dead on. What about the fit

myplot=function(p) {
	p0=10^p
	np=as.numeric(p0) # strip names
	yout <-ode(y = c(Epo=np[2],EpoR=np[1],EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), 
			times = tpts, func = "derivsEpo", parms = p0, 
			dllname = "myelo",initfunc = "parmsEpo",
			nout = 3, outnames = c("y1", "y2","y3")) 
	with(epo,matplot(time,epo[,-1]))
	matlines(yout[,1],yout[,-(1:7)])
}

myplot(ssolm$par) # the fit looks as good by eye as the Table 1 fit 
c(fopt(thetahat,epo,tpts,Sig),fopt(ssolm$par,epo,tpts,Sig)) # but the Table 1 fit is actually better.

# Conclusion: between roughly equally good fits yet much distance between the 
# optimum parameter values in Table 1 and those found here, and between there also being many NaN 
# in the CI, one suspects that the model may be over-parameterized, and must therefore question the
# the tightness of the CI given in Table 1. Oddly, such tight CI are also found below if we start 
# close to the Table 1 values. 


# mle2 is an optim wrapper in bbmle (Ben Bolker's MLE package). It has several conveniences. The main one
# needed here is his likelihood profile method. 
require(bbmle) 
nLL<-function(Bmax,Epo0,Kd,kde,kdi,ke,kex,kon,kt,scale,d,ts,IW)  	# tpts and Sig still passed globablly
{	p0=10^c(Bmax,Epo0,Kd,kde,kdi,ke,kex,kon,kt,scale)
	np=as.numeric(p0) # strip names for optimized IC parameters. Bmax is both a param and EpoR0
	yout <- ode(y = c(Epo=np[2],EpoR=np[1],EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), times = ts, 
			func = "derivsEpo", parms = p0, dllname = "myelo",initfunc = "parmsEpo",
			nout = 3, outnames = c("y1", "y2","y3")) 
	E=cbind(rep(yout[,8],each=3),rep(yout[,9],each=3),rep(yout[,10],each=3))
	res=(E-d[,2:4])/IW
#				res[1:3,2]=0 # take out of sse since ic = zero => no headway
	sum(res^2)+2*(np[2]-2100)^2/44100  # 44100=210^2  soft constraint = prior on Epo IC 
}


fit0 <- mle2(nLL,start=as.list(strt),control=list(maxit=1e4),
#		method="Nelder-Mead",  # same as optim default, so should be the same,
		method="L-BFGS-B", # brings SSE down by ~100 
		data=list(d=epo,ts=tpts,IW=Sig))
(s=summary(fit0))   # same trouble with CI, 
myplot(coef(s)[,1])  # fit looks fine
cbind(mle2=coef(s)[,1],optim=ssolm$par,thetahat) # close to optim values. Not sure why they differ.
c(fopt(thetahat,epo,tpts,Sig),fopt(ssolm$par,epo,tpts,Sig),fopt(coef(s)[,1],epo,tpts,Sig)) 
# Table 1 fit is still best

# now start at Table 1 point estimate
fit0 <- mle2(nLL,start=as.list(thetahat),control=list(maxit=1e4),
		method="L-BFGS-B", 
		data=list(d=epo,ts=tpts,IW=Sig))
(s=summary(fit0))   # trouble with CIs is now gone, 
myplot(coef(s)[,1])  # fit looks fine
cbind(mle2=coef(s)[,1],optim=ssolm$par,thetahat) # close to optim values. Not sure why they differ.
c(fopt(thetahat,epo,tpts,Sig),fopt(ssolm$par,epo,tpts,Sig),fopt(coef(s)[,1],epo,tpts,Sig)) 
# now better than Table 1 fit 
myplot(coef(s)[,1])  # fit looks great, but the std devs seem too tight given 
# all of the previous CI trouble. 

# uncomment the next block to run the profile likelihoods.

#system.time(p2<-profile(fit0,alpha=0.05,maxsteps=20)) 
## this takes ~900 seconds = 15 minutes = too long => use multicores
## it also gives a max step warning about Kd, but the Kd profile below was computed
#plot(p2)
#confint(p2) # CI still seem too tight, and indeed, they are relative to Table 1


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

# As prep for parallelization, if I do this loop on the parameters
#system.time({
#			myprof=function(k) profile(fit0,which=k,alpha=0.05,maxsteps=20)
#			L=lapply(1:10,myprof)
#		})
# it crashes with "hit maximum number of steps (Kd)"
# after 880 secs with max steps=20 and after 1400 secs with max steps =100

# For kicks I tried it with parameter limits set to those seen in Figure 3 of Chaos 2010 
#lowLim=c(Bmax=2.7,Epo0=3.22,Kd=1.4,kde=-1.95,kdi=-3.1,ke=-1.205,kex=-2.8,kon=-4.25,kt=-1.84,scale=2.13)
#upLim=c(Bmax=2.95,Epo0=3.42,Kd=3,kde=-1.82,kdi=-2.5,ke=-1.145,kex=-2.2,kon=-3.96,kt=-1.66,scale=2.33)
#fit0 <- mle2(nLL,start=as.list(thetahat), 
#		method="L-BFGS-B",
#		lower=lowLim,
#		upper=upLim,
#		data=list(d=epo,ts=tpts,IW=Sig))
#(s=summary(fit0))   
#myplot(coef(s)[,1]) 
#
#system.time(p2<-profile(fit0,alpha=0.05,maxsteps=20)) # fails because it finds a lower optimum
#system.time(p2<-profile(fit0,alpha=0.05,maxsteps=20,prof.lower=lowLim,prof.upper=upLim)) # same here
#system.time({myprof=function(k) profile(fit0,which=k,alpha=0.05,prof.lower=lowLim,prof.upper=upLim,maxsteps=20)
#			 L=lapply(1:10,myprof)}) #this also fails with maxstep error on Kd
#
# this gives an error, and may not work with RHS in C
#fit0 <- mle2(nLL,start=as.list(thetahat[-3]), 
#		fix=list(Kd=2.583), 
#		method="Nelder-Mead", 
#		data=list(x=epo))
# also tried to bypass this with objective function that ends as 
# res^2+(np[2]-2178)^2+(np[1]-644.7)^2 +(np[3]-383.4)^2  # i.e. constraint = strong priors
# but this didn't fix things. 
#
# So I hard coded (HC) the RHS with these "data" values fixed within it.
c(Epo0=2100,Kd=164,Bmax=516)  #From the Excel data file of Raue et al 2010. 
# These differ somewhat from the optimum values in thetahat in Table 1
#Epo0        Kd      Bmax 
#2178.4370  383.4184  644.7391 


nLLHC<-function(kde,kdi,ke,kex,kon,kt,scale,d,ts,IW)  	
{	p0=10^c(kde,kdi,ke,kex,kon,kt,scale)
	yout <- ode(y = c(Epo=2100,EpoR=516,EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), times = ts, 
			func = "derivsEpoHC", parms = p0, dllname = "myelo",initfunc = "parmsEpoHC",
			nout = 3, outnames = c("y1", "y2","y3")) 
	E=cbind(rep(yout[,8],each=3),rep(yout[,9],each=3),rep(yout[,10],each=3))
	res=(E-d[,2:4])/IW
	sum(res^2) 
}

thetahatHC=c(kde  = -1.884,kdi=-2.730, ke= -1.177, kex=-2.447,kon=-4.091,kt=-1.758,scale= 2.210)  
fit0 <- mle2(nLLHC,start=as.list(thetahatHC), 
		method="L-BFGS-B", # other methods are "Nelder-Mead" and "SANNN" 
		data=list(d=epo,ts=tpts,IW=Sig))
(s=summary(fit0))   
#				Estimate Std. Error  z value     Pr(z)    
#		kde   -1.8446003  0.0119058 -154.933 < 2.2e-16 ***
#		kdi   -2.6776102  0.0426771  -62.741 < 2.2e-16 ***
#		ke    -1.1959008  0.0085628 -139.662 < 2.2e-16 ***
#		kex   -2.6041774  0.0826259  -31.518 < 2.2e-16 ***
#		kon   -3.9358728  0.0093823 -419.500 < 2.2e-16 ***
#		kt    -1.6261053  0.0116032 -140.143 < 2.2e-16 ***
#		scale  2.2073982  0.0011478 1923.220 < 2.2e-16 ***


myplotHC=function(p) {
	p0=10^p
	yout <-ode(y = c(Epo=2100,EpoR=516,EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), 
			times = tpts, func = "derivsEpoHC", parms = p0, 
			dllname = "myelo",initfunc = "parmsEpoHC",
			nout = 3, outnames = c("y1", "y2","y3")) 
	with(epo,matplot(time,epo[,-1]))
	matlines(yout[,1],yout[,-(1:7)])
}
myplotHC(coef(s)[,1])  # fit still looks good
#system.time(p2<-profile(fit0,alpha=0.05,maxsteps=20)) # and we're down to just 140 seconds now
#plot(p2)
#confint(p2) # 

#let's try with IC at neighboring whole integers
thetahatHC=c(kde  = -2,kdi=-3, ke= -1, kex=-2,kon=-4,kt=-2,scale= 2)  
fit0 <- mle2(nLLHC,start=as.list(thetahatHC), 
		method="L-BFGS-B", 
		data=list(d=epo,ts=tpts,IW=Sig))
(s=summary(fit0))
# 		 		Estimate Std. Error  z value     Pr(z)    
#       kde   -1.8446277  0.0119055 -154.939 < 2.2e-16 ***
#		kdi   -2.6777016  0.0426833  -62.734 < 2.2e-16 ***
#		ke    -1.1959091  0.0085616 -139.684 < 2.2e-16 ***
#		kex   -2.6041994  0.0826086  -31.525 < 2.2e-16 ***
#		kon   -3.9358831  0.0093821 -419.510 < 2.2e-16 ***
#		kt    -1.6261252  0.0116025 -140.153 < 2.2e-16 ***
#		scale  2.2073988  0.0011478 1923.226 < 2.2e-16 ***

#yes, optimum now looks robust to changing initial parameter values!

# one last step, let's try to cut the time in half with parallel processing. 
# For this we need to get into the lapply mindset, i.e. doing things like this
system.time({myprof=function(k) profile(fit0,which=k,alpha=0.05,maxsteps=20)
			L=lapply(1:7,myprof)})  # also takes 140 secs, as expected, so should be right
par(mfrow=c(2,4))
lapply(L,plot)  # looks good
L


#install.packages("snow")
require(snow)
ptype="SOCK"  # sockets
cpusPerHost=c("localhost" = 2)
cpus=sum(cpusPerHost)
hosts=names(cpusPerHost)
strn=rep(hosts,times=cpusPerHost)
cl <- makeCluster(strn, type=ptype, verbose=TRUE) # this one can work either way
clusterCall(cl, function() getwd())
clusterCall(cl, function() Sys.info()[c("nodename","machine")])
clusterEvalQ(cl, require(myelo))
clusterEvalQ(cl, require(bbmle))
clusterExport(cl, ls())
clusterCall(cl, function() ls(.GlobalEnv))
#clusterExport(cl, c("fit0","myprof","thetahatHC","nLLHC","tpts","Sig","f"))
#clusterCall(cl, function() ls(.GlobalEnv))
clusterCall(cl, function() myprof)
clusterCall(cl, function() fit0)
clusterCall(cl, function() epo)
clusterCall(cl, function() dyn.load(f))
clusterCall(cl, function() profile(fit0,which=1,alpha=0.05,maxsteps=20))
system.time(PL<-clusterApply(cl, 1:7, myprof))  # takes ~95 secs, so not quite cut in half from 140
par(mfrow=c(2,4))                               # but 4 jobs to one cpu and 3 to the other, so expect this
lapply(PL,plot)          # on a quad proc the time was 46 secs, i.e. 3 cpus with 2 to do, 1 with 1.

# how much does load balancing the 7 jobs help. 
system.time(LB<-clusterApplyLB(cl, 1:7, myprof))  # takes ~88 secs, so jobs each take about the same time
par(mfrow=c(2,4))
lapply(LB,plot) 

system.time(PL<-parLapply(cl, 1:7,myprof)) # also in snow, also takes 90 secs 
par(mfrow=c(2,4))                             
lapply(PL,plot) 

parLapply(cl, 1:20, get("+"), 3) # use get for symbol operators
parSapply(cl, 1:20, get("+"), 3) 
stopCluster(cl)


