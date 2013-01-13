library(myelo) 

# first just fit the PK data

rsPK<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				dGt= -kt*Gt      # G-CSF in tissue, versus blood below
				dGb=  kt*Conv*Gt -kb*Gb # - sig*Nb*Gb^2/(kG+Gb^2)  # see if data are rich enough to estimate these params
				list(c(dGt,dGb))
			})
}

# 5 ug/kg injection => 350 ug / 70 kg adult => 350 ug/ 5000 mL of adult blood 
#                   => 350e3 ng/ 5e3 mL = 70 ng/mL expected if all goes to  blood from tissue
# So, if we keep Gt in ug/kg units, an initial estimate of the conversion factor Conv is 70/5 = 14
# But some drug may never make it to the blood, so we will let Conv be estimated from the data
# and the ratio relative to 14 will be viewed as the bioavailability. 


library(myelo)
graphics.off()
(df=auwera[!is.na(auwera[,2]),1:3])
with(df,matplot(hours,df[,2:3],ylab="GCSF in blood plasma"))
title("1 = 5 ug/kg           2 = 10 ug/kg")
# Looking at the G-CSF kinetics, the rising edge time constant is perhaps ~3 hours and the falling edge is ~6
# so initial kt and kb estimates might be 1/3 and 1/6, respectively

times <- c(0:40)
yout5 <- ode(c(Gt=5,Gb=0),	times = times, func = rsPK,parms = c(kt=1/3,kb=1/6,Conv=14)) 
yout10 <- ode(c(Gt=10,Gb=0),	times = times, func = rsPK,parms = c(kt=1/3,kb=1/6,Conv=14)) 
lines(times,yout5[,3])
lines(times,yout10[,3])

# so it looks like less than half is bioavailable, but let's see what a fit says
# since some is being removed from the blood as it comes in. 

# First, to automate plotting based on log10 parameter estimates, define myplot as
myplot=function(p) {
	p0=10^p
	np=as.numeric(p0) # strip names
	yout5 <- ode(c(Gt=5,Gb=0),times = times, func = rsPK, parms = c(kt=np[1],kb=np[2],Conv=np[3])) 
	yout10 <- ode(c(Gt=10,Gb=0),times = times, func = rsPK, parms = c(kt=np[1],kb=np[2],Conv=np[3])) 
	with(df,matplot(hours,df[,2:3],ylab="GCSF in blood plasma"))
	title("1 = 5 ug/kg           2 = 10 ug/kg")
	lines(times,yout5[,3])
	lines(times,yout10[,3])
}

strt=log10(c(kt=1/3,kb=1/6,Conv=14))  # let the log10 values be our starting points
myplot(strt)   # show that new function works fine

# now use bbmle to fit the PK model to the PK data of Auwera
require(bbmle) 
nLL<-function(kt,kb,Conv)  	
{	p0=10^c(kt,kb,Conv)
	np=as.numeric(p0) # strip names for optimized IC parameters.
	yout5 <- ode(c(Gt=5,Gb=0),times = d$hours, func = rsPK, parms = c(kt=np[1],kb=np[2],Conv=np[3])) 
	yout10 <- ode(c(Gt=10,Gb=0),times = d$hours, func = rsPK, parms = c(kt=np[1],kb=np[2],Conv=np[3])) 
	E=cbind(yout5[,c(1,3)],yout10[,3])
	res=E-d
	sum(res^2)
}


fit0 <- mle2(nLL,start=as.list(strt),control=list(maxit=1e4),
#		method="Nelder-Mead",  # same as optim default, so should be the same,
		method="L-BFGS-B", # brings SSE down by ~100 
		data=list(d=df))
(s=summary(fit0))   # same trouble with CI, 
myplot(coef(s)[,1])  # fit looks fine
(estimates=cbind(mle=10^coef(s)[,1]))
estimates["Conv",1]/14  # indeed, 58% (i.e. more than half) is bioavailable


##########################################################
##########################################################

# what about increases in the removal rate with increases in ANC? Is the data rich enough to see this from it?


rsPK2<-function(Time, State, Pars,Nb) {
	with(as.list(c(State, Pars)), {
				dGt= -kt*Gt      # G-CSF in tissue, versus blood below
				dGb=  kt*Conv*Gt -kb*Gb - sig*Nb(Time)*Gb^2/(kG+Gb^2)  # see if data are rich enough to estimate these params
				list(c(dGt,dGb))
			})
}

(dfANC=auwera[!is.na(auwera[,8]),c(1,8,9)])
Nb5=with(dfANC,approxfun(hours,anc5,rule=2))
Nb10=with(dfANC,approxfun(hours,anc10,rule=2))

myplot2=function(p) {
	p0=10^p
	np=as.numeric(p0) # strip names
	yout5 <- ode(c(Gt=5,Gb=0),times = times, func = rsPK2, parms = c(kt=np[1],kb=np[2],Conv=np[3],sig=np[4],kG=np[5]),Nb=Nb5) 
	yout10 <- ode(c(Gt=10,Gb=0),times = times, func = rsPK2, parms = c(kt=np[1],kb=np[2],Conv=np[3],sig=np[4],kG=np[5]),Nb=Nb10) 
	with(df,matplot(hours,df[,2:3],ylab="GCSF in blood plasma"))
	title("1 = 5 ug/kg           2 = 10 ug/kg")
	lines(times,yout5[,3])
	lines(times,yout10[,3])
}


# in Brooks they have kG=10 (ug/ml)^2  but our units are ng/mL so kG= (3106 ng/mL)^2, ie. nowhere close to data at ~30 ng/mL and less 
# so we'll assume units were wrong somewhere and just try 10 as an initial value. 
# for sigma I wanted a starting value that would make the flux out by ANC's start with a magnitude nearly that of renal clearance
strt2=log10(c(kt=1/5,kb=1/5,Conv=8,sig=1e-4,kG=10))  # let the log10 values be our starting points
myplot2(strt2)   # new function and parameter starting point are OK


nLL2<-function(kt,kb,Conv,sig,kG)  	
{	p0=10^c(kt,kb,Conv,sig,kG)
	np=as.numeric(p0) # strip names for optimized IC parameters.
	yout5 <- ode(c(Gt=5,Gb=0),times = times, func = rsPK2, parms = c(kt=np[1],kb=np[2],Conv=np[3],sig=np[4],kG=np[5]),Nb=Nb5) 
	yout10 <- ode(c(Gt=10,Gb=0),times = times, func = rsPK2, parms = c(kt=np[1],kb=np[2],Conv=np[3],sig=np[4],kG=np[5]),Nb=Nb10) 
	E=cbind(yout5[,c(1,3)],yout10[,3])
	res=E-d
	sum(res^2)
}


fit2 <- mle2(nLL2,start=as.list(strt2),control=list(maxit=1e4),
		method="Nelder-Mead",  
#		method="L-BFGS-B",  # stuck at IC irrespective of optim algoritm used
		data=list(d=df))
(s2=summary(fit2))   
myplot2(coef(s2)[,1])  # fit looks like initial estimate plot. 
(estimates=cbind(mle2=10^coef(s2)[,1])) # yes, estimates stay at IC, so overparameterized


##########################################################
##########################################################
# so try with one less parameter, i.e. with
rsPK3<-function(Time, State, Pars,Nb) {
	with(as.list(c(State, Pars)), {
				dGt= -kt*Gt      # G-CSF in tissue, versus blood below
				dGb=  kt*Conv*Gt -kb*Gb - sig*Nb(Time)*Gb 
				list(c(dGt,dGb))
			})
}

myplot3=function(p) {
	p0=10^p
	np=as.numeric(p0) # strip names
	yout5 <- ode(c(Gt=5,Gb=0),times = times, func = rsPK3, parms = c(kt=np[1],kb=np[2],Conv=np[3],sig=np[4]),Nb=Nb5) 
	yout10 <- ode(c(Gt=10,Gb=0),times = times, func = rsPK3, parms = c(kt=np[1],kb=np[2],Conv=np[3],sig=np[4]),Nb=Nb10) 
	with(df,matplot(hours,df[,2:3],ylab="GCSF in blood plasma"))
	title("1 = 5 ug/kg           2 = 10 ug/kg")
	lines(times,yout5[,3])
	lines(times,yout10[,3])
}

strt3=log10(c(kt=1/5,kb=1/5,Conv=8,sig=1e-4)) 
myplot3(strt3)   # new function and parameter starting point are OK


nLL3<-function(kt,kb,Conv,sig)  	
{	p0=10^c(kt,kb,Conv,sig)
	np=as.numeric(p0) # strip names for optimized IC parameters.
	yout5 <- ode(c(Gt=5,Gb=0),times = times, func = rsPK3, parms = c(kt=np[1],kb=np[2],Conv=np[3],sig=np[4]),Nb=Nb5) 
	yout10 <- ode(c(Gt=10,Gb=0),times = times, func = rsPK3, parms = c(kt=np[1],kb=np[2],Conv=np[3],sig=np[4]),Nb=Nb10) 
	E=cbind(yout5[,c(1,3)],yout10[,3])
	res=E-d
	sum(res^2)
}


fit3 <- mle2(nLL3,start=as.list(strt3),control=list(maxit=1e4),
#		method="Nelder-Mead", 
		method="L-BFGS-B", 
		data=list(d=df))
(s3=summary(fit3))  
myplot3(coef(s3)[,1]) 
(estimates=cbind(mle2=10^coef(s2)[,1])) # estimates still stay at IC => still overparameterized 



