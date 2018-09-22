# use COPASI to save parmar sbml as XPPAUT *.ODE and do rest by hand + find/replace
library(myelo)   # load definition of function parmar17
library(deSolve)
library(rodeoExt)
(parameters=parmarPars)
(ic=c(Tf=2.05684e-08, FeDuo=4.52271e-07, FeRest=5.64928e-06, FeBM=6.50966e-07, NTBI=4.33046e-11,           #states are 
     FeLiver=2.32394e-06, FeSplee=2.65354e-07, Hepcidi=2.99023e-11, Fe2Tf=1.75823e-08, FeRBC=3.00042e-05)) #numbers in Moles
TfTot=5.0309999999999992e-08 # total transferin number in Moles

sum(parameters[1:8]) #23.2 gm mouse
graphics.off()
with(as.list(parameters),plot(c(vDietL,vDiet,vDietH),c(vHepSynL,vHepSyn,vHepSynH)))
out=ode(y = ic, times = seq(-30,0,1), func = parmar17, parms = parameters) 
(N=length(ic))
n=dim(out)[1]
X0=out[n,2:(N+1)]
names(X0)<-names(ic)
parameters["vHepSyn"]=5*parameters["vHepSyn"]
out2   <- ode(y=X0, times=seq(1, 360, by = 1), func = parmar17, parms = parameters)
out=rbind(out,out2)
head(out)
?plot.deSolve
graphics.off()
quartz(width=8,height=7)
# par(mfrow=c(2,4))
vars=c("FePlas","FeDuo_c","FeBM_c","FeRBC_c","FeLiver_c","FeSplee_c","FeRest_c","Hepcidi_c")
plot(out,which=vars,xlab="Days",ylab="Concentration (M)")

dput(out)
dput(out2)
graphics.off()
plot(out[ , 1], out[ , "FePlas"], type = "l", xlab = "time", ylab = "Iron in Duodenum",
     main = "Parmar/Mendes Fe Model", log="y", lwd = 2,col=1)
plot(out[ , 1], out[ , "FeDuo"], type = "l", xlab = "time", ylab = "Iron in Duodenum",
     main = "Parmar/Mendes Fe Model", log="y", lwd = 2,col=1)

plot(out)

tail(out)  # Mismatch with Table 3: Nm/Np = 1.9 < 2.3, Nm = 2.1e9 < 2.6e9,=one problem, c(2.3/1.9,2.6/2.1)

X0=out[700,2:15] # start at steady state, then hit with chemo dropping all cells to 1%
(eventdat <- data.frame(var = names(X0), time = 0, value = 5, method = "mult"))
out   <- ode(X0, -20:300, fokas91, pars,events = list(data = eventdat))
tail(out) 


# now do normal state
pars=c(	fb=0.8, # blast fraction remaining active per division (nb = 3 is assumed)
		fp=0.6, # progenitor fraction remaining active per division (np=2 is assumed)
		fm=0.5, # myelocyte fraction remaining active per division  (nm=2 is assumed)
		T= 20, # hours per stage
		q=72e6/16  # flux into the blast pool per hour per kg of marrow
)  # 

X0=rep(0,14)
names(X0)<-c(paste("A",1:7,sep=""),paste("M",1:7,sep=""))
times <- seq(1, 700, by = 1)
out   <- ode(X0, times, fokas91, pars)
tail(out)  # Mismatch with Table 3: Nm/Np = 1.9 < 2.3, Nm = 2.1e9 < 2.6e9,=one problem, c(2.3/1.9,2.6/2.1)

X0=out[700,2:15] # start at steady state, then hit with chemo dropping all cells to 1%
(eventdat <- data.frame(var = names(X0), time = 0, value = 0.01, method = "mult"))
out   <- ode(X0, -20:300, fokas91, pars,events = list(data = eventdat))
tail(out) 

lines(out[ , 1], out[ , "Ntot"],lwd = 2,col=2)
legend(150,1e8, c("CML","Normal"), col = 1:2, bty="n",lwd=2)
# recoveries from chemo in 7 days as shown here seem reasonable
