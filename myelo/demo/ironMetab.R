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
# ?plot.deSolve
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
