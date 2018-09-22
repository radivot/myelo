library(myelo)  
library(deSolve)
library(rodeoExt)
(parameters=c(sn=0.073,dn=0.04,de=.06,dc=0.2, kn=0.001, eta=100, 
              alfn=.41, alfe=0.2,Cmax=3e5,rc=0.03,ge=0.005,gc=0.005)) 
# units are cells/uL and days
(ic=c(Tn=1510,Te=20,C=1e4)) 

graphics.off()
out=ode(y = ic, times = seq(-30,0,1), func = moore04, parms = parameters) 
out #validate that IC is SS
# now add in one bad guy
(ic=c(C1=C1ss,C2=C2ss,C3=C3ss,L1=1,L2=0,L3=0)) # 1 LSC start is in 3B
out2   <- ode(y=ic, times=seq(1, 500, by = 2), func = stiehl15, parms = parameters)
out=rbind(out,out2)
graphics.off()
quartz(width=8,height=7)
# The following plot should be the same as the grey bundle in Fig. 3A, but it rises a little early (50% blasts at ~175 days)
plot(out,which=c("B"),xlab="Days",ylab="Marrow Blasts (%)")  # instead of 50% at ~200 days


