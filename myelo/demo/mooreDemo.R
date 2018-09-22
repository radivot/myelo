library(myelo)  
library(deSolve)
library(rodeoExt)
(parameters=c(sn=0.073,dn=0.04,de=.06,dc=0.2, kn=0.001, eta=100, 
              alfn=.41, alfe=0.2,Cmax=3e5,rc=0.03,ge=0.005,gc=0.005)) 
# units are cells/uL and days
(ic=c(Tn=1510,Te=20,C=1e4)) 

graphics.off()
out=ode(y = ic, times = seq(0,365,20), func = moore04, parms = parameters) 
out 
graphics.off()
quartz(width=8,height=7)
plot(out,which=c("T"),xlab="Days",ylab="Cell/uL")  
plot(out,which=c("C"),xlab="Days",ylab="Cell/uL")  


