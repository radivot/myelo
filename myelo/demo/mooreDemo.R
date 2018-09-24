library(myelo)  
library(deSolve)
library(rodeoExt)
#Table 1
(parameters=c(sn=0.073, dn=0.04, de=0.06, dc=0.2,  kn=0.001, eta=100, alfn=.41, alfe=0.2, Cmax=3e5,   rc=0.03, ge=0.005, gc=0.005)) 
#Table 3
(parameters6=c(sn=0.37, dn=0.23, de=0.30, dc=0.024,kn=0.062,eta=720,alfn=0.14,alfe=0.98,Cmax=230000,rc=0.0057,ge=0.057, gc=0.0034)) 
(parameters7=c(sn=0.29, dn=0.35, de=0.40, dc=0.012,kn=0.066,eta=140,alfn=0.39,alfe=0.65,Cmax=160000,rc=0.011, ge=0.079, gc=0.058)) 
(parameters8=c(sn=0.071,dn=0.050,de=0.12, dc=0.68, kn=0.063,eta=43, alfn=0.56,alfe=0.53,Cmax=190000,rc=0.23,  ge=0.0077,gc=0.047)) 

# units are cells/uL and days
(ic=c(Tn=1510,Te=20,C=1e4)) 
graphics.off()
out=ode(y = ic, times = seq(0,750,20), func = moore04, parms = parameters6) 
graphics.off()
quartz(width=8,height=7)
plot(out,which=c("C"),xlab="Days",ylab="Cell/uL")  
out=ode(y = ic, times = seq(0,750,1), func = moore04, parms = parameters7) 
plot(out,which=c("C"),xlab="Days",ylab="Cell/uL")  
out=ode(y = ic, times = seq(0,100,1), func = moore04, parms = parameters8) 
plot(out,which=c("C"),xlab="Days",ylab="Cell/uL")  


