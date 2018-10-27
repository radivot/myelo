library(myelo)  
library(deSolve)
library(tidyverse)
(M=as.data.frame(matrix(   # Table 3
c(1,2.521e2,0.051,3.647e-6,6.610e4,3.624e5,
  2,1.133e3,0.026,2.405e-8,3.831e4,3.055e5,
  3,4.205e2,0.054,4.224e-7,1.617e4,3.133e5,
  4,5.691e3,0.181,8.499e-6,1.206e3,1.090e4,
  5,4.594e3,0.038,5.723e-9,1.841e3,3.401e4,
  6,2.853e3,0.058,1.358e-9,7.143e3,7.576e4),ncol=6,byrow=T)))
names(M)=c("id","inh1","d","mu","ymin","ymax")
(M=M%>%mutate(s=120*d,eps=1/(ymin*ymax),alf=(ymin+ymax)*eps*d))
M
(d=M%>%group_by(id,ymin,ymax)%>%nest())
#Table 1 universal params: units are cells/mL and days
(fixPars=c(a0=0.0027,b1=0.0247,r=0.07775,K=41.677,a1=1.35e5,
              d2=0.0375,beta=3,x=1.5e8))
(ic=c(Y1=41.677,Y2=1.5e8,Z=120))
mymap=function(x) c(fixPars,unlist(x))
(d=d%>%mutate(pars=map(data,mymap)))
d$pars[[1]] # has all params needed below
d1=0.00225
besse18<-function(Time, State, Pars) {
  with(as.list(c(Time, State, Pars)),{
    if (Time>0) a1=a1/inh1
    # dY1 = r*Y1*(1-Y1/K)                         -  mu*Y1*Z
    dY1 = r*Y1*(1-Y1/K) - d1*Y1                   -  mu*Y1*Z
    dY2 = a1*Y1         - d2*Y2                 -  mu*Y2*Z
    dZ  = s             -  d*Z                  + alf*Y2*Z/(1+eps*Y2^2)
    list(c(dY1,dY2,dZ),c(ratio=log10(100*beta*Y2/(Y2+2*x))))
  })
}

fsim=function(x) ode(y = ic, times = seq(0,10*365,1), func = besse18, parms = x) 
d=d%>%mutate(out=map(pars,fsim))
head(d$out[[1]])
fday2mo=function(x) {x[,1]=x[,1]/30.5; x}
d=d%>%mutate(out=map(out,fday2mo))
head(d$out[[1]])
for (i in 1:6) plot(d$out[[i]],which=c("ratio"),xlab="Months",ylab="log10(BCR-ABL%)",main=paste("Patient",i))  
