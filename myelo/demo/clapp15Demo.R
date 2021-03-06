library(myelo)  
library(deSolve)
library(tidyverse)
#make Table S2
(M=as.data.frame(matrix(c(1,9.944,131.016,0.187,4.021e-8,6.001e4,1.754e5,
2,33.268,148.517,0.131,1.515e-8,1.443e4,4.521e4,
3,4.612,92.3215,0.031,9.964e-7,4.994e4,5.598e5,
4,1.456,545.150,0.099,1.504e-8,3.765e4,2.759e5,
5,1.872,1700.274,0.128,3.082e-7,2.482e4,6.513e4,
6,9.652,43.752,0.238,1.350e-7,1.050e5,8.637e5,
7,5.771,155.963,0.019,4.057e-8,4.846e4,3.695e5,
8,591.591,14.568,0.040,2.371e-7,3.132e3,2.228e4,
9,486.315,226.000,0.075,2.879e-8,3.536e2,1.684e3,
10,50.988,79.645,0.005,1.271e-6,1.182e3,5.482e4,
11,30.208,359.979,0.371,2.263e-7,4.959e3,1.353e4,
12,1.5201,265.6435,0.015,2.748e-8,2.031e4,9.352e5),ncol=7,byrow=T)))
names(M)=c("id","inh1","inh2","dz","mu","ymin","ymax")
(M=M%>%mutate(sz=120*dz,eps=1/(ymin*ymax),alf=(ymin+ymax)*eps*dz))
(d=M%>%group_by(id,ymin,ymax)%>%nest())
#load sup Table 1 info  # units are cells/mL and days
(fixPars=c(a0=0.0027,b1=0.0247,r=0.08,K=4.2872,a1=24.0005,a2=899.9820,
              d1=0.00225,d2=0.006,d3=0.0375,beta=3,x=1.5e8))
(ic=c(Y0=37.5,Y1=4.1667,Y2=1.6667e4,Y3=1.5e8,Z=120))
mymap=function(x) c(fixPars,unlist(x))
(d=d%>%mutate(pars=map(data,mymap)))
d$pars[[1]] # has all params needed below
# tools::showNonASCIIfile("R/clapp15.R")
# clapp15<-function(Time, State, Pars) {
#   with(as.list(c(Time, State, Pars)),{
#     if (Time>0) {a1=a1/inh1; a2=a2/inh2}
#     dY0 = b1*Y1 - a0*Y0                         -  mu*Y0*Z/(1+eps*Y3^2)
#     dY1 = a0*Y0 - b1*Y1 - d1*Y1 + r*Y1*(1-Y1/K) -  mu*Y1*Z/(1+eps*Y3^2)
#     dY2 = a1*Y1         - d2*Y2                 -  mu*Y2*Z/(1+eps*Y3^2)
#     dY3 = a2*Y2         - d3*Y3                 -  mu*Y3*Z/(1+eps*Y3^2)
#     dZ  = sz            - dz*Z                  + alf*Y3*Z/(1+eps*Y3^2)
#     list(c(dY0,dY1,dY2,dY3,dZ),c(ratio=log10(100*beta*Y3/(Y3+2*x))))
#   }) 
# }

fsim=function(x) ode(y = ic, times = seq(0,10*365,10), func = clapp15, parms = x) 
d=d%>%mutate(out=map(pars,fsim))
head(d$out[[1]])
fday2mo=function(x) {x[,1]=x[,1]/30.5; x}
d=d%>%mutate(out=map(out,fday2mo))
head(d$out[[1]])
for (i in 1:12) plot(d$out[[i]],which=c("ratio"),xlab="Months",ylab="log10(BCR-ABL%)",main=paste("Patient",i))  
for (i in 1:12) plot(d$out[[i]],which=c("Z"),xlab="Months",ylab="T-Cells/mL",main=paste("Patient",i))  
for (i in 1:12) plot(d$out[[i]],which=c("Y3"),xlab="Months",ylab="CML Cells/mL",main=paste("Patient",i),log="y")  
