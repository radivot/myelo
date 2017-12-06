library(myelo)  
library(deSolve)
library(rodeoExt)
# start with C  in ss without L 
# s=1/(1+kc*C3)
# 0=(2*a1c*s-1)*p1c*C1 => s=1/2a1c = 1/(1+kc*C3ss) => kc=(2*a1-1)/C3ss
C3ss=4e8 # neuts in blood per Kg body. At 6L/70Kg, (70/6)*4e8 = 4.6e9/L = ~5k/uL, check
a1c=0.85 # from bone marrow transplantation data ([15] in supplement)  patients need ~15 days to
# engraft to 5e8 neutrophils per L of blood (4e7 per kg ) after infusion
# of 5e6 immature cells per kg of body weight. 8 fold increase => 3 doublings in 15 days.
#Guessing other parameters were set as below and a1c was then tuned to match this
(kc=(2*a1c-1)/C3ss)
# 0=2*(1-a1c*s)*p1c*C1ss + (2*a2c*s-1)*p2c*C2ss
# and s=1/(2a1c) from 1 yields
# 1*p1c*C1ss = (1-a2c/a1c)*p2c*C2ss 
# => p2c = p1c*(C1ss/C2ss)/(1-a2c/a1c) ####Eq. 6
# 0=2*(1-a2c*s)*p2c*C2 - d3c*C3
# and s=1/(2a1c) from 1 yields
# 2*(1-(1/2)(a2c/a1c))*p2c*C2ss = d3c*C3ss 
# => (2-a2c/a1c)*p2c*C2ss = d3c*C3ss 
# => (2-a2c/a1c)*p2c/d3c = C3ss/C2ss 
# => (2-a2c/a1c)*p1c*(C1ss/C2ss)/(1-a2c/a1c) = d3c*C3ss/C2ss 
# => (2-a2c/a1c)/(1-a2c/a1c) = d3c*C3ss/(C1ss*p1c) 
# => (2-a2c/a1c)= (1-a2c/a1c)*d3c*C3ss/(C1ss*p1c) 
# => (2-d3c*C3ss/(C1ss*p1c))= a2c/a1c*(1-d3c*C3ss/(C1ss*p1c)) 
# => 
d3c=2.3 #per day based on T1/2 of 7 hours
p1c=0.1 # HSC and committed progenitor average of dividing ~ once per week
C1ss=1e7 #in marrow per kg (this is an average of HSC and very early committed progenitors)
(x=d3c*C3ss/(C1ss*p1c))
(a2c=a1c*(2-x)/(1-x)) #eq 10
C2ss=2e9 #in marrow per kg
(p2c = p1c*(C1ss/C2ss)/(1-a2c/a1c))


(parameters=c(a1c=a1c,    a2c=a2c,  p1c=p1c,    p2c=p2c,    d3c=d3c,      kc=kc,
          a1l=1.6*a1c,a2l=1*a2c,p1l=2.0*p1c,p2l=0.5*p2c,d3l=0.25*d3c, kl=kc)) 

(ic=c(C1=C1ss,C2=C2ss,C3=C3ss,L1=0,L2=0,L3=0)) # 100 LSC start is in 3B

graphics.off()
out=ode(y = ic, times = seq(-30,0,1), func = stiehl15, parms = parameters) 
out #validate that IC is SS
# now add in one bad guy
(ic=c(C1=C1ss,C2=C2ss,C3=C3ss,L1=1,L2=0,L3=0)) # 1 LSC start is in 3B
out2   <- ode(y=ic, times=seq(1, 500, by = 2), func = stiehl15, parms = parameters)
out=rbind(out,out2)
graphics.off()
quartz(width=8,height=7)
# The following plot should be the same as the grey bundle in Fig. 3A, but it rises a little early (50% blasts at ~175 days)
plot(out,which=c("B"),xlab="Days",ylab="Marrow Blasts (%)")  # instead of 50% at ~200 days

#now try to recreat Fig 3B
(ic=c(C1=C1ss,C2=C2ss,C3=C3ss,L1=100,L2=0,L3=0)) # 1 LSC start is in 3B
(parameters=c(a1c=a1c,    a2c=a2c,  p1c=p1c,    p2c=p2c,    d3c=d3c,      kc=kc,
              a1l=2.0*a1c,a2l=1*a2c,p1l=1.6*p1c,p2l=0.5*p2c,d3l=0.25*d3c, kl=kc)) 
# Problem with this might be that a1l is bigger than 1

out2   <- ode(y=ic, times=seq(1, 500, by = 2), func = stiehl15, parms = parameters)
out=rbind(out,out2)
graphics.off()
quartz(width=8,height=7)
plot(out,which=c("B"),xlab="Days",ylab="Marrow Blasts (%)")


