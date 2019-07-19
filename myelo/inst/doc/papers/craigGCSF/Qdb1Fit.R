rm(list=ls())
library(tidyverse)  
library(deSolve)
library(bbmle)
library(myelo)
(x0=craigIC[c(1,8)])
(parsQ=craigPars[c("Qss","Aqss","tauS","fQ","the2","s2")])
parsQ["kapDel"]=craigPars["kapss"]+craigPars["kapDel"]
parsQ

attach(as.list(parsQ))
fbeta=function(Q) fQ/(1+(Q/the2)^s2)
betaSS=fbeta(Qss)
(kapDel=(Aqss-1)*betaSS)
detach(as.list(parsQ))

(StrTimes=seq(0,80,14))
(StpTimes=StrTimes+5)
nc=length(StpTimes)
(events=tibble(var=rep("Aq",2*nc),
                    time=sort(c(StrTimes,StpTimes)),
                    value=rep(c(0.0*parsQ["Aqss"],parsQ["Aqss"]),nc),
                    method=rep("rep",2*nc)))
events2=events
events2$time=events2$time+150
(eventsdat=as.data.frame(bind_rows(events,events2)))
(f=file.path(system.file(paste("libs",Sys.getenv("R_ARCH"),sep=""), package = "myelo"),
             paste("myelo",.Platform$dynlib.ext,sep="")))
dyn.load(f)

times <- seq(-20,500,by=0.01)
# (parsQdb1=c(fQ = 2*betaSS, thresh = 0.1, betaSS = betaSS, kapDel = kapDel))
(parsQdb1=c(fQ = 0.1, thresh = 0.1, betaSS = betaSS, kapDel = kapDel))
(x0=c( Q = 1.10216835127605, Aq = 1.5116))
yout=dede(x0,times = times, func = "derivsQdb1",	parms = parsQdb1,
                          dllname = "myelo",initfunc = "parmsQdb1",
                          events=list(data=eventsdat),method="lsoda",  
                          nout = 1, outnames = c("beta"))   
D=data.frame(yout)
head(D,2)
tail(D,2)
d=D%>%select(time,Q,Aq)%>%gather(key="Lab",value="Value",-time)
gx=xlab("Days")
tc=function(sz) theme_classic(base_size=sz)
sbb=theme(strip.background=element_blank())
d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb
# ggsave("~/Results/myelo/Qdb1_2x6.pdf",width=5, height=5) 
head(d)
d=d%>%filter(Lab=="Q",time%%1==0)%>%select(t=time,Q=Value)
d%>%ggplot(aes(x=t,y=Q))+geom_line(size=1)+gx+tc(14)
d[1:100,] #plenty of samples on each plateau
d$Q=d$Q+rnorm(dim(d)[1],sd=0.01)
d[1:100,] 
d%>%ggplot(aes(x=t,y=Q))+geom_line(size=1)+gx+tc(14) #barely see the noise. 
# ggsave("~/Results/myelo/Qdb1_2x6ns.pdf",width=5, height=5) 
# now fit Q0, betaSS, and beta0 to this new simulated data
# should get Q = 1.10216835127605,  betaSS = 0.04303868, and beta0=0.1

head(d)
head(eventsdat)
Aqss=1.5116
nLL<-function(Q0,beta0,betaSS,sigma) { # pass these globally: d,eventsdat,Aqss
  (x0=c( Q = Q0, Aq = Aqss))
  kapDel=(Aqss-1)*betaSS
  parsQdb1=c(fQ = beta0, thresh = 0.1, betaSS = betaSS, kapDel = kapDel)
  out=dede(x0,times = d$t, func = "derivsQdb1",	parms = parsQdb1,
       dllname = "myelo",initfunc = "parmsQdb1",
       events=list(data=eventsdat),method="lsoda",  
       nout = 1, outnames = c("beta"))  
  # -sum(dnorm(log(d$Q), mean=log(out[,"Q"]),sigma,log=TRUE))
  -sum(dnorm(d$Q, mean=out[,"Q"],sigma,log=TRUE)) #gets sigma estimate right at 0.002
}

IC=c(Q0=1.1,beta0=0.1,betaSS=0.04303868,sigma=0.002) # leaves at initial optimum
IC=c(Q0=.6,beta0=0.07,betaSS=0.03,sigma=0.01) # recreates opt above used to make data
fit=c(1:4)
summary(M<-mle2(nLL,method="Nelder-Mead",
                 fixed=as.list(IC)[-fit],
                 start=as.list(IC),data=list(d=d),
                 control = list(maxit=50000, parscale=IC[fit]) ) )
# 
# Coefficients:
#   Estimate Std. Error  z value     Pr(z)    
#   Q0     1.1021e+00 3.6840e-04 2991.637 < 2.2e-16 ***
#   beta0  9.9462e-02 2.9259e-04  339.931 < 2.2e-16 ***
#   betaSS 4.3025e-02 2.2473e-05 1914.491 < 2.2e-16 ***
#   sigma  1.9685e-03 6.0940e-05   32.303 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# -2 log L: -5013.127 


psim=function(Q0,beta0,betaSS,sigma) {
  (x0=c( Q = Q0, Aq = Aqss))
  kapDel=(Aqss-1)*betaSS
  parsQdb1=c(fQ = beta0, thresh = 0.1, betaSS = betaSS, kapDel = kapDel)
  data.frame(dede(x0,times = d$t, func = "derivsQdb1",	parms = parsQdb1,
           dllname = "myelo",initfunc = "parmsQdb1",
           events=list(data=eventsdat),method="lsoda",  
           nout = 1, outnames = c("beta"))) 
}


signif(coef(summary(M))[,1:2],3)
round(coef(M),5)
p=as.list(coef(M))
out=do.call(psim,p)
d$ei=out[,"Q"]
gy=ylab(quote(paste(10^6," HSC/(kg body mass)")))
d%>%ggplot(aes(x=t,y=Q))+geom_point(size=0.3)+geom_line(aes(x=t,y=ei),size=.3)+gx+gy+tc(14)
ggsave("~/Results/myelo/Qdb1_2x6fit.pdf",width=4, height=4) 


