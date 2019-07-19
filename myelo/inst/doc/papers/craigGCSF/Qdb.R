rm(list=ls())
library(tidyverse)  
library(deSolve)
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

(parsQdb=c(tauS=2.8, fQ = 2*betaSS, thresh = 0.1, betaSS = betaSS, kapDel = kapDel))

times <- seq(-20,500,by=0.01)
(x05=c( Q = 1.10216835127605,   S1 = 0.0330752837276347, 
        S2 = 0.0330752837276349, S3 = 0.0330752837276351, S4 = 0.0330752837276354, 
        Aq = 1.5116))

system.time(yout5 <- dede(x05,times = times, func = "derivsQdb",	parms = parsQdb,
                          dllname = "myelo",initfunc = "parmsQdb",
                          events=list(data=eventsdat),method="lsoda",  
                          nout = 1, outnames = c("beta"))    )

D5=data.frame(yout5)
head(D5,2)
tail(D5,2)
gx=xlab("Days")
sbb=theme(strip.background=element_blank())
cc=coord_cartesian(xlim=c(-2,125))
tc=function(sz) theme_classic(base_size=sz)
d5=D5%>%select(time,Q,Aq)%>%gather(key="Lab",value="Value",-time)
d5%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb#+cc
ggsave("~/Results/myelo/Qdb2x6.pdf",width=5, height=5)


(parsQdb1=c(fQ = 2*betaSS, thresh = 0.1, betaSS = betaSS, kapDel = kapDel))
times <- seq(-20,500,by=0.01)
(x0=c( Q = 1.10216835127605, Aq = 1.5116))

system.time(yout <- dede(x0,times = times, func = "derivsQdb1",	parms = parsQdb1,
                         dllname = "myelo",initfunc = "parmsQdb1",
                         events=list(data=eventsdat),method="lsoda",  
                         nout = 1, outnames = c("beta"))    )

D=data.frame(yout)
head(D,2)
tail(D,2)
d=D%>%select(time,Q,Aq)%>%gather(key="Lab",value="Value",-time)
d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb#+cc
ggsave("~/Results/myelo/Qdb1_2x6.pdf",width=5, height=5) #flat steps => slightly more killing (better)

### Qdb1 is twice as fast
# 
# library(rbenchmark)
# benchmark("db" = {
#   dede(x05,times = times, func = "derivsQdb",	parms = parsQdb,
#        dllname = "myelo",initfunc = "parmsQdb",
#        events=list(data=eventsdat),method="lsoda",  
#        nout = 1, outnames = c("beta")) 
# },
# "db1"={
#   dede(x0,times = times, func = "derivsQdb1",	parms = parsQdb1,
#        dllname = "myelo",initfunc = "parmsQdb1",
#        events=list(data=eventsdat),method="lsoda",  
#        nout = 1, outnames = c("beta"))    
# },
# "db1ode"={
#   ode(x0,times = times, func = "derivsQdb1",	parms = parsQdb1,
#        dllname = "myelo",initfunc = "parmsQdb1",
#        events=list(data=eventsdat),method="lsoda",  
#        nout = 1, outnames = c("beta"))    
# },
# 
# replications = 25,
# columns = c("test", "replications", "elapsed",
#             "relative", "user.self", "sys.self")
# )
# 
# 
