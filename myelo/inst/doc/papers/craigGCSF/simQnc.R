library(tidyverse)  
library(deSolve)
library(myelo)
(x0=craigIC[c(1,8)])
(parsQ=craigPars[c("Qss","Aqss","tauS","fQ","the2","s2")])
parsQ["kapDel"]=craigPars["kapss"]+craigPars["kapDel"]
parsQ
(StrTimes=seq(0,80,14))
(StpTimes=StrTimes+5)
nc=length(StpTimes)
(events=tibble(var=rep("Aq",2*nc),
                    time=sort(c(StrTimes,StpTimes)),
                    value=rep(c(0.8*parsQ["Aqss"],parsQ["Aqss"]),nc),
                    method=rep("rep",2*nc)))
(eventsdat=as.data.frame(events))

(f=file.path(system.file(paste("libs",Sys.getenv("R_ARCH"),sep=""), package = "myelo"),
             paste("myelo",.Platform$dynlib.ext,sep="")))
dyn.load(f)

# times <- seq(-20,1750,by=0.01) 
# system.time(yout <- dede(x0,times = times, func = "derivsCraig16Qnc",	parms = parsQ,
#                          dllname = "myelo",initfunc = "parmsCraig16Qnc",
#                          events=list(data=eventsdat),method="lsoda",  #method made no diff
#                          nout = 2, outnames = c("Qts","Aqts"))    )
# D=data.frame(yout)
# dput(tail(D,1))


x0=c(Q = 1.10216835127603, Aq = 1.5116)

times <- seq(-20,175,by=0.01) 
system.time(yout <- dede(x0,times = times, func = "derivsCraig16Qnc",	parms = parsQ,
                         dllname = "myelo",initfunc = "parmsCraig16Qnc",
                         events=list(data=eventsdat),method="lsoda",  #method made no diff
                         nout = 2, outnames = c("Qts","Aqts"))    )

D=data.frame(yout)
tail(D,2)
d=D%>%select(time,Q,Aq)%>%gather(key="Lab",value="Value",-time)
tc=function(sz) theme_classic(base_size=sz)
gx=xlab("Days")
sbb=theme(strip.background=element_blank())
cc=coord_cartesian(xlim=c(-2,125))
d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb#+cc
ggsave("~/Results/myelo/Qnc.pdf",width=5, height=5)



(events1=tibble(var=rep("Aq",2*nc),
               time=sort(c(StrTimes,StpTimes)),
               value=rep(c(0.0*parsQ["Aqss"],parsQ["Aqss"]),nc),
               method=rep("rep",2*nc)))
events2=events1
events2$time=events2$time+150
(eventsdat2=as.data.frame(bind_rows(events1,events2)))

x0=c(Q = 1.10216835127603, Aq = 1.5116)
system.time(yout <- dede(x0,times = seq(-20,500,by=0.01), func = "derivsCraig16Qnc",	parms = parsQ,
                         dllname = "myelo",initfunc = "parmsCraig16Qnc",
                         events=list(data=eventsdat2),method="lsoda",  #method made no diff
                         nout = 2, outnames = c("Qts","Aqts"))    )

D=data.frame(yout)
tail(D,2)
d=D%>%select(time,Q,Aq)%>%gather(key="Lab",value="Value",-time)
d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb#+cc
ggsave("~/Results/myelo/Qnc2x6.pdf",width=5, height=5)






#################   find steady state of S1 to S4
# (x05=c(x0[1],S1=0,S2=0,S3=0,S4=0,x0[2]))
# system.time(yout5 <- dede(x05,times = seq(-20,1000,by=0.01), func = "derivsCraig16Q5nc",	parms = parsQ,
#                           dllname = "myelo",initfunc = "parmsCraig16Q5nc",
#                           events=list(data=eventsdat),
#                           method="lsoda",
#                           nout = 1, outnames = c("beta"))    )
# 
# D5=data.frame(yout5)
# head(D5,2)
# dput(tail(D5,1))
# (x05=c(x0[1],S1=0.03307528,S2=0.03307528,S3=0.03307528,S4=0.03307528,x0[2]))
(x05=c(Q = 1.10216835127605, S1 = 0.0330752837276347, 
       S2 = 0.0330752837276349, S3 = 0.0330752837276351, S4 = 0.0330752837276354, 
       Aq = 1.5116))

system.time(yout5 <- dede(x05,times = times, func = "derivsCraig16Q5nc",	parms = parsQ,
                         dllname = "myelo",initfunc = "parmsCraig16Q5nc",
                         events=list(data=eventsdat),method="lsoda",  
                         nout = 1, outnames = c("beta"))    )

D5=data.frame(yout5)
head(D5,2)
tail(D5,2)
d5=D5%>%select(time,Q,Aq)%>%gather(key="Lab",value="Value",-time)
d5%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb
ggsave("~/Results/myelo/Q5nc.pdf",width=5, height=5)


# (x10=c(x0[1],S1=0.03307528,S2=0.03307528,S3=0.03307528,S4=0.03307528,
#        S5=0.03307528,S6=0.03307528,S7=0.03307528,S8=0.03307528,S9=0.03307528,x0[2]))
# system.time(yout10 <- dede(x10,times = seq(-20,1000,by=0.01), func = "derivsCraig16Q10nc",	parms = parsQ,
#                           dllname = "myelo",initfunc = "parmsCraig16Q10nc",
#                           events=list(data=eventsdat),
#                           method="lsoda",
#                           nout = 1, outnames = c("beta"))    )
# 
# D10=data.frame(yout10)
# head(D10,2)
# tail(D10,2)
# dput(tail(D10,1))

x10=c(Q = 1.10216835127602, S1 = 0.0147001261011713, 
S2 = 0.0147001261011713, S3 = 0.0147001261011713, S4 = 0.0147001261011713, 
S5 = 0.0147001261011714, S6 = 0.0147001261011714, S7 = 0.0147001261011714, 
S8 = 0.0147001261011714, S9 = 0.0147001261011715, Aq = 1.5116)

system.time(yout10 <- dede(x10,times = times, func = "derivsCraig16Q10nc",	parms = parsQ,
                          dllname = "myelo",initfunc = "parmsCraig16Q10nc",
                          events=list(data=eventsdat),method="lsoda",  
                          nout = 1, outnames = c("beta"))    )

D10=data.frame(yout10)
head(D10,2)
tail(D10,2)
dput(tail(D10,1))

d10=D10%>%select(time,Q,Aq)%>%gather(key="Lab",value="Value",-time)
d10%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb
ggsave("~/Results/myelo/Q10nc.pdf",width=5, height=5)

# head(d10)
# d$n=1
# d5$n=5
# d10$n=10
# D=bind_rows(d,d5,d10)%>%filter(Lab=="Q")
# D$n=factor(D$n)
# D$n=factor(D$n,levels=c(10,5,1))
# D%>%ggplot(aes(x=time,y=Value,col=n))+geom_line(size=.3)+gx+tc(14)+sbb

library(rbenchmark)
benchmark("delayQ" = {
  dede(x0,times = times, func = "derivsCraig16Qnc",	parms = parsQ,
                           dllname = "myelo",initfunc = "parmsCraig16Qnc",
                           events=list(data=eventsdat),method="lsoda",  #
                           nout = 2, outnames = c("Qts","Aqts"))    
},
"Q5"={dede(x05,times = times, func = "derivsCraig16Q5nc",	parms = parsQ,
                                     dllname = "myelo",initfunc = "parmsCraig16Q5nc",
                                     events=list(data=eventsdat),method="lsoda",  
                                     nout = 1, outnames = c("beta"))    
},
"Q10"={dede(x10,times = times, func = "derivsCraig16Q10nc",	parms = parsQ,
           dllname = "myelo",initfunc = "parmsCraig16Q10nc",
           events=list(data=eventsdat),method="lsoda",  
           nout = 1, outnames = c("beta"))    
},
replications = 25,
columns = c("test", "replications", "elapsed",
            "relative", "user.self", "sys.self")
)


