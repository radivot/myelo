library(tidyverse)  # figure 9 using C code
library(deSolve)
library(myelo)
(x0=craigIC)
(gtimes=as.numeric(t(outer(seq(0,80,14),4:13,"+"))))
n=length(gtimes)
(eventG=tibble(var=rep("Gs",n),
                   time=gtimes,
                   value=rep(craigPars16[["F300"]]*300e3/craigPars16[["Vd300"]],n),
                   method=rep("add",n)))
(ctimes=seq(0,80,14))
nc=length(ctimes)
dose=4*craigPars16[["BSA"]]*1e3# ug of chemo per injection (D=4 mg/m2)
dose=1*craigPars16[["BSA"]]*1e3# Q does better, but N spikes more
dose=.1*craigPars16[["BSA"]]*1e3# Q change too small
dose=.4*craigPars16[["BSA"]]*1e3# Q change too small
dose=.56*craigPars16[["BSA"]]*1e3# Q on the money, ANC spike up to 50 => need some other way to fix spike
(eventCspike=tibble(var=rep("Cp",nc),
                    time=ctimes,
                    value=rep(dose,nc),
                    method=rep("add",nc)))
(eventdat=as.data.frame(bind_rows(eventG,eventCspike)%>%arrange(time)))

(f=file.path(system.file(paste("libs",Sys.getenv("R_ARCH"),sep=""), package = "myelo"),
             paste("myelo",.Platform$dynlib.ext,sep="")))
dyn.load(f)

times <- seq(-20,85,by=0.01) 
system.time(yout <- dede(x0,times = times, func = "derivsCraig16",	parms = craigPars,
                         dllname = "myelo",initfunc = "parmsCraig16",
                         # events=list(data=eventdat),method="lsodar",
                         events=list(data=eventdat),method="lsoda",  #method made no diff
                         # events=list(data=eventdat),method="lsode",
                         # events=list(data=eventdat),method="lsodes",
                         # events=list(data=eventdat),method="radau",
                         nout = 14, outnames = c("ANC","Qts","Cpts","Qtn","G1tn","Cptn","G1tnm","Cptnm",
                                                 "dAn","dTn","eta","etanm","etan","tGam"))    )

D=data.frame(yout)
tail(D,2)
d=D%>%select(time,ANC,Q)%>%gather(key="Lab",value="Value",-time)
tc=function(sz) theme_classic(base_size=sz)
gx=xlab("Days")
sbb=theme(strip.background=element_blank())
cc=coord_cartesian(xlim=c(-2,85))
d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
ggsave("~/Results/myelo/ancQfig9n10.pdf",width=5, height=5)

d=D%>%select(time,ANC)%>%gather(key="Lab",value="Value",-time)
cc=coord_cartesian(xlim=c(-2,30))
d%>%ggplot(aes(x=time,y=Value))+ylab("Absolute Neutrophil Count (1000 cells/uL)")+geom_line(size=1)+gx+tc(14)+sbb+cc
ggsave("~/Results/myelo/ancfig9.pdf",width=5, height=5)


d=D%>%select(time,Nr,G1,G2)%>%gather(key="Lab",value="Value",-time)
cc=coord_cartesian(xlim=c(-2,30))
d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
ggsave("~/Results/myelo/NrG1n2fig10.pdf",width=5, height=5)

d=D%>%select(time,ANC,Nr,An,dTn,Tn,tGam,Cptnm)%>%gather(key="Lab",value="Value",-time)
cc=coord_cartesian(xlim=c(-2,30))
d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
ggsave("~/Results/myelo/spikeCntrl.pdf",width=5, height=8)


D=data.frame(yout)
head(D,2)
D=D%>%select(time:Nr,ANC,G1:)%>%gather(key="Lab",value="Value",-time)
cc=coord_cartesian(xlim=c(-2,85))
D%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
ggsave("~/Results/myelo/Counts85C.pdf",width=5, height=8)


system.time(yout <- dede(x0,times = times, func = craig16,	parms = craigPars,
                         events=list(data=eventdat),method="lsodar"))
D=data.frame(yout)
head(D,2)
D=D%>%select(time:Nr,ANC,G1:An,Cp)%>%gather(key="Lab",value="Value",-time)
cc=coord_cartesian(xlim=c(-2,85))
D%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
ggsave("~/Results/myelo/Counts85R.pdf",width=5, height=8)


# D=data.frame(yout)
# tail(D,2)
# D=D%>%select(time,Nr:Cp,dAn:tGam)
# d=D%>%gather(key="Lab",value="Value",-time)
# tc=function(sz) theme_classic(base_size=sz)
# gx=xlab("Days")
# sbb=theme(strip.background=element_blank())
# d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
# ggsave("~/Results/myelo/dAnProb85R.pdf",width=5, height=8)
# # ggsave("~/Results/myelo/dAnProb30R.pdf",width=5, height=8)
# 
# 
# # cc=coord_cartesian(xlim=c(-12,85))#clips high errorbars
# #now try redoing using C code RHS
# 
# 
# 
# 
# 
# # plot(yout)
# myPlotDelay=function(yout,cc) {
#   D=data.frame(yout)
#   D=D%>%select(time:Q,G1,Cp,Qts:Cptnm)
#   names(D)
#   d=D%>%gather(key="Lab",value="Value",-time)
#   tail(D)
#   head(D)
#   d=d%>%mutate(Lab=factor(Lab,levels=c("Q","Qts","Qtn","G1","G1tn","G1tnm","Cp","Cpts","Cptn","Cptnm")))
#   # table(d$Lab,useNA="always")
#   # head(d)
#   tc=function(sz) theme_classic(base_size=sz)
#   gx=xlab("Days")
#   sbb=theme(strip.background=element_blank())
#   g=d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
#   print(g)
# }
# cc=coord_cartesian(xlim=c(-2,185))#clips high errorbars
# myPlotDelay(yout,cc)
# 
# myPlot=function(yout,cc) {
#   D=data.frame(yout)
#   D%>%filter(time>-.1,time<1)
#   head(D)
#   tail(D)
#   d=D%>%select(time:Cp,Cs1,ANC)%>%gather(key="Lab",value="Value",-time)%>%
#     mutate(Lab=factor(Lab,levels=c("Q","Nr","N","G1","G2","Tn","An","Aq","Cp","Cs1","ANC")))
#   tc=function(sz) theme_classic(base_size=sz)
#   gx=xlab("Days")
#   sbb=theme(strip.background=element_blank())
#   g=d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
#   # g=d%>%ggplot(aes(x=time,y=Value))+facet_wrap(Lab~.,ncol=2,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
#   print(g)
# }
# myPlot(yout,cc)
# 
# 
# 
# # thetahat=craigPars16
# # p0=10^thetahat