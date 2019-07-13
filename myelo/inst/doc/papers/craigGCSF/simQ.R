library(tidyverse)  
library(deSolve)
library(myelo)
(x0=craigIC[c(1,8:12,14)])
(ctimes=seq(0,80,14))
nc=length(ctimes)
dose=4*craigPars16[["BSA"]]*1e3# ug of chemo per injection (D=4 mg/m2)
dose=.56*craigPars16[["BSA"]]*1e3# Q on the money, ANC spike up to 50 => need some other way to fix spike
(eventCspike=tibble(var=rep("Cp",nc),
                    time=ctimes,
                    value=rep(dose,nc),
                    method=rep("add",nc)))
(eventdat=as.data.frame(eventCspike%>%arrange(time)))

(f=file.path(system.file(paste("libs",Sys.getenv("R_ARCH"),sep=""), package = "myelo"),
             paste("myelo",.Platform$dynlib.ext,sep="")))
dyn.load(f)
parsQ=craigPars[c("Qss","Aqss","tauS","fQ","the2","s2","kapss","kapDel",
                  "k12","k21","k24","k42","k13","k31","kelC","hQ")]
times <- seq(-20,85,by=0.01) 
system.time(yout <- dede(x0,times = times, func = "derivsCraig16Q",	parms = parsQ,
                         dllname = "myelo",initfunc = "parmsCraig16Q",
                         events=list(data=eventdat),method="lsoda",  #method made no diff
                         nout = 2, outnames = c("Qts","Cpts"))    )

D=data.frame(yout)
tail(D,2)
d=D%>%select(time,Q,Cp)%>%gather(key="Lab",value="Value",-time)
tc=function(sz) theme_classic(base_size=sz)
gx=xlab("Days")
sbb=theme(strip.background=element_blank())
cc=coord_cartesian(xlim=c(-2,85))
d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
ggsave("~/Results/myelo/Q.pdf",width=5, height=5)

