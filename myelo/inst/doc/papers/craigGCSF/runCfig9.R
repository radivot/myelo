library(tidyverse)  # figure 9 using C code
library(deSolve)
library(myelo)

(x0=c(Q=craigPars16[["Qss"]],Nr=craigPars16[["Nrss"]],N=craigPars16[["Nss"]],
      G1=craigPars16[["G1ss"]],G2=craigPars16[["G2ss"]],
      Tn=craigPars16[["tauNP"]]+craigPars16[["aNM"]],  #here time Tn is to reserves
      An=craigPars16[["ANss"]],Aq=craigPars16[["AQss"]],Cp=0,Cf=0,Cs1=0,Cs2=0,Gs=0,Ic=0,Ig=0))

(gtimes=as.numeric(t(outer(seq(0,80,14),4:13,"+"))))
n=length(gtimes)
(eventG=tibble(var=rep("Gs",n),
                   time=gtimes,
                   value=rep(craigPars16[["F300"]]*300e3/craigPars16[["Vd300"]],n),
                   # value=rep(craigPars16[["F750"]]*750e3/craigPars16[["Vd750"]],n),
                   method=rep("add",n)))

(ctimes=seq(0,80,14))
nc=length(ctimes)
dose=4*craigPars16[["BSA"]]*1e3# ug of chemo per injection (D=4 mg/m2)
(eventCspike=tibble(var=rep("Cp",nc),
                    time=ctimes,
                    value=rep(dose,nc),
                    method=rep("add",nc)))
(eventdat=as.data.frame(bind_rows(eventG,eventCspike)%>%arrange(time)))

#### First run using R code RHS
times <- seq(-12,85,by=.1) # 14 secs
system.time(yout <- dede(x0,times = times, func = craig16,	parms = craigPars16,
                         events=list(data=eventdat),method="lsodar"))

myPlot=function(yout,cc) {
  D=data.frame(yout)
  D%>%filter(time>-.1,time<1)
  head(D)
  tail(D)
  d=D%>%select(time:Cp,Cs1,ANC)%>%gather(key="Lab",value="Value",-time)%>%
    mutate(Lab=factor(Lab,levels=c("Q","Nr","N","G1","G2","Tn","An","Aq","Cp","Cs1","ANC")))
  tc=function(sz) theme_classic(base_size=sz)
  gx=xlab("Days")
  sbb=theme(strip.background=element_blank())
  g=d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
  # g=d%>%ggplot(aes(x=time,y=Value))+facet_wrap(Lab~.,ncol=2,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
  print(g)
}
cc=coord_cartesian(xlim=c(-3,85))#clips high errorbars
myPlot(yout,cc)

#now try redoing using C code RHS
(f=file.path(system.file(paste("libs",Sys.getenv("R_ARCH"),sep=""), package = "myelo"),
             paste("myelo",.Platform$dynlib.ext,sep="")))
dyn.load(f)
getLoadedDLLs()
getDLLRegisteredRoutines("myelo")
getDLLRegisteredRoutines("myelo")


system.time(yout <- dede(x0,times = times, func = "derivsCraig16",	parms = craigPars16,
                         dllname = "myelo",initfunc = "parmsCraig16",
                         events=list(data=eventdat),method="lsodar", # I worry about events here
                         nout = 1, outnames = c("ANC"))    )
myPlot(yout,cc)

# thetahat=craigPars16
# p0=10^thetahat