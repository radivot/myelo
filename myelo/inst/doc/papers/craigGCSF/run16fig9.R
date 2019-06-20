############################ craig16 Model ##################
library(tidyverse)
library(deSolve)
library(myelo)
craigPars16

# craigPars16[["kapmin"]] =0.0052359  # value 1 in Table 2 = Fig 11a => WORSE! higher ANC peaks, lower HSC Nadirs

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
(delt=round(1/24,2)) # 1-hour chemo infusions
dose=4*craigPars16[["BSA"]]*1e3# ug of chemo per injection (D=4 mg/m2)
infusionRate=dose/delt # this feeds straight into dCp
(eventCon=tibble(var=rep("Ic",nc),
                     time=ctimes,
                     value=rep(infusionRate,nc),
                     method=rep("rep",nc)))

(eventCoff=tibble(var=rep("Ic",nc),
                      time=ctimes+delt,
                      value=rep(0,nc),
                      method=rep("rep",nc)))
(eventdat=as.data.frame(bind_rows(eventG,eventCon,eventCoff)%>%arrange(time)))

# times <- seq(-12,30,by=.01)
# times <- seq(-12,30,by=.01) # 15 secs
times <- seq(-15,85,by=.01) # 42 secs
# yout <- dede(x0,times = times, func = craig16,	parms = craigPars16,
#              events=list(data=eventdat),method="lsodar")
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
# ggsave("~/Results/myelo/craig16fig9kpminLow.png",height=9,width=6.5) #WORSE! higher ANC peaks, lower HSC Nadirs
ggsave("~/Results/myelo/craig16fig9.png",height=9,width=6.5)


(eventCspike=tibble(var=rep("Cp",nc),
                    time=ctimes,
                    value=rep(dose,nc),
                    method=rep("add",nc)))

(eventdat=as.data.frame(bind_rows(eventG,eventCspike)%>%arrange(time)))

times <- seq(-12,85,by=.01) # 40 secs
times <- seq(-12,85,by=.1) # 18 secs
# times <- seq(-12,85,by=1) # 18 secs
system.time(yout <- dede(x0,times = times, func = craig16,	parms = craigPars16,
             events=list(data=eventdat),method="lsodar"))
myPlot(yout,cc)
ggsave("~/Results/myelo/craig16fig9spike.png",height=9,width=6.5)


(eventG=tibble(var=rep("Gs",n),
               time=gtimes,
               value=rep(craigPars16[["F750"]]*750e3/craigPars16[["Vd750"]],n),
               method=rep("add",n)))
(eventdat=as.data.frame(bind_rows(eventG,eventCspike)%>%arrange(time)))

yout <- dede(x0,times = times, func = craig16,	parms = craigPars16,
             events=list(data=eventdat),method="lsodar")

myPlot(yout,cc)
ggsave("~/Results/myelo/craig16fig9spikeChemoG750.png",height=9,width=6.5)

