############################ craig16 Model ##################
library(tidyverse)
library(deSolve)
library(myelo)
craigPars16
(x0=c(Q=craigPars16[["Qss"]],Nr=craigPars16[["Nrss"]],N=craigPars16[["Nss"]],
      G1=craigPars16[["G1ss"]],G2=craigPars16[["G2ss"]],
      Tn=craigPars16[["tauNP"]]+craigPars16[["aNM"]],  #here time Tn is to reserves
      An=craigPars16[["ANss"]],Aq=craigPars16[["AQss"]],Cp=0,Cf=0,Cs1=0,Cs2=0,Gs=0,Ic=0,Ig=0))
(eventdat=data.frame(var="Gs",
                     time=0,
                     value=craigPars16[["F750"]]*750e3/craigPars16[["Vd750"]],
                     method="add"))
times <- seq(-15,2,by=.01)
yout <- dede(x0,times = times, func = craig16,	parms = craigPars16,
             events=list(data=eventdat),method="lsodar")

plotQtoG2=function(yout,cc) {
  D=data.frame(yout)
  d=D%>%select(time:G2)%>%gather(key="Lab",value="Value",-time)%>%
    mutate(Lab=factor(Lab,levels=c("Q","Nr","N","G1","G2")))
  tc=function(sz) theme_classic(base_size=sz)
  gx=xlab("Days")
  sbb=theme(strip.background=element_blank())
  g=d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
  print(g)
}
cc=coord_cartesian(xlim=c(-.1,2))#clips high errorbars
plotQtoG2(yout,cc)
ggsave("~/Results/myelo/craig16fig5b.png",height=6,width=6.5)

(eventdat=data.frame(var="G1",
                     time=0,
                     value=750e3/craigPars16[["Vd750"]],
                     method="add"))
yout <- dede(x0,times = times, func = craig16,	parms = craigPars16,
             events=list(data=eventdat),method="lsodar")
plotQtoG2(yout,cc)
ggsave("~/Results/myelo/craig16fig5aSpike.png",height=6,width=6.5)

yout[,"G1"]=log10(yout[,"G1"])
plotQtoG2(yout,cc)
ggsave("~/Results/myelo/craig16fig5aSpikeLog.png",height=6,width=6.5)

Tinf=25/(24*60) #turn off at t=25 minutes of infusion
(eventdat=data.frame(var=rep("Ig",2),
                     time=c(0,Tinf), #turn off at t=Tinf days
                     value=c(750e3/(Tinf*craigPars16[["Vd750"]]),0),
                     method=rep("rep",2)))
yout <- dede(x0,times = times, func = craig16,	parms = craigPars16,
             events=list(data=eventdat),method="lsodar")
plotQtoG2(yout,cc)
ggsave("~/Results/myelo/craig16fig5aSquareWave.png",height=6,width=6.5)
yout[,"G1"]=log10(yout[,"G1"])
plotQtoG2(yout,cc)
ggsave("~/Results/myelo/craig16fig5aSquareWaveLog.png",height=6,width=6.5)


