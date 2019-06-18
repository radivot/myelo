############################ craig15 Model ##################
library(tidyverse)
library(deSolve)
library(myelo)
craigPars15
(oldSS<-craigPars15[c("Qss","Nss")])# Q in 1e6, N in 1e9 
x0=c(Q=craigPars15[["Qss"]],Nr=craigPars15[["Nrss"]],N=craigPars15[["Nss"]],
     Gs=0,G=craigPars15[["Gss"]],
     Tn=craigPars15[["tauNP"]]+craigPars15[["aNM"]]+craigPars15[["tauNr"]],
     Cp=0,Cf=0,Cs1=0,Cs2=0,Aq=craigPars15[["AQss"]],An=craigPars15[["ANss"]])
x0
times <- seq(-50,200,by=1)
yout <- dede(x0,times = times, func = craig15,	parms = craigPars15)

D=data.frame(yout)
head(D)
d=D%>%select(time:N)%>%gather(key="Cell",value="Counts",-time)%>%
  mutate(Cell=factor(Cell,levels=c("Q","Nr","N")))
tc=function(sz) theme_classic(base_size=sz)
gx=xlab("Days")
sbb=theme(strip.background=element_blank())
d%>%ggplot(aes(x=time,y=Counts))+facet_grid(Cell~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb
ggsave("~/Results/myelo/craig15ss.pdf",height=6,width=6.5)

