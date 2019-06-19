############################ craig16 Model ##################
library(tidyverse)
library(deSolve)
library(myelo)
craigPars16
x0=c(Q=craigPars16[["Qss"]],Nr=craigPars16[["Nrss"]],N=craigPars16[["Nss"]],
     G1=craigPars16[["G1ss"]],G2=craigPars16[["G2ss"]],
     Tn=craigPars16[["tauNP"]]+craigPars16[["aNM"]],  #here time Tn is to reserves
     An=craigPars16[["ANss"]])
times <- seq(-50,200,by=1)
(eventdat=data.frame(var="G1",time=25,value =150,method="add"))
# (eventdat <- data.frame(var = c("N"),
#                         time = c(25) ,
#                         value = c(2),
#                         method = c("rep")))
yout <- dede(x0,times = times, func = craig16,	parms = craigPars16,
             events=list(data=eventdat),method="lsodar")
D=data.frame(yout)
head(D)
d=D%>%select(time:N)%>%gather(key="Cell",value="Counts",-time)%>%
  mutate(Cell=factor(Cell,levels=c("Q","Nr","N")))
tc=function(sz) theme_classic(base_size=sz)
gx=xlab("Days")
sbb=theme(strip.background=element_blank())
d%>%ggplot(aes(x=time,y=Counts))+facet_grid(Cell~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb
ggsave("~/Results/myelo/craig16noChemo.pdf",height=6,width=6.5)

