############################ ColijnO5b  ##################
library(deSolve)
library(myelo)
colijnPars05b
times <- seq(-200,5000,by=0.1)
yout <- dede(c(Q=colijnPars05b[["Qss"]],N=colijnPars05b[["Nss"]],R=colijnPars05b[["Rss"]],P=colijnPars05b[["Pss"]]),
             times = times, func = colijn05b,	parms = colijnPars05b)
plot(yout)
plot(yout,select="Q")
plot(yout,select="N")
plot(yout,select="R")
plot(yout,select="P")
(SS<-colijnPars05b[c("Qss","Nss","Rss","Pss")])
x0=as.vector(tail(yout,1))[2:5]
nms=colnames(yout)
names(x0)=nms[2:5]
(newSS=x0)
colijnPars05b[["Qss"]]=x0[1]
colijnPars05b[["Nss"]]=x0[2]
colijnPars05b[["Rss"]]=x0[3]
colijnPars05b[["Pss"]]=x0[4]

times <- seq(-200,2000,by=0.1)
yout <- dede(c(Q=colijnPars05b[["Qss"]],N=colijnPars05b[["Nss"]],R=colijnPars05b[["Rss"]],P=colijnPars05b[["Pss"]]),
             times = times, func = colijn05b,	parms = colijnPars05b)
plot(yout,select="Q")
plot(yout,select="N")
plot(yout,select="R")
plot(yout,select="P")

# kick N as before
(eventdat <- data.frame(var = c("N"),
                        time = c(25) ,
                        value = c(20),
                        method = c("rep")))

times <- seq(-200,200,by=0.1)
yout <- dede(c(Q=colijnPars05b[["Qss"]],N=colijnPars05b[["Nss"]],
               R=colijnPars05b[["Rss"]],P=colijnPars05b[["Pss"]]),
             times = times, func = colijn05b,	parms = colijnPars05b,
             events=list(data=eventdat),method="lsodar")
             
D=data.frame(yout)
library(tidyverse)
tc=function(sz) theme_classic(base_size=sz)
gy=ylab("Neutrophil Counts")
gx=xlab("Days")
cc=coord_cartesian(xlim=c(-3,50))#clips high errorbars
D%>%ggplot(aes(x=time,y=N))+geom_line(size=1)+gx+gy+tc(14)+cc
ggsave("~/Results/myelo/colijn05bJumpN.png",height=6,width=6.5)

