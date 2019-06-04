############################ ColijnO5a  ##################
library(deSolve)
library(myelo)
colijnPars05a
times <- seq(-200,5000,by=0.1)
yout <- dede(c(Q=colijnPars05a[["Qss"]],N=colijnPars05a[["Nss"]],R=colijnPars05a[["Rss"]],P=colijnPars05a[["Pss"]]),
             times = times, func = colijn05a,	parms = colijnPars05a)
plot(yout)
plot(yout,select="Q")
plot(yout,select="N")
plot(yout,select="R")
plot(yout,select="P")
(SS<-colijnPars05a[c("Qss","Nss","Rss","Pss")])
x0=as.vector(tail(yout,1))[2:5]
nms=colnames(yout)
names(x0)=nms[2:4]
newSS=x0
colijnPars05a[["Qss"]]=x0[1]
colijnPars05a[["Nss"]]=x0[2]
colijnPars05a[["Rss"]]=x0[3]
colijnPars05a[["Pss"]]=x0[4]

times <- seq(-200,2000,by=0.1)
yout <- dede(c(Q=colijnPars05a[["Qss"]],N=colijnPars05a[["Nss"]],R=colijnPars05a[["Rss"]],P=colijnPars05a[["Pss"]]),
             times = times, func = colijn05a,	parms = colijnPars05a)
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
yout <- dede(c(Q=colijnPars05a[["Qss"]],N=colijnPars05a[["Nss"]],
               R=colijnPars05a[["Rss"]],P=colijnPars05a[["Pss"]]),
             times = times, func = colijn05a,	parms = colijnPars05a,
             events=list(data=eventdat),method="lsodar")
             
D=data.frame(yout)
library(tidyverse)
tc=function(sz) theme_classic(base_size=sz)
gy=ylab("Neutrophil Counts")
gx=xlab("Days")
D%>%ggplot(aes(x=time,y=N))+geom_line(size=1)+gx+gy+tc(14) 
ggsave("~/Results/myelo/colijn05aJumpN.png",height=6,width=6.5)

