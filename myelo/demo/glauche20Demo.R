library(myelo)  
library(deSolve)
library(tidyverse)
glauchePars20
(ic=c(X=0,Y=1,Z=0))
(d=glauchePars20%>%group_by(id)%>%nest())
d$data[[1]]
glauche20<-function(Time, State, Pars) {
  with(as.list(c(Time, State, Pars)),{
    dX = -pxy*X + pyx*Y 
    dY =  pxy*X - pyx*Y + py*Y*(1-Y/Ky)   -  m*Y*Z 
    dZ =  rz    -   a*Z                   + pz*Y*Z/(Kz^2+Y^2)  
    list(c(dX,dY,dZ),c(ratio=2+log10(Y/Ky)))
    # list(c(dX,dY,dZ),c(ratio=2+log10(Y/(Y+2*(Ky-Y))))) #more correct, but also more clutter
  })
}

fsim=function(x) {
  ic[3]=x$rz/x$a
  ode(y = ic, times = seq(0,10*12,1), func = glauche20, parms = x)
}
d=d%>%mutate(out=map(data,fsim))
head(d$out[[1]])
for (i in 1:21) plot(d$out[[i]],which=c("ratio"),xlab="Months",ylab="2+log10(Y/Ky)",main=paste("Patient",i))  
as_tibble(d$out[[1]])%>%select(time,ratio)
d=d%>%mutate(D=map(out,function(x) as_tibble(x)%>%select(time,ratio)
                   %>%mutate(time=as.numeric(time),ratio=as.numeric(ratio))))
pdf("~/Results/CML/hahnel.pdf",width=4,height=4)
# png("~/Results/CML/hahnel.png") #only last plot shows up
# ?plot.deSolve
# par(mfrow=c(6,2))
for (i in 1:21) plot(d$out[[i]],which=c("ratio"),xlab="Months",ylab="2+log10(Y/Ky)",main=paste("Patient",i))  
dev.off()
dd=d%>%select(id,D)
dd=dd%>%unnest(cols=D)
dd$id=as_factor(dd$id)
tc=function(sz) theme_classic(base_size=sz)
gx=xlab("Months")
gy=ylab("2+log10(Y/Ky)")
sbb=theme(strip.background=element_blank())
dd%>%ggplot(aes(x=time,y=ratio))+facet_wrap(id~.,ncol=3)+geom_line(size=1)+gx+gy+tc(14)+sbb 
ggsave("~/Results/CML/hahnel2.png",width=4,height=6)

