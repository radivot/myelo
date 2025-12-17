library(myelo)  
library(deSolve)
library(tidyverse)
glauchePars20
(ic=c(x=0,y=1,z=0))
(d=glauchePars20%>%group_by(id)%>%nest())
d$data[[1]]
hahnXYZ<-function(Time, State, Pars) {
  with(as.list(c(Time, State, Pars)),{
    dX = -pxy*x + pyx*y 
    dY =  pxy*x - pyx*y + py*y*(1-y/Ky)   -  m*y*z 
    dZ =  rz    -   a*z                   + pz*y*z/(Kz^2+y^2)  
    list(c(dX,dY,dZ))
  })
}

fsim=function(x) {
  ic[3]=x$rz/x$a
  ode(y = ic, times = seq(0,10*12,1), func = hahnXYZ, parms = x)
}
d=d%>%mutate(out=map(data,fsim))
head(d$out[[1]])

(dd=as_tibble(d$out[[1]])%>%gather(key=Cell,value=CellNumber,-time))

d=d%>%mutate(D=map(out,function(x) as_tibble(x)%>%gather(key=Cell,value=CellNumber,-time)
                 %>%mutate(time=as.numeric(time)/12,CellNumber=as.numeric(CellNumber))))
dd=d%>%select(id,D)
dd=dd%>%unnest(cols=D)
dd$id=as_factor(dd$id)
tc=function(sz) theme_classic(base_size=sz)
gx=xlab("Years")
gy=ylab("Number of Cells")
# sy=scale_y_log10()
sy=scale_y_log10(breaks=c(1,1e3,1e6),labels=c(1,quote(paste(10^3)),quote(paste(10^6))))
ltp=theme(legend.position="top",legend.key.height=unit(.7,'lines'))
ltb=theme(legend.margin=margin(0,0,0,0),legend.title=element_blank())
sbb=theme(strip.background=element_blank())
tms=scale_x_continuous(breaks=c(0,1,5,10))# times
# ybks=scale_y_continuous(breaks=c(1,1e3,1e6),labels=c(1,1000,quote(paste(10^6))))
dd%>%ggplot(aes(x=time,y=CellNumber,col=Cell))+facet_wrap(id~.,ncol=3)+
  geom_line(size=1)+gx+gy+tc(14)+sbb+sy +ltp+ltb+tms
ggsave("~/Results/twoCities/hahnXYZ.pdf",width=6,height=7)
# ggsave("~/Results/twoCities/hahnXYZ.png",width=6,height=7)

