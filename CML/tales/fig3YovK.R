library(myelo)  
library(deSolve)
library(tidyverse)
glauchePars20
(ic=c(X=0,Y=1,Z=0))
(d=glauchePars20%>%group_by(id)%>%nest())
d$data[[1]]
hahnYovK<-function(Time, State, Pars) {
  with(as.list(c(Time, State, Pars)),{
    dX = -pxy*X + pyx*Y 
    dY =  pxy*X - pyx*Y + py*Y*(1-Y/Ky)   -  m*Y*Z 
    dZ =  rz    -   a*Z                   + pz*Y*Z/(Kz^2+Y^2)  
    list(c(dX,dY,dZ),c(ratio=2+log10(Y/Ky)))
  })
}

fsim=function(x) {
  ic[3]=x$rz/x$a
  ode(y = ic, times = seq(0,10*12,1), func = hahnYovK, parms = x)
}
d=d%>%mutate(out=map(data,fsim))
head(d$out[[1]])
as_tibble(d$out[[1]])%>%select(time,ratio)
d=d%>%mutate(D=map(out,function(x) as_tibble(x)%>%select(time,ratio)
                   %>%mutate(time=as.numeric(time)/12,ratio=as.numeric(ratio))))
dd=d%>%select(id,D)
dd=dd%>%unnest(cols=D)
dd$id=as_factor(dd$id)
tc=function(sz) theme_classic(base_size=sz)
gx=xlab("Years")
gy=ylab("2+log10(y/Ky)")
tms=scale_x_continuous(breaks=c(0,1,5,10))# times
sbb=theme(strip.background=element_blank())
dd%>%ggplot(aes(x=time,y=ratio))+facet_wrap(id~.,ncol=3)+geom_line(size=1)+gx+gy+tc(14)+sbb+tms 
ggsave("~/Results/twoCities/hahnYovK.pdf",width=6,height=7)
# ggsave("~/Results/twoCities/hahnYovK.png",width=4,height=6)
library(WriteXLS)
tb1=glauchePars20%>%select(-CessT,-TKI)%>%rename(Patient=id) #%>%head()
WriteXLS(list(hahnel=tb1), ExcelFileName="~/Results/twoCities/hahnelPars.xlsx")
