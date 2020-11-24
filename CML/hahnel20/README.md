# CML models


## Model-based inference and classification of immunological control mechanisms from TKI cessation and dose reduction in CML patients. Hahnel et al, Cancer Research (2020).


This model captures quiescent (x) and dividing (y) CML cells interacting with anti-CML immune cells (z). 
![](../../docs/glauche.png)

The differential equations of this model are:
![](../../docs/hahnelDEQs.png)

The following code runs the model for 10 years (120 months) from an initial condition of y(0)=1. 

```
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
  })
}

fsim=function(x) {
  ic[3]=x$rz/x$a
  ode(y = ic, times = seq(0,10*12,1), func = glauche20, parms = x)
}
d=d%>%mutate(out=map(data,fsim))
head(d$out[[1]])
d=d%>%mutate(D=map(out,function(x) as_tibble(x)%>%select(time,ratio)
                   %>%mutate(time=as.numeric(time),ratio=as.numeric(ratio))))
dd=d%>%select(id,D)
dd=dd%>%unnest(cols=D)
dd$id=as_factor(dd$id)
tc=function(sz) theme_classic(base_size=sz)
gx=xlab("Months")
gy=ylab("2+log10(Y/Ky)")
sbb=theme(strip.background=element_blank())
dd%>%ggplot(aes(x=time,y=ratio))+facet_grid(id~.)+geom_line(size=1)+gx+gy+tc(14)+sbb 
ggsave("~/Results/CML/hahnel.png",width=4,height=12)
```

Running this code generates the following plot
![](../../docs/hahnel.png)



