## Alpha and Betas in Hahnel et al. Cancer Research (2020).

Goal: Use BCR-ABL percentages above detection and before cessation to demo 
population level modeling of bi-exponential BCR-ABL 
decays often observed initially with TKI use. 

The model of interest is y=Ae<sup>αt</sup>+Be<sup>βt</sup>. When the first term dominates, the slope on a natural 
log scale is α and when the second term dominates the slope is β. On a log10 scale, these
values are log(10)=2.3 times smaller, as is the residual error of the fits. We will plot 
data and fits on a log10 scale, but we will perform fits on a natural log scale. Taking logs lets small y gain weight needed to estimate β. 


First, as in Hahnel et al, maxLik is used to fit the model.
```
rm(list=ls()) 
library(tidyverse)
library(myelo)
head(d<-hahnelFigS2)
(d=d%>%filter(UL==0)%>%select(-UL)) # use only non-zero measurements to keep it simple
names(d)=c("TIME","DV","ID")  #change names to NONMEM style
d=d%>%mutate(LDV=log10(DV),LNDV=log(DV)) 
library(maxLik)
loglik <- function(param) {
  lA <- log(param[1])
  alpha <- param[2]
  lB <- log(param[3])
  beta <- param[4]
  mu <- log(exp(lA+alpha*t)+exp(lB+beta*t))
  sigma <- param[5]
 sum(dnorm(x, mu, sigma, log=TRUE))
} #indirect passing of x and t via global is a weakness of maxLik
(d21=d%>%group_by(ID)%>%summarize(Tf=max(TIME)))
d21$Tpred=lapply(d21$Tf,function(x) seq(0,x))
L=NULL
for (i in 1:21) {   # this weakness drove my use here of a for loop
  print(i)
  d1=d%>%filter(ID==i)
  print(d1)
  t=d1$TIME
  x=d1$LNDV
  N=length(x)
  L[[i]]=summary(maxLik(loglik, start=c(A=100,alpha=-1,B=1,beta=-0.05,sigma=0.5)))  
}
(L=lapply(L,coef))
(P=lapply(L,function(x) x[,1]))
d21$P=P
d21
simY <- function(param,tpred,id) {
  A <- param[1]
  alpha <- param[2]
  B <- param[3]
  beta <- param[4]
  data.frame(t=tpred,y=log10(A*exp(alpha*tpred)+B*exp(beta*tpred)),ID=id)
}
(D=mapply(simY,P,d21$Tpred,1:21,SIMPLIFY = F))
library(data.table)
(D=data.table(bind_rows(D)))
names(D)=c("TIME","LDV","ID")
mkDF=function(x) data.frame(TIME=100,LDV=2.5,a=round(x[["alpha"]],3),b=round(x[["beta"]],3))
(DTx=bind_rows(lapply(P,mkDF)))
DTx$ID=1:21
DTxb=DTx
DTxb$LDV=1
d%>%ggplot(aes(x=TIME,y=LDV))+facet_wrap(ID~.,ncol=4)+geom_point(size=1)+
   geom_line(data=D)+geom_text(aes(label=a),data=DTx,parse=T)+geom_text(aes(label=b),data=DTxb,parse=T)
ggsave("../docs/alphaNbetaDataNfits.png",width=7,height=8)
```

![](../../../docs/alphaNbetaDataNfits.png)
which shows α and β values ~2.3-fold higher than log10 scale α and β slopes in Fig S2. 
Fits of patients 4, 10 and 18 could be better.

Using the R package bbmle, similar fits can be obtained as follows. 
```
library(bbmle)
dn=d%>%group_by(ID)%>%nest()
fitBi=function(d1)  {
  coef(summary(mle2(LNDV~dnorm(mean=log((exp(lA+alpha*TIME) + exp(lB+beta*TIME))),sd=sigma),
               method="Nelder-Mead", 
               start=list(lA=log(100),alpha=-1,lB=log(1),beta=-0.05,sigma=0.5),data=d1,
               control = list(maxit=500))))  
}
dn=dn%>%mutate(M=map(data,fitBi))
dn$M[[20]]
getPars=function(x){
  x=x[,1] 
  x[c(1,3)]=exp(x[c(1,3)])
  names(x)[c(1,3)]=c("A","B")
  x
}
dn=dn%>%mutate(P=map(M,getPars))

(D=mapply(simY,dn$P,d21$Tpred,1:21,SIMPLIFY = F))
(D=data.table(bind_rows(D)))
names(D)=c("TIME","LDV","ID")
(DTx=bind_rows(lapply(dn$P,mkDF)))
DTx$ID=1:21
DTxb=DTx
DTxb$LDV=1
d%>%ggplot(aes(x=TIME,y=LDV))+facet_wrap(ID~.,ncol=4)+geom_point(size=1)+
  geom_line(data=D)+geom_text(aes(label=a),data=DTx,parse=T)+geom_text(aes(label=b),data=DTxb,parse=T)
ggsave("../docs/alphaNbetaDataNfitsBB.png",width=7,height=8)

```

![](../../../docs/alphaNbetaDataNfitsBB.png)
which shows better fits of patients 4, 10 and 18.

Differences in parameter estimates can be drastic: in the table below, compare A and alpha of 
pts 1, 3, 6, and 7, and B and beta of pt 13. 
```
pML<-cbind(d21[,1],bind_rows(P))
pML$method="maxLik"
pBB<-cbind(d21[,1],bind_rows(dn$P))
pBB$method="bbmle"
bind_rows(pML,pBB)%>%arrange(ID)%>%filter(ID%in%c(1,3,6,7,13))
   ID            A      alpha           B         beta     sigma method
1   1 4.222862e+01 -0.9970453 0.016054450 -0.033776730 0.6272613 maxLik
2   1 8.070877e+04 -3.8916201 0.024939337 -0.041162706 0.7423625  bbmle
3   3 2.298123e+01 -1.0438443 0.469013963 -0.130820280 0.4787884 maxLik
4   3 4.144447e+02 -2.0406933 0.500721898 -0.133096927 0.4802431  bbmle
5   6 1.171979e+03 -2.3246318 0.027987216 -0.050265072 0.6637285 maxLik
6   6 4.024998e+05 -3.7792521 0.028233184 -0.050498802 0.6635207  bbmle
7   7 5.140817e-01 -0.1885761 0.005353326 -0.011377815 0.7660208 maxLik
8   7 1.855883e+03 -2.6750067 0.032174230 -0.033431420 0.9228490  bbmle
9  13 1.044945e+02 -0.5194102 0.001098639  0.003243757 1.6680958 maxLik
10 13 1.611426e+02 -0.5611876 0.051454394 -0.055396728 2.5671772  bbmle
```



