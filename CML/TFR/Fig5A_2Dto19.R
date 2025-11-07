## Fig5A_2Dto19.R  PY and deaths only up to end of 2019 to make sure excess OC are not from COVID
graphics.off();rm(list=ls())
setwd("~/cml/TFR")
library(tidyverse)
library(SEERaBomb)
library(mgcv)
source("setup.R")
system.time(load("~/data/seer25/cml25.RData")) ### made by Fig1A_mkSEER.R  33541
(d=d%>%filter(agedx<90)) # 32807 cases
(d=d%>%filter(yrdx>=2005))# 21149 cases 
(d=d%>%filter(yrdx<2020))# 17203 cases 
# (d=d%>%filter(agedx>20)) # 16706 cases
table(d$yrdx)
table(d$agedx)
d=d|>mutate(FUend=yrdx+surv)
d=d|>mutate(excess=FUend>2020)
d=d|>mutate(excessPY=ifelse(excess,FUend-2020,0))
d=d|>mutate(nsurv=ifelse(excess,surv-excessPY,surv))
d=d|>mutate(cod=ifelse(!excess,CODS,"alive"))
d
# # A tibble: 17,203 × 12
# agedx sex     yrdx   surv status cancer CODS  FUend excess excessPY  nsurv cod  
# <dbl> <fct>  <dbl>  <dbl>  <dbl> <fct>  <chr> <dbl> <lgl>     <dbl>  <dbl> <chr>
# 1    87 Male    2005 0.375       1 CML    OC    2005. FALSE     0     0.375  OC   
# 2    87 Female  2009 0.292       1 CML    OC    2009. FALSE     0     0.292  OC   
# 3    81 Female  2005 5.38        1 CML    HM    2010. FALSE     0     5.38   HM   
# 4    89 Female  2006 0.125       1 CML    OC    2006. FALSE     0     0.125  OC   
# 5    88 Male    2008 1.79        1 CML    HM    2010. FALSE     0     1.79   HM   
# 6    87 Male    2015 0.125       1 CML    HM    2015. FALSE     0     0.125  HM   
# 7    71 Female  2011 6.04        1 CML    OC    2017. FALSE     0     6.04   OC   
# 8    89 Male    2017 3.38        1 CML    OC    2020. TRUE      0.375 3      alive
# 9    89 Female  2012 0.0417      1 CML    OC    2012. FALSE     0     0.0417 OC   
# 10   54 Female  2013 9.62        0 CML    alive 2023. TRUE      2.62  7      alive

# yrs=2005:2019
# ages=0.5:70.5
# PYM=matrix(0,ncol=length(yrs),nrow=length(ages))
# colnames(PYM)=yrs
# rownames(PYM)=ages
# (PYin=structure(c(3.5, 11.25,5.2, 51.5, 58.5,0.75, 2005, 2007,2016),.Dim = c(3L,3L), 
#                 .Dimnames=list(c("1","2","3"),c("py", "ageL", "year"))))
# fillPYM(PYin, PYM)
yrs=1975:2019; ages=0.5:125.5
# yrs=2005:2019; ages=0.5:125.5
Z=matrix(0,ncol=length(yrs),nrow=length(ages))
colnames(Z)=yrs
rownames(Z)=ages
head(Z)# all zeros initially
# survCut=0
brks=c(0,0.5,1,2,3,4,5,6,8,10,12) 
(binS<-levels(cut(brks+0.1,breaks=c(brks,100))))
L=NULL  # initiate big list L
L$brks=brks
L$binS=binS
#these two files are made by Fig5_3mkSplines.R
load("~/data/seer25/Gh.RData")
load("~/data/seer25/Go.RData")
(L$G[["HM"]]=Gh)
(L$G[["OC"]]=Go)
L$d=d
L

for (i in c("HM","OC"))
  for (bin in binS) 
    for (j in c("Male","Female")) { 
      # j="Male"; bin="(0,0.5]"; i="HM"
      (binIndx=getBinInfo(bin,binS)["index"])
      (bin=binS[binIndx])
      (LL=getBinInfo(bin,binS)["LL"])
      head(d<-L$d%>%filter(sex==j))
      d=d%>%mutate(survC=cut(nsurv,breaks=c(-1,brks,100),include.lowest = TRUE)) 
      d=d%>%filter(nsurv>LL) 
      d=d%>%mutate(py=getPY(nsurv,bin,binS,brks)) # getpy leaves zeros when surv end is left of LL
      d=d%>%filter(py>0)  #get rid of such rows upfront
      d=d%>%mutate(ageL=agedx+LL) 
      d$year=floor(d$yrdx+LL)
      PYin=d%>%select(py,ageL,year)
      if(dim(PYin)[1]==0) binMidPnt=LL else binMidPnt=LL+sum(PYin$py)/dim(PYin)[1]/2
      PYin=as.matrix(PYin)
      PY=Z+0 # add 0 to make sure not a lazy copy of just the pointer 
      fillPYM(PYin, PY)
      dO=d%>%filter(survC==bin) 
      # (nDead=sum(dO$CODS==i)) 
      (nDead=sum(dO$cod==i)) 
      N=1e4 ## shrink by N big then multiply by it later, so each PY of 1 lands in one bin
      # head(Od<-dO%>%mutate(py=(CODS==i)/N)%>%select(py,ageL,year))
      head(Od<-dO%>%mutate(py=(cod==i)/N)%>%select(py,ageL,year))
      head(Oin<-as.matrix(Od))
      O=Z+0
      fillPYM(Oin, O)
      O=O*N
      head(dO<-reshape2:::melt(O,value.name="Obs"))
      head(dP<-reshape2:::melt(PY,value.name="PY"))
      head(dd<-left_join(dO,dP))
      names(dd)[1:2]<-c("age","year")
      (D=dd%>%filter(PY>0))
      D$sex=j
      D=D%>%rename(num=Obs,denom=PY)
      D$Eus=as.numeric(exp(predict(L$G[[i]],D)))
      D$cause=i
      # tibble(D)
      D=D%>%filter(age<90,age>20) 
      # sum(D$num) 
      L[[i]][[bin]][[j]]$mid=binMidPnt
      L[[i]][[bin]][[j]]$D=D
    }

str(L[["OC"]])
str(L[["HM"]])
getMid=function(x) mean(x$Male$mid,x$Female$mid)
sapply(L[["HM"]],getMid)
(mids=sapply(L[["OC"]],getMid)) # same PY at risk for each death type
poolMF=function(x) bind_rows(x$Male$D,x$Female$D)
(Doc=lapply(L[["OC"]],poolMF)) 
(Dhm=lapply(L[["HM"]],poolMF)) 
sumM=function(x) x%>%summarize(O=sum(num),PY=sum(denom),E=sum(Eus),cause=toupper(cause[1]))
(D=tibble(t=mids,int=binS,hm=Dhm,oc=Doc))
D=D%>%mutate(hm1=map(hm,sumM))
D=D%>%mutate(oc1=map(oc,sumM))
D=D%>%select(-hm,-oc)
Dh=D%>%select(-oc1)
Do=D%>%select(-hm1)
(Dh=Dh%>%unnest(hm1))
(Do=Do%>%unnest(oc1))
D=bind_rows(Do,Dh)
(D=D%>%mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY))
D=D%>%mutate(cause=if_else(cause=="OC","Other Cause","Heme Malignancy"))
D=D%>%mutate(cause=as_factor(cause)) # to get Other Cause first

leg=theme(legend.margin=margin(0,0,0,0),legend.title=element_blank(),
          legend.position="top") #redfine leg without size increase
EARbrks=c(0,0.01,0.05,0.1)
ccE=coord_cartesian(ylim=c(0,0.08))
ghp01=geom_hline(yintercept=c(0.01),col="gray")
gh0=geom_hline(yintercept=0)

D|>ggplot(aes(x=t,y=EAR,col = cause))+ghp01+geE+geom_point(size=0.7)+gyE +gl+ccE+
    scale_y_continuous(minor_breaks=NULL,breaks=EARbrks) +tc(13)+gh0+leg+gx+
   scale_color_manual(values = c("Other Cause" = "black", "Heme Malignancy" = "gray60"))
ggsave(file=paste0("outs/Fig5A_2Dto19.png"),height=3,width=4) 

