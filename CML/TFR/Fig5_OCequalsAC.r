## remake Fig5 with OC mrt assumed to be equal to AC, i.e. with HM mort assumed to be zero. 
# The idea here is that for CML patient deaths by HM, O is 20-50 fold bigger than E, so EAR = (O-E)/PY is ~O/PY.
# This avoids 2D surface fitting. Measured all-cause mortality rates can now be used directly across covid years,
# so there is no longer any reason to clip survival time strips at the end of 2019. 
graphics.off();rm(list=ls())
setwd("~/cml/mort")
library(tidyverse)
library(SEERaBomb)
source("setup.R")
load("~/data/mrt/mrtUSA.RData")#loads US mortality data
mrt[[1]]
str(mrt)
F=reshape2:::melt(mrt[[1]],value.name="rate")
M=reshape2:::melt(mrt[[2]],value.name="rate")
F$sex="Female"
M$sex="Male"
bh<-bind_rows(F,M)|>tibble()
names(bh)[1:2]=c("age","year")
bh=bh|>mutate(age=as.numeric(as.character(age))+0.5)|>filter(!is.na(age))
bh  #background hazard as a long tibble

load("~/data/seer25/cml25.RData") ### made by Fig2A_mkSEER.R  
d=d%>%filter(agedx<90) 
(d=d%>%filter(yrdx>=2005)) 
yrs=1975:2022; ages=0.5:125.5
Z=matrix(0,ncol=length(yrs),nrow=length(ages))
colnames(Z)=yrs
rownames(Z)=ages
head(Z)# all zeros initially
brks=c(0,0.5,1,2,3,4,5,6,8,10,12) 
(binS<-levels(cut(brks+0.1,breaks=c(brks,100))))
L=NULL  # initiate big list L
L$brks=brks
L$binS=binS
L$d=d
L

for (i in c("HM","OC"))
  for (bin in binS) 
    for (j in c("Male","Female")) { 
      # j="Male"; bin="(0,0.5]"; i="OC"
      (binIndx=getBinInfo(bin,binS)["index"])
      (bin=binS[binIndx])
      (LL=getBinInfo(bin,binS)["LL"])
      head(d<-L$d%>%filter(sex==j))
      d=d%>%mutate(survC=cut(surv,breaks=c(-1,brks,100),include.lowest = TRUE)) 
      d=d%>%filter(surv>LL) 
      d=d%>%mutate(py=getPY(surv,bin,binS,brks)) # getpy leaves zeros when surv end is left of LL
      d=d%>%filter(py>0)  #get rid of such rows upfront
      d=d%>%mutate(ageL=agedx+LL) 
      d$year=floor(d$yrdx+LL)
      PYin=d%>%select(py,ageL,year)
      if(dim(PYin)[1]==0) binMidPnt=LL else binMidPnt=LL+sum(PYin$py)/dim(PYin)[1]/2
      PYin=as.matrix(PYin)
      PY=Z+0 # add 0 to make sure not a lazy copy of just the pointer 
      fillPYM(PYin, PY)
      dO=d%>%filter(survC==bin) 
      (nDead=sum(dO$CODS==i)) 
      N=1e4 ## shrink by N big then multiply by it later, so each PY of 1 lands in one bin
      head(Od<-dO%>%mutate(py=(CODS==i)/N)%>%select(py,ageL,year))
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
      D$cause=i
      (D=D|>tibble())
      D=left_join(D,bh)
      if (i=="HM") D$Eus=0  
      if (i=="OC")  D=D|>mutate(Eus=PY*rate)
      D=D%>%rename(num=Obs,denom=PY)
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
(mids=sapply(L[["HM"]],getMid)) # same PY at risk for each death type
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
ggsave(file=paste0("outs/Fig5OCequalsAC.png"),height=3,width=4) 
