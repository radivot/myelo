## Fig5_3mkSplines.R   Step 3 of 3 before running Fig5.R   
# make surface splines of SEER HM and OC mortality data made in Step 2
setwd("~/ccf/CMLR/seer25")
graphics.off();rm(list=ls())#clear plots and environment 
library(mgcv)
library(rgl)
library(tidyverse)
load("~/data/seer25/d2.RData") #made by mkSEERmrts.R, gets dh and do

## Heme Malignancy (HM) cause of death, mortality vs age x year x sex
D=dh
(D=D%>%filter(age>20))
D$sex=as_factor(D$sex)
(Df=D%>%filter(year<2020)) # data for fitting
summary(Gh<-gam(num ~ sex+s(age,year,by=sex)+ti(age,year)+offset(log(denom)),family=poisson(),data=Df)) 
save(Gh,file="~/data/seer25/Gh.RData")  

# rest below is just to see how well the 2D HM splines fit the 2D data
D$E=exp(predict(Gh,D))
head(D<-D%>%mutate(Eincid=E/denom))
head(D<-D%>%mutate(incid=num/denom))
head(D<-D%>%mutate(incid=if_else(incid==0,0.0001,incid)))
with(D,plot3d(year,age,log10(incid),col=ifelse(sex=="Female","red","blue"),xlab="Year",ylab="Age",zlab="log10(M)",alpha=1,size=4))
(Ages=seq(min(D$age),max(D$age)))
(Years=seq(min(D$year),max(D$year)))
head(nD<-expand.grid(Ages,Years))
names(nD)<-c("age","year")
nD$denom=1
nDf=nD%>%mutate(sex="Female")
nDm=nD%>%mutate(sex="Male")
nD=bind_rows(nDf,nDm)
head(nD)
dim(nD)  # = 99*35
nD$E=exp(predict(Gh,nD))
M=reshape2::acast(nD%>%filter(sex=="Female")%>%select(year,age,E), year~age, value.var="E")
surface3d(Years,Ages,log10(M),col="red",alpha=0.5) # M for Matrix
M=reshape2::acast(nD%>%filter(sex=="Male")%>%select(year,age,E), year~age, value.var="E")
surface3d(Years,Ages,log10(M),col="blue",alpha=0.5)
clear3d(type="lights")
light3d(theta = 0, phi = 70) #expand X window and take snapshot into powerpoint (from there export as pdf)


##  other cause (OC) mortality
D=do
(D=D%>%filter(age>20))
D$sex=as_factor(D$sex)
(Df=D%>%filter(year<2020)) # data for fitting
summary(Go<-gam(num ~ sex+s(age,year,by=sex)+ti(age,year)+offset(log(denom)),family=poisson(),data=Df))  
save(Go,file="~/data/seer25/Go.RData")

# rest below is just to see how well the OC splines fit the data
D$E=exp(predict(Go,D))
head(D<-D%>%mutate(Eincid=E/denom))
head(D<-D%>%mutate(incid=num/denom))
head(D<-D%>%mutate(incid=if_else(incid==0,0.0001,incid)))
with(D,plot3d(year,age,log10(incid),col=ifelse(sex=="Female","red","blue"),xlab="Year",ylab="Age",zlab="log10(M)",alpha=1,size=4))
(Ages=seq(min(D$age),max(D$age)))
(Years=seq(min(D$year),max(D$year)))
head(nD<-expand.grid(Ages,Years))
names(nD)<-c("age","year")
nD$denom=1
nDf=nD%>%mutate(sex="Female")
nDm=nD%>%mutate(sex="Male")
nD=bind_rows(nDf,nDm)
head(nD)
dim(nD)  # = 99*35
nD$E=exp(predict(Go,nD))
M=reshape2::acast(nD%>%filter(sex=="Female")%>%select(year,age,E), year~age, value.var="E")
surface3d(Years,Ages,log10(M),col="red",alpha=0.5) # M for Matrix
M=reshape2::acast(nD%>%filter(sex=="Male")%>%select(year,age,E), year~age, value.var="E")
surface3d(Years,Ages,log10(M),col="blue",alpha=0.5)
clear3d(type="lights")
light3d(theta = 0, phi = 70) 
