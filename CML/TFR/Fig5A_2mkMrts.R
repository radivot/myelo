## Fig5_2mkMrts.R   Step 2 of 3 before running Fig5.R   
# make SEER OC and HM mortality data objects used in Step 3
setwd("~/ccf/CMLR/seer25")
rm(list=ls())#clear environment 
library(tidyverse) 
nms=c("sex","age","COD","year","num") #SEER: Mortality - All COD, Aggregated Total U.S. (1969-2023)
head(b <-read_csv("~/data/seer25/csvs/cod4.txt",skip=1, col_names=nms)) 
# 0=All Causes of Death
# 1="    Lymphoma"
# 2="    Myeloma"
# 3="    Leukemia"
nms=c("sex","age","year","denom") #SEER:  Populations - Total U.S. (1969-2023)
head(p <-read_csv("~/data/seer25/csvs/pops.txt",skip=1, col_names=nms))
d=left_join(b,p)
(d=d|>mutate(rate=num/denom,year=year+1975,sex=ifelse(sex==0,"Male","Female")))
#load Ages, a list of two sexes of age-year matrices of PY-weighted 5-year bin age-group midpoints 
load("~/data/mrt/Ages.RData") ## Use Ages made in Fig5_mkAges.R
getAge=function(age,year,sex) Ages[[sex]][[as.character(age),as.character(year)]]
(d=d%>%mutate(age=pmap_dbl(list(age,year,sex),getAge))) #update ages to set up Poisson regressions
(da=d%>%filter(COD==0)%>%mutate(COD="AC")%>%select(-rate))
(dh=d%>%filter(COD>0)%>%mutate(COD="HM")%>%select(-rate))
(dh=dh%>%group_by(COD,sex,age,year,denom)%>%summarize(num=sum(num),.groups="drop"))
(da=da[,names(dh)]%>%arrange(sex,age,year))
(do=da%>%mutate(COD="OC",num=num-dh$num))
save(do,dh,file="~/data/seer25/d2.RData")  ## Used in Fig5_3mkSplines.R

