## Fig1A_mkSEER25.R   makes binaries of SEER CMLincidence data for Fig1A.R and Fig5.R 
library(tidyverse)
# Go into SEER*stat, and for SEER 8, 12 and 17, using plus, make these case-listing columns
# Patient_ID
# Agerecodewithsingle_ages_and_90
# Sex
# Race_recode_White_Black_Other
# Year_of_diagnosis
# Site_recode_ICD_O_3_WHO_2008
# Histologic_Type_ICD_O_3
# Combined_Summary_Stage_2004
# Survival_months
# COD_to_site_recode
# Radiation_recode
# Chemotherapy_recode_yes_no_unk
nms=c("id","agedx","sex","race","yrdx","who","histo3","surv","COD","rad","chemo")
system.time(d8 <-read_csv("~/data/seer25/csvs/s8_25.txt", col_names=nms,skip=1)) 
system.time(d12<-read_csv("~/data/seer25/csvs/s12_25.txt",col_names=nms,skip=1))
system.time(d17<-read_csv("~/data/seer25/csvs/s17_25.txt",col_names=nms,skip=1)) # 2 secs (uses multi-cores)
d=bind_rows(d8,d12,d17) # need patient ID to make each row unique
(d=distinct(d))
d=d%>%mutate(histo3=8000+histo3)
d=d%>%mutate(yrdx=yrdx+1800,status=as.numeric(COD>0),surv=(surv+0.5)/12)
d=d%>%filter(surv<200)  # use 200 years to clip off surv = 9999 months = 833 years
d=d%>%mutate(sex=ifelse(sex==1,"Male","Female"))
d=d%>%mutate(sex=as_factor(sex))
d$seqnum=1 # didn't fetch from SEER*stat so stick first cancer in there (SEERaBomb wants the field present)
head(d) #have who and histo3 to map to cancer types
# devtools::install_github("radivot/SEERaBomb",subdir="SEERaBomb") # get 2019.5, not CRAN's 2019.2
library(SEERaBomb)
d=mapCancsW(d)   
mapCancsW  # starts with who (Site_recode_ICD_O_3_WHO_2008) and over-rides with histo3 defs    
d=d%>%rename(cancer=cancerW)
d%>%filter(cancer%in%c("CMLw")) #empty, good => histo3 covered all who==78 CMLs
head(d<-d%>%filter(cancer%in%c("CML"))) #so zoom in on just histo3 CML 
table(d$who)  # all who = 78 good too, nothing added by histo3 => both on same page
table(d$histo3,d$yrdx) # CML is all 9863 for yrdx up to 1993
table(d$COD) # below, not only 790 AML > 273 lung => misclassification, but also ALL & CLL > pancreatic
# from s8_25.sas           among CML patients
# 33 = "Pancreas"            62 deaths 
# 39 = "Lung and Bronchus"   273 
# 67 = "Hodgkin Lymphoma"      5
# 70 = "Non-Hodgkin Lymphoma"   107
# 73 = "Myeloma"                 31
# 74 = "Acute Lymphocytic Leukemia" 119
# 75 = "Chronic Lymphocytic Leukemia" 209  ***  high => many are really CML
# 76 = "Other Lymphocytic Leukemia"    29
# 77 = "Acute Myeloid Leukemia"       790 **** high => really CML    
# 80 = "Acute Monocytic Leukemia"      11
# 78 = "Chronic Myeloid Leukemia"      6297 expect to be high 
# 89 = "Other Myeloid/Monocytic Leukemia" 519
# 83 = "Other Acute Leukemia"           302
# 85 = "Aleukemic, Subleukemic and NOS"  1006 **** high => really CML    
#### so pool it all together as HM, and any excess HM deaths are then due to CML
mapCOD2HMOC=function(D){
  COD=D$COD #start with vec of integers. Map to a vec of Strings
  CODS=character(dim(D)[1]) #initialize string (S) version of COD to all ""
  CODS[(COD>=67)&(COD<=85)|(COD==89)]="HM"
  CODS[(COD<=66)|(COD==86)|(COD>=90)]="OC"
  CODS[COD==0]="alive"
  D$CODS=CODS 
  D
}
d=mapCOD2HMOC(d)
table(d$CODS)
# alive    HM    OC 
# 16279  9425  7837 
d=d%>%select(-seqnum,-COD,-who,-id,-race,-histo3,-rad,-chemo) # zap what is not needed in this paper
d
system.time(save(d,file="~/data/seer25/cml25.RData")) #0.03 secs  33,541 CML cases

