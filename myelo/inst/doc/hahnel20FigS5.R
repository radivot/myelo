# The program WebPlotDigitizer (see https://automeris.io/WebPlotDigitizer/) was
# used to convert  TFR (post stopping TKI) data in Figure S5 of Hahnel et al into CSV files. 
### Note: you shouldn't need to run this It is here just to document how
### I went from  csv outputs of WebPlotDigitizer to an R dataframe. 

myfun=function(j) {
  m=FALSE;u=FALSE
  if(file.exists(f<-paste0("CML/dresden/s5csv/pt",j,"m.csv"))) {
    (dm=read.csv(f,header=F))
    names(dm)=c("Months","Prct")
    dm$Pt=j
    dm$UL=0
    m=TRUE
    # j=2
  }
  if(file.exists(f<-paste0("CML/dresden/s5csv/pt",j,"u.csv"))) {
    (du=read.csv(f,header=F))
    names(du)=c("Months","Prct")
    du$Pt=j
    du$UL=1
    u=TRUE
  }
  if(m&u)  d=rbind(dm,du) else
    if (m) d=dm else d=du
    d
}
L=lapply(1:21,myfun)
library(tidyverse)
hahnelFigS5=bind_rows(L)
save(hahnelFigS5,file="CML/dresden/data/hahnelFigS5.rda")


