# The program WebPlotDigitizer (see https://automeris.io/WebPlotDigitizer/) was
# used to convert Figure S2 data in Hahnel et al into CSV files. For each
# patient (pt) there is at least measured (m) data and for all but 3, upper
# limits (u) of what they might be when left censored.

### Note: you shouldn't need to run this It is here just to document how
### I went from  csv outputs of WebPlotDigitizer to an R dataframe. 

myfun=function(j) {
  (dm=read.csv(paste0("CML/dresden/csv/pt",j,"m.csv"),header=F))
  names(dm)=c("Months","Prct")
  dm$Pt=j
  dm$UL=0
  # j=2
  if(file.exists(f<-paste0("CML/dresden/csv/pt",j,"u.csv"))) {
    (du=read.csv(f,header=F))
    names(du)=c("Months","Prct")
    du$Pt=j
    du$UL=1
    dm=rbind(dm,du)
  }
  dm
}
L=lapply(1:21,myfun)
library(tidyverse)
hahnelFigS2=bind_rows(L)
save(hahnelFigS2,file="CML/dresden/data/hahnelFigS2.rda")


