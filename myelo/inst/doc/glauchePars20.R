# Tom Hahnel , Christoph Baldow , Artur C. Fassoni , Joelle Guilhot , Francois Guilhot , Susanne
# Saussele , Satu Mustjoki , Stefanie Jilg , Philipp J. Jost , Stephanie Dulucq , Francois-Xavier
# Mahon , Ingo Roeder , Ingmar Glauche
# Inferring immunological control mechanisms from TKI dose alterations in CML patients
# in http://dx.doi.org/10.1101/722546

# In Table S1 of the supplement of this paper we find: 
pt1 = c(0.01,0.05,13.90,2210.00,362.0, 85.90)
pt2 = c(0.00,0.05, 1.95,1940.00,508.0, 85.70)
pt3 = c(0.88,0.18, 6.48,   0.00, 4.29, 72.90)
pt4 = c(0.00,0.00, 2.31,53800.0, 9570, 95.20)
pt5 = c(0.00,0.01, 2.38,  176.0,28.70, 60.30)
pt6 = c(0.70,0.21, 2.07, 311000, 6480, 53.00)
pt7 = c(0.00,0.00, 1.86,   0.01,26.30,108.00)
pt8 = c(0.00,0.00, 1.86, 696.00,74.70, 84.70)
pt9 = c(0.00,0.07, 2.56, 325.00,85.70, 60.70)
pt10= c(0.00,0.04, 2.77,    473,115.0, 70.50)
pt11= c(0.00,0.00, 2.57,2820.00,318.0, 50.10)
pt12= c(0.04,0.00, 4.30,19800.0,462.0,141.00)
pt13= c(0.11,0.00, 3.23, 132000, 2120,145.00)
pt14= c(0.01,0.00, 2.13,  54000, 1430,165.00)
pt15= c(0.59,0.06, 3.62,  82600, 1300,131.00)
pt16= c(0.59,0.00,17.70,  57600, 1550, 60.70)
pt17= c(0.00,0.05, 2.36,  26800, 5700, 95.60)
pt18= c(0.00,0.05, 3.32,   0.00, 6.84, 81.30)
pt19= c(0.00,0.00, 1.86,  14.70, 1740, 92.60)
pt20= c(0.01,0.04, 2.38,   4420,  419, 99.00)
pt21= c(0.02,0.04, 7.01,   2530,  212, 69.30)
library(tidyverse)
pts=rbind(pt1,pt2,pt3,pt4,pt5,pt6,pt7,pt8,pt9,pt10,pt11,pt12,pt13,pt14,pt15,pt16,pt17,pt18,pt19,pt20,pt21)
colnames(pts)=c("pyx","pxy","TKI","pz","Kz","CessT")
pts=as.data.frame(pts)
pts$pzOvKz=round(pts$pz/(pts$Kz)^2,2)
pts$py=1.658 # cells/month
pts$Ky=1e6 # cells
pts$a=2 # 1/month (apoptosis rate constant)
pts$rz=200 # cells/month
pts$m=1e-4 # 1/(cell*month)
pts$id=rownames(pts)
(x=pts%>%arrange(pzOvKz))
(x=pts%>%arrange(pz))
hist(x$pzOvKz,50)
pts$pzOvKz=NULL
glauchePars20=pts%>%select(id,everything())
save(glauchePars20,file="glauchePars20.rda")
as_tibble(glauchePars20%>%select(-TKI,-CessT))%>%print(n=21)
library(tools)
showNonASCIIfile("inst/doc/glauchePars20.r") 
