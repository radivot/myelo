# Tom Hahnel , Christoph Baldow , Artur C. Fassoni , Joelle Guilhot , Francois Guilhot , Susanne
# Saussele , Satu Mustjoki , Stefanie Jilg , Philipp J. Jost , Stephanie Dulucq , Francois-Xavier
# Mahon , Ingo Roeder , Ingmar Glauche
# Inferring immunological control mechanisms from TKI dose alterations in CML patients
# in http://dx.doi.org/10.1101/722546

# In Table S1 of the supplement of this paper we find: 
pt1 = c(0.01,0.05,13.90,2210.00,362.0, 85.90, 742834)
pt2 = c(0.00,0.05, 1.95,1940.00,508.0, 85.70, 109001)
pt3 = c(0.88,0.18, 6.48,   0.00, 4.29, 72.90,  33678)
pt4 = c(0.00,0.00, 2.31,53800.0, 9570, 95.20, 217698)
pt5 = c(0.00,0.01, 2.38,  176.0,28.70, 60.30, 312895)
pt6 = c(0.70,0.21, 2.07, 311000, 6480, 53.00,   2740)
pt7 = c(0.00,0.00, 1.86,   0.01,26.30,108.00,  11949)
pt8 = c(0.00,0.00, 1.86, 696.00,74.70, 84.70,  61731)
pt9 = c(0.00,0.07, 2.56, 325.00,85.70, 60.70, 739824)
pt10= c(0.00,0.04, 2.77,    473,115.0, 70.50, 559804)
pt11= c(0.00,0.00, 2.57,2820.00,318.0, 50.10, 878259)
pt12= c(0.04,0.00, 4.30,19800.0,462.0,141.00, 863636)
pt13= c(0.11,0.00, 3.23, 132000, 2120,145.00, 975440)
pt14= c(0.01,0.00, 2.13,  54000, 1430,165.00,1045239)
pt15= c(0.59,0.06, 3.62,  82600, 1300,131.00, 599823)
pt16= c(0.59,0.00,17.70,  57600, 1550, 60.70, 546945)
pt17= c(0.00,0.05, 2.36,  26800, 5700, 95.60, 996734)
pt18= c(0.00,0.05, 3.32,   0.00, 6.84, 81.30,1630508)
pt19= c(0.00,0.00, 1.86,  14.70, 1740, 92.60, 637927)
pt20= c(0.01,0.04, 2.38,   4420,  419, 99.00, 628117)
pt21= c(0.02,0.04, 7.01,   2530,  212, 69.30, 309978)
library(tidyverse)
pts=rbind(pt1,pt2,pt3,pt4,pt5,pt6,pt7,pt8,pt9,pt10,pt11,pt12,pt13,pt14,pt15,pt16,pt17,pt18,pt19,pt20,pt21)
colnames(pts)=c("pyx","pxy","TKI","pz","Kz","CessT","Y0")
pts=as.data.frame(pts)
pts$py=1.658 # cells/month
pts$Ky=1e6 # cells
pts$a=2 # 1/month (apoptosis rate constant)
pts$rz=200 # cells/month
pts$m=1e-4 # 1/(cell*month)
pts$id=rownames(pts)
(x=pts%>%arrange(pz))
# Classes
# A: 6  HM-like no stable subclinical steady state  (all 6 recurred), no bands in hahnel S5
A=c(2,3,7,9,18,19)  ;length(A)
# B  8  N-like very stable subclinical steady states 
B=c(6, 12,13,14,15,16, 20, 21);length(B)#table S2 (wide bands in S5),REB figure 2 (never form CML if start over) 
# C  7  HF-like labile stable subclinical steady states (4 of 7 recurred)
C=c(1,4,5,8,10,11,17); length(C)  # difference, narrow bands, could go either way
# pts$abc="A";pts$abc[B]="B";pts$abc[C]="C"
pts$grp="A_hm";pts$grp[B]="B_n";pts$grp[C]="C_hf"
glauchePars20=pts%>%select(id,grp,everything())
glauchePars20=glauchePars20%>%mutate(A=rz-a*py/m,B=pz*py/m,C=Kz^2*(rz-a*py/m) )
glauchePars20=glauchePars20%>%mutate(Uss=round((-B-sqrt(B^2-4*A*C))/(2*A)),
                                     Sss=round((-B+sqrt(B^2-4*A*C))/(2*A),1),
                                     lGap=round(log10(Uss)-log10(Sss),2) )
glauchePars20=glauchePars20%>%mutate(lGap=ifelse(is.nan(lGap),-1,lGap))%>%select(-A,-B,-C)
glauchePars20
glauchePars20=as_tibble(glauchePars20)
glauchePars20
glauchePars20%>%ggplot(aes(x=1:21,y=lGap,col=grp))+geom_point() # separable in diff of logs, 20 and 21 are more C like

save(glauchePars20,file="glauchePars20.rda")
# library(tools)
# showNonASCIIfile("inst/doc/glauchePars20.r") 
