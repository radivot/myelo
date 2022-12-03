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
# higher res vals from Tom Hahnel's email 
pts$pyx=c(0.0126, 0.00033, 0.878, 0.000271, 0.000181, 0.702, 0.00106, 
          0.000215, 0.000378, 0.000635, 0.000302, 0.0379, 0.107, 0.00985, 
          0.585, 0.586, 0.00122, 0.000457, 4.54e-05, 0.00895, 0.0224)
pts$pxy=c(0.0483, 0.0452, 0.184, 1.68e-06, 0.00615, 0.205, 0.000406, 
          4.26e-05, 0.0732, 0.0447, 0.00292, 0.00276, 4.56e-05, 1.1e-06, 
          0.056, 3.22e-05, 0.0452, 0.051, 2.86e-06, 0.0387, 0.0353)
pts$pz=c(2210, 1940, 1.5e-07, 53800, 176, 311000, 0.00585, 696, 325, 
         473, 2820, 19800, 132000, 54000, 82600, 57600, 26800, 6.57e-07, 
         14.7, 4420, 2530)
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
# dX = -pxy*X + pyx*Y#                              Xss=(pyx/pxy)*Yss
# dY =  pxy*X - pyx*Y + py*Y*(1-Y/Ky) -  m*Y*Z#   Zss=(py/m)(1-Yss/Ky)
# dZ =  rz    -   a*Z +  pz*Y*Z/(Kz^2+Y^2)
# rz -a*(py/m)(1-Yss/Ky) + pz*Yss*(py/m)(1-Yss/Ky)/(Kz^2+Yss^2)
# assuming 1-Yss/Ky is roughly 1 (i.e. Yss << Ky), yields
# rz -a*(py/m) + pz*Yss*(py/m)/(Kz^2+Yss^2)
# (rz -a*(py/m))(Kz^2+Yss^2) + pz*Yss*(py/m)
# (rz -a*(py/m))*Kz^2 + pz*Yss*(py/m) +(rz -a*(py/m))*Yss^2
glauchePars20=glauchePars20%>%mutate(A=rz-a*py/m,B=pz*py/m,C=Kz^2*(rz-a*py/m) )
glauchePars20=glauchePars20%>%mutate(U2=round((-B-sqrt(B^2-4*A*C))/(2*A)),
                                     S2=round((-B+sqrt(B^2-4*A*C))/(2*A)),
                                     lGap2=round(log10(U2)-log10(S2),2) )
glauchePars20=glauchePars20%>%mutate(lGap2=ifelse(is.nan(lGap2),-Inf,lGap2))%>%select(-A,-B,-C)
glauchePars20
glauchePars20=as_tibble(glauchePars20)
glauchePars20
glauchePars20%>%ggplot(aes(x=1:21,y=lGap2,col=grp))+geom_point() # separable in diff of logs, 20 and 21 are more C like



# NOT assuming 1-Yss/Ky = 1 yields cubic polynomial in Yss 
# rz -a*(py/m)(1-Yss/Ky) + pz*Yss*(py/m)(1-Yss/Ky)/(Kz^2+Yss^2)
# rz*(Kz^2+Yss^2) -a*(py/m)(1-Yss/Ky)*(Kz^2+Yss^2) + pz*Yss*(py/m)(1-Yss/Ky)
# Yss^0 * (rz - a*(py/m))*Kz^2 +  
# Yss^1 * ((a*(py/m)*Kz^2)/Ky  +  pz*(py/m))  +  
# Yss^2 * (rz - a*(py/m)  - pz*(py/m)/Ky)  +  
# Yss^3 * (a*(py/m)/Ky)  
library(RConics)
glauchePars20=glauchePars20%>%mutate(A=a*(py/m)/Ky,
                                     B=rz - a*(py/m)  - pz*(py/m)/Ky,
                                     C=(a*(py/m)*Kz^2)/Ky  +  pz*(py/m),
                                     D=(rz - a*(py/m))*Kz^2)

fcub1=function(x) {
  inp=c(x$A,x$B,x$C,x$D)
  rt=cubic(inp)[1]
  if(is.complex(rt)) return(NaN) else round(return(rt))
}


fcub2=function(x) {
  inp=c(x$A,x$B,x$C,x$D)
  rt=cubic(inp)[2]
  if(is.complex(rt)) return(NaN) else round(return(rt))
}

fcub3=function(x) {
  inp=c(x$A,x$B,x$C,x$D)
  rt=cubic(inp)[3]
  if(is.complex(rt)) return(NaN) else round(return(rt))
}

(dn=glauchePars20%>%group_by(id)%>%nest())
dn$data[1]
dn=dn%>%mutate(S3=round(map_dbl(data,fcub1))) 
dn=dn%>%mutate(U3=round(map_dbl(data,fcub2))) 
dn=dn%>%mutate(U3=round(map_dbl(data,fcub3))) 
dn=dn%>%mutate(lGap3=round(log10(U3)-log10(S3),2)) 
glauchePars20=dn%>%unnest(data)%>%select(-A,-B,-C,-D)
glauchePars20=glauchePars20%>%mutate(lGap3=ifelse(is.nan(lGap3),-Inf,lGap3))
glauchePars20
glauchePars20%>%ggplot(aes(x=1:21,y=lGap3,col=grp))+geom_point() # same look as lGap2
glauchePars20%>%select(S2,S3,U2,U3)%>%print(n=21) #biggest diff is  pt6 (highest Uss not negligible relative to 1e6)


# Fig. S7 immune window
# what about Z balance point of death = activation
# dZ =  rz    -   a*Z +  pz*Y*Z/(Kz^2+Y^2)
#i.e.
# pz*Y/(Kz^2+Y^2) = a
# pz*Y= a(Kz^2+Y^2)
# a*Y^2 - pz*Y + a*Kz^2

glauchePars20=glauchePars20%>%mutate(A=a,B=-pz,C=Kz^2*a) 
glauchePars20=glauchePars20%>%mutate(Ymin=round((-B-sqrt(B^2-4*A*C))/(2*A)),
                                     Ymax=round((-B+sqrt(B^2-4*A*C))/(2*A)),
                                     lGap=round(log10(Ymax)-log10(Ymin),2) )
glauchePars20=glauchePars20%>%mutate(lGap=ifelse(is.nan(lGap),-Inf,lGap))%>%select(-A,-B,-C)
glauchePars20%>%ggplot(aes(x=1:21,y=lGap,col=grp))+geom_point() # also same look as lGap2
glauchePars20%>%select(Ymin,S2,S3,Ymax,U2,U3)%>%print(n=21)
glauchePars20%>%select(Ymin,S2,S3,Ymax,U2,U3)%>%as.data.frame

#Initial conditions: assume quasi-equilibrium
glauchePars20=glauchePars20%>%mutate(X0=(pyx/pxy)*Y0,humpVal=pz*Y0/(Kz^2+Y0^2),Z0=rz/(a - humpVal))
glauchePars20%>%select(pz,Kz,Y0,X0,a,humpVal,Z0)%>%arrange(desc(humpVal))%>%as.data.frame
# patient 6 is in the window where activation is high, so equilibrium is not possible, so set to pre-CML steady state
glauchePars20=glauchePars20%>%mutate(Z0=ifelse(Z0>0,Z0,rz/a))
glauchePars20=glauchePars20%>%mutate(Z0=signif(Z0,3),X0=signif(X0,3))

glauchePars20%>%select(pz,Kz,Y0,X0,a,humpVal,Z0)%>%arrange(desc(humpVal))%>%as.data.frame

glauchePars20%>%select(pz,Kz,Y0,X0,a,Z0)%>%arrange(desc(X0))%>%as.data.frame

glauchePars20=glauchePars20%>%select(-humpVal)
glauchePars20



glauchePars20=glauchePars20%>%select(id:CessT,py:m,X0,Y0,Z0,Ymin,S2,S3,Ymax,U2,U3,lGap,lGap2,lGap3)
glauchePars20=glauchePars20%>%select(-lGap,-lGap2)%>%as.data.frame
glauchePars20

save(glauchePars20,file="glauchePars20.rda")


# library(tools)
# showNonASCIIfile("inst/doc/glauchePars20.r") 
