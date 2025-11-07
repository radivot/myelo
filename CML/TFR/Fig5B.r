library(ggplot2)
library(dplyr)
library(survival)
library(survminer)
load("~/data/seer25/cml25.RData") ### made by Fig1A_mkSEER.R 
range(d$agedx) # 90 is 90+
(d=d|>filter(yrdx>=2005)) #all ages
d=d|>mutate(ageC=cut(agedx,seq(0,100,10),include.lowest = TRUE, right = FALSE))
(L=split.data.frame(d,d$ageC))
myplt=function(df){
  fit=survfit(Surv(surv, status) ~ sex, data = df)
  names(fit$strata) <- gsub("sex=", "", names(fit$strata))
  ggsurvplot(fit,legend.title="", title=paste0("Dx Age ",df$ageC[1]), data=d)+xlab("Years")
}
(Lg=lapply(L,myplt))
# Lg[[1]] #replots 1st age group

library(SEERaBomb)
load("~/data/mrt/mrtUSA.RData")#loads US mortality data mrt
(Lf=lapply(L,function(x) msd(x,mrt,brkst=c(0,5))|>foldD()))
Db=bind_rows(Lf,.id="ageGroup")|>filter(int=="(5,100]")
(D5=lapply(L,function(x) x|>filter(surv>5)))
Db$OC=sapply(D5,function(x) sum(x$CODS=="OC"))
Db$HM=sapply(D5,function(x) sum(x$CODS=="HM"))
Db=Db|>filter(ageGroup!="[90,100]")
Lf=lapply(L,function(x) msd(x,mrt,brkst=c(0,1))|>foldD())
(Db=bind_rows(Lf,.id="ageGroup")|>filter(int=="(1,100]"))
(D1=lapply(L,function(x) x|>filter(surv>1)))
Db$OC=sapply(D1,function(x) sum(x$CODS=="OC"))
(Db=Db|>filter(ageGroup!="[90,100]"))
(Db=Db|>mutate(EARoc = (OC - E)/PY, LLoc = EARoc - 1.96 * sqrt(OC)/PY, ULoc = EARoc + 1.96 * sqrt(OC)/PY, 
            RRoc = OC/E, rrLoc = qchisq(0.025,2 * OC)/(2 * E), rrUoc = qchisq(0.975, 2 * OC + 2)/(2*E)) |> filter(E > 0))
#   ageGroup int         O       E     PY     t     EAR       LL     UL    RR   rrL   rrU    OC   EARoc      LLoc    ULoc  RRoc rrLoc rrUoc
# 1 [0,10)   (1,100]     5   0.130   545.  4.65 0.00894 0.000895 0.0170 38.6  12.5  90.1      2 0.00343 -0.00165  0.00852 15.4   1.87 55.8 
# 2 [10,20)  (1,100]    31   2.83   3074.  4.96 0.00916 0.00561  0.0127 11.0   7.44 15.5     10 0.00233  0.000316 0.00435  3.53  1.69  6.50
# 3 [20,30)  (1,100]   117  12.6    8059.  4.56 0.0130  0.0103   0.0156  9.30  7.69 11.1     35 0.00278  0.00134  0.00422  2.78  1.94  3.87
# 4 [30,40)  (1,100]   217  32.7   13344.  4.42 0.0138  0.0117   0.0160  6.65  5.79  7.59    78 0.00340  0.00210  0.00470  2.39  1.89  2.98
# 5 [40,50)  (1,100]   396 102.    20158.  4.63 0.0146  0.0127   0.0165  3.90  3.52  4.30   161 0.00295  0.00171  0.00418  1.58  1.35  1.85
# 6 [50,60)  (1,100]   720 254.    24551.  4.35 0.0190  0.0169   0.0211  2.84  2.64  3.06   438 0.00751  0.00584  0.00919  1.73  1.57  1.90
# 7 [60,70)  (1,100]  1119 443.    21156.  3.88 0.0320  0.0289   0.0351  2.53  2.38  2.68   758 0.0149   0.0123   0.0174   1.71  1.59  1.84
# 8 [70,80)  (1,100]  1344 594.    12339.  3.30 0.0608  0.0550   0.0666  2.26  2.14  2.39   860 0.0216   0.0169   0.0262   1.45  1.35  1.55
# 9 [80,90)  (1,100]   996 495.     4286.  2.57 0.117   0.102    0.131   2.01  1.89  2.14   662 0.0389   0.0272   0.0507   1.34  1.24  1.44
source("setup.R")
EARbrks=c(0,0.01,0.05)
ccE=coord_cartesian(ylim=c(0,0.085))
ghp01=geom_hline(yintercept=c(0.01),col="gray")
gh0=geom_hline(yintercept=0,col="gray")
(Db$age=seq(5,85,by=10))
geEoc=geom_errorbar(aes(ymin=LLoc,ymax=ULoc),width=0.2)#for absolute risks
Db|>ggplot(aes(x=age,y=EARoc))+ghp01+geEoc+geom_point(size=0.7)+gyE+ ggtitle("Excess OC Risks at >1y after Dx")+gl+#ccE+
  scale_y_continuous(minor_breaks=NULL,breaks=EARbrks) +tc(13)+gh0+leg+xlab("Decade of Age at Dx")
ggsave(file=paste0("outs/Fig5B.pdf"),height=3,width=4) 

