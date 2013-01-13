# 1x/wk fails, 2x/wk and greater works
library(myelo)
(f=file.path(system.file(paste("libs",Sys.getenv("R_ARCH"),sep=""), package = "myelo"),
					paste("myelo",.Platform$dynlib.ext,sep="")))
dyn.load(f)

(ff=file.path(system.file("demo", package = "myelo"),
					"SphaseAgentFuncs.R"))
# The model used here can be viewed as a barebones RL of just the malignant cells, only tracking the marrow
source(ff)

# fix S phase rates for 8 hour S-phase and
# let 0X => factor of 10 growth in 28 days 
# 10 =exp (28*24/tau) => tau = 292 
MDS=c(LI=8/24,Ts=8,Tc=24)
graphics.off()
windows(width=12,height=6)
par(mfcol=c(2,4),mar=c(5.1,4.1,2.1,2.1),cex.lab=1.6,cex.axis=1.6) 

# now fix Teff based on 1X => factor of 10 growth in 56 days
# Here Teff = effective time in cell per dose


TeffL=c(30,40,50,65)
TeffL=c(1,2,3,10)
TeffL=c(5,8,12,15)
TeffL=c(4,5,6,7)
TeffL=c(2,3,4,5)
#Teff=.10 #
drug=2

for (i in 1:length(TeffL))
{
  Teff=TeffL[i]
  
# define dosing schedules
# un per week
  
  upw <- data.frame(var = "drug",  
      time = sort(c(24*c(0,7,14,21),c(24*c(0,7,14,21)+Teff))),
      value = c(drug,0),  
      method = "rep"  )
  upw
# deux per week 
# 24 hrs apart
  dpw24 <- data.frame(var = "drug",  
      time = sort(c(24*c(0,1,7,8,14,15,21,22),c(24*c(0,1,7,8,14,15,21,22)+Teff))),
      value = c(drug,0),  
      method = "rep"  )
  dpw24
  if (Teff>=24) dpw24=dpw24[-c(2:3,6:7,10:11,14:15),]
  dpw24
  
  dpw48 <- data.frame(var = "drug",  
      time = sort(c(24*c(0,2,7,9,14,16,21,23),c(24*c(0,2,7,9,14,16,21,23)+Teff))),
      value = c(drug,0),  
      method = "rep")
  dpw48
  if (Teff>=48) dpw48=dpw48[-c(2:3,6:7,10:11,14:15),]
  dpw48
  
  
# 96 hrs apart
  dpw96 <- data.frame(var = "drug",  
      time = sort(c(24*c(0,4,7,11,14,18,21,25),c(24*c(0,4,7,11,14,18,21,25)+Teff))),
      value = c(drug,0),  
      method = "rep" )
  dpw96
  
# trois per week
  tpw4848 <- data.frame(var = "drug",  
      time = sort(c(24*c(0,2,4,7,9,11,14,16,18,21,23,25),c(24*c(0,2,4,7,9,11,14,16,18,21,23,25)+Teff))),
      value = c(drug,0),  
      method = "rep" )
  tpw4848
  
  tpw4848
  if (Teff>=48) tpw4848=tpw4848[-c(2:5,8:11,14:17,20:23),]
  tpw4848
  
  MDACC <- data.frame(var = "drug",  
      time = sort(c(24*(0:4),24*(0:4)+Teff)),
      value = c(drug,0),  
      method = "rep")
  MDACC
  
  if (Teff>=24) MDACC=MDACC[-c(2:9),]
  MDACC
  
  
  
#debug(mkg)  
  upwM=mkg(MDS,upw,Teff)
  dpw24M=mkg(MDS,dpw24,Teff)
  dpw48M=mkg(MDS,dpw48,Teff)
  dpw96M=mkg(MDS,dpw96,Teff)
  tpwM=mkg(MDS,tpw4848,Teff)
  
  lims=range(c(upwM$totals,dpw24M$totals,dpw48M$totals,dpw96M$totals,tpwM$totals))
#  lims=c(1,3e12)
  plot(upwM$days ,upwM$totals,xlab="Days",col="red",ylab="Leukemia Cells",log="y",type="l",ylim=lims
      ,main=paste("Teff = ",Teff))
  lines(dpw24M$days,dpw24M$totals,lty=1,col="orange")
  lines(dpw48M$days,dpw48M$totals,lty=1,col="green")
  lines(dpw96M$days,dpw96M$totals,lty=1,col="blue")
  lines(tpwM$days,tpwM$totals,lty=1,col="violet")
#legend("topleft",c("24/144","48/120","96/72","48/48/72"),lty=1,bty="n",col=c("red","green","blue","violet"))
  if (i==2)  legend("bottomleft",c("168","24/144","48/120","96/72","48/48/72"),lty=1,bty="n",
        col=c("red","orange","green","blue","violet")) # this plot cropped is Fig. 10 in the R01 
  
  #  debug(mkTmCrs)  
  mkTmCrs(MDS,Teff)
#  title(main=paste("Teff = ",Teff,"h  Ts = 8 h   Tc = 292 h"))
  
  
}

dyn.unload(f)
