# use COPASI to save parmar sbml as XPPAUT *.ODE and do rest by hand + find/replace
library(myelo)   # load definition of function parmar19
library(deSolve)
library(rodeoExt) #provides rbind for deSolve class matrices

ic=c(NTBI=5.2e-11,Tf=1.50934e-08,FeDuo=3.85469e-07,FeBM=4.45824e-07,FeRest=7.91257e-06,FeLiver=1.51048e-06,
     Fe1Tf=1.05653e-08,FeSplee=1.52679e-07,EPO=3.15471e-15,Hepcidi=2.9912e-11,FeRBC=1.817e-05,Fe2Tf=2.46513e-08) 
(parameters=parmarPars19)

out=ode(y = ic, times = seq(-30,0,1), func = parmar19, parms = parameters) 
(N=length(ic))
n=dim(out)[1]
X0=out[n,2:(N+1)]
names(X0)<-names(ic)
parameters["ksHepci"]=0
out2   <- ode(y=X0, times=seq(1, 365, by = 1), func = parmar19, parms = parameters)
out=rbind(out,out2)
head(out)
graphics.off()
quartz(width=8,height=7)
vars=c("NTBI_c","Tf_c","FeDuo_c","FeBM_c","FeRBC_c","FeLiver_c","FeSplee_c","FeRest_c","Hepcidi_c")
plot(out,which=vars,xlab="Days",ylab="Concentration (M)") #screen capture to png

###########
library(tidyverse)
ltb=theme(legend.margin=margin(0,0,0,0),legend.title=element_blank())
ltp=theme(legend.position="top",legend.direction="horizontal")
tc=function(sz) theme_classic(base_size=sz)
sbb=theme(strip.background=element_blank())
gy=ylab("Concentration in M")
gx=xlab("Days")
D=data.frame(out)
head(D)
D=D%>%select(time,NTBI_c:Hepcidi_c)
d=D%>%gather(key=variable,value=Concentration,-time)
d%>%ggplot(aes(x=time,y=Concentration))+facet_wrap(~variable,scales="free",nrow=3)+geom_line(size=1)+gx+gy+tc(14)+ltb+ltp+sbb
ggsave("~/Results/myelo/parmar19.png",height=6,width=6.5)

library(CoRC)
unloadAllModels()
(L=loadExamples())
runTC(model=L[[1]])  #runs fine
path="~/ccf/jarek/grants/msb/iron/parmar19sup/cps/"
(m0=loadModel(paste0(path,"IronMousePV3.cps")))
runTC(model=m0) # no results
(m1=loadModel(paste0(path,"IronMousePV3_Hemochromatosis.cps")))
runTC(model=m1) # no results

# saveSBML("~/Results/myelo/brus.sbml",level=2,version=4,model=L[[1]],overwrite=T)
# saveSBML("~/Results/myelo/m0.sbml",level=2,version=4,model=m0,overwrite=T)
