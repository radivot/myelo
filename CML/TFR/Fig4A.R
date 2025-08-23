# Fig4A.R
graphics.off();rm(list=ls())
setwd("~/cml/TFR") 
library(tidyverse)
library(pomp)
library(subplex)
library(myelo)
source("setup.R")
### get parameters from myelo
head(d <- glauchePars20) 
(pars=d[3,]%>%select(pyx:Z0,-CessT))
(nms=names(pars))
(pars=unlist(pars))
(CessT=unlist(d[3,"CessT"])) #72.9/12  = 6.075 years 

# get data from myelo
head(d1<-hahnelFigS2)
head(d2<-hahnelFigS5)
head(dh<-bind_rows(d1,d2),2)
head(dh<-dh%>%mutate(R=Prct/100,time=Months)%>%arrange(time))
head(dh<-dh%>%filter(Pt==3))
head(dh<-dh%>%select(UL,R,time))

# Utki is unit step drop of doses to 0 at/by CessT
(covs0=data.frame(time= c(-1,CessT, CessT+1, 121),
                  Utki = c(1,1, 0,0))) 

rinit=Csnippet("X = X0;
                Y = Y0;
                Z = Z0;")

skel=vectorfield(Csnippet("DX = -pxy*X + pyx*Y;
                           DY = pxy*X - pyx*Y + py*Y*(1-Y/Ky) -  m*Y*Z - TKI*Utki*Y;
                           DZ = rz + pz*Y*Z/(pow(Kz,2.0)+pow(Y,2.0))  -  az*Z;"))   

hahn=dh%>%pomp(times="time",t0=0, 
               rinit=rinit, 
               skeleton=skel,
               covar = covariate_table(covs0, order="linear", times = "time"),
               paramnames=nms,
               statenames=c("X","Y","Z"),
               params=pars)


head(D1<-hahn%>%trajectory(times=0:120,format="data.frame")%>%mutate(R=Y/(Y+2.0*(pars["Ky"]-Y)))) 
D1%>%ggplot(aes(x=time/12,y=R)) + theme_grey(base_size=14) +
  geom_hline(yintercept=0.001,col="gray",linetype="dotted")+
  geom_hline(yintercept=0.0001,col="gray",linetype="dashed")+
  geom_line()+scale_y_log10()+expand_limits(y=c(1e-6,1))+xlab("Years")+ ylab("")+  
  geom_vline(xintercept=CessT/12,col="gray")+
  geom_point(aes(x=time/12,y=R,shape=factor(UL)),data=dh%>%mutate(R=ifelse(UL==1,0,R)))+
  theme(legend.position="none") +ylab("BCR::ABL1/ABL1")+scale_shape_manual(values = c(19, 2))

# get data ready for Poisson regression 
dh=dh%>%mutate(Denom=ifelse(UL>0,3/R,1e5))  # R = 3/read)
(dx=dh%>%filter(UL>0))
(GlobDenom=mean(dx$Denom)) # 89028.47
(dh=dh%>%mutate(Num=ifelse(UL>0,0,ceiling(R*GlobDenom))))  #use ceiling to avoid zeros 
(dh=dh%>%mutate(Denom=ifelse(UL>0,3/R,GlobDenom)))  # R = 3/read)
(CessT=CessT -dh$time[1])  #  69.98464/12 = 5.83 years
dh$time=dh$time-dh$time[1]
dh
Iend=95 # time at which immunity is lost
Irestrt=105 # time at which immunity is regained, at 104 it returns to immuno-setpoint (doesn't make it to watershed)
drop=0  # full stop between CessT and  CessT + 1
(covs0=data.frame(time= c(-1,   CessT, CessT+1,    Iend-1,  Iend, Irestrt,Irestrt+1, 121),
                  Utki = c(1,       1,    drop,     drop,  drop,   drop,    drop,   drop),
                  Upz =   c(1,      1,       1,        1,     0,      0,       1,     1), 
                  Denom=  c(GlobDenom, NA,  NA,       NA,    NA,     NA,      NA,     NA)))  
(covs=bind_rows(covs0,dh%>%select(time,Denom))%>%arrange(time)%>%fill(Utki,Upz,Denom))# fill NAs 

dmeas <- Csnippet("lik = dpois_raw(round(Num),round(Denom*Y/(Y+2*(Ky-Y))),give_log);
                    lik = (give_log) ? lik : exp(lik);")   

# fixed params
(R0=dh$R[1])    # ba=Y/(Y+2(K-Y)) => baY+ ba2K-ba2Y=Y => Y(1+ba)=ba2K 
(pars["Y0"]=R0*2e6/(1+R0)) #27953.25
pars["Kz"]=632.5  
pars["pz"]=4373.5
#fitted params
pars["pxy"]=0.11
pars["pyx"]=2.91
pars["TKI"]=14.03

rinit=Csnippet("X = pyx*Y0/pxy;
                Y = Y0;
                Z = rz/(az-pz*Y0/(pow(Kz,2.0)+pow(Y0,2.0)));")

skel=vectorfield(Csnippet("DX = -pxy*X + pyx*Y;
                           DY = pxy*X - pyx*Y + py*Y*(1-Y/Ky) -  m*Y*Z - TKI*Utki*Y;
                           DZ = rz + pz*Upz*Y*Z/(pow(Kz,2.0)+pow(Y,2.0))  -  az*Z;"))   

dh=dh%>%select(-Denom) # drop Denom from dataset else name conflict with covs
hahn=dh%>%pomp(times="time",t0=0, 
               rinit=rinit, 
               skeleton=skel,
               dmeasure=dmeas,
               covar = covariate_table(covs, order="linear", times = "time"),
               partrans=parameter_trans(log=c("pyx","pxy","TKI")),
               paramnames=nms,
               statenames=c("X","Y","Z"),
               params=pars) 

## check plot of initial parameter estimates to see they are in the ballpark
head(D0<-hahn%>%trajectory(times=0:120,format="data.frame")%>%mutate(R=Y/(Y+2.0*(pars["Ky"]-Y)))) 
D0%>%ggplot(aes(x=time/12,y=R)) + 
  theme_grey(base_size=14) +
  geom_line()+scale_y_log10()+expand_limits(y=c(1e-6,1))+xlab("Years")+ ylab("")+ 
  geom_point(aes(x=time/12,y=R,shape=factor(UL)),data=dh%>%mutate(R=ifelse(UL==1,0,R)))+
  theme(legend.position="none") +ylab("BCR::ABL1/ABL1")+scale_shape_manual(values = c(19, 2))

##### optimize estimates of pyx, pxy and TKI
(lpars=c(log(pars[c(1:3)])))
estS<-c("pyx","pxy","TKI")
(ofun=hahn%>%traj_objfun(est=estS))
(fit=subplex(lpars,fn=ofun,control=list(abstol=1e-5,maxit=2e6),hessian=TRUE))  # -2LL  209.7345
tpars=pars
tpars[c(1:3)]=c(exp(fit$par[c(1:3)]))
hahn=hahn%>%pomp(params=tpars)
head(D1<-hahn%>%trajectory(times=0:120,format="data.frame")%>%mutate(R=Y/(Y+2.0*(tpars["Ky"]-Y)))) 
#          X         Y        Z time .id           R
# 1 688039.1 27953.253 108.4821    0   1 0.014174742  # get X0 and Z0 from here
tpars["X0"]=D1[1,"X"]
tpars["Z0"]=D1[1,"Z"]
dput(c(tpars,CessT=CessT))
# c(pyx = 2.891981305438, pxy = 0.117493735003685, TKI = 13.3427401229456, 
#   pz = 4373.5, Kz = 632.5, py = 1.658, Ky = 1e+06, az = 2, rz = 200, 
#   m = 1e-04, X0 = 688039.114941631, Y0 = 27953.2531178926, Z0 = 108.482082597056, 
#   CessT = 69.9846394984326)  # these values go into Table 1
bks=c(0.01,0.1,1)
sy=scale_y_log10(breaks=bks,labels=bks)
D1%>%ggplot(aes(x=time/12,y=100*R))+gl+tc(14)+gx+
  geom_hline(yintercept=0.1,col="gray",linetype="dotted")+
  geom_hline(yintercept=0.01,col="gray",linetype="dashed")+
  geom_vline(xintercept=CessT/12,col="gray")+
  geom_vline(xintercept=(Iend-1)/12,col="gray")+
  geom_vline(xintercept=(Irestrt)/12,col="gray")+
  expand_limits(y=c(1e-4,100))+ #xlab("Years")+
  # scale_y_log10(breaks=c(0.00001,0.001,0.1))+expand_limits(y=c(1e-6,1))+xlab("Years")+ ylab("BCR-ABL/ABL")+
  geom_point(aes(x=time/12,y=100*R,shape=factor(UL)),data=dh%>%mutate(R=ifelse(UL==1,0,R)))+
  # theme(legend.position="none") +ylab("BCR::ABL1/ABL1")+scale_shape_manual(values = c(19, 2))
  theme(legend.position="none") + gyIS + sy+ scale_shape_manual(values = c(19, 2))
ggsave("outs/Fig4A.pdf",width=4,height=4.5)
(ic=as.numeric(D1[round(CessT)+1,nms<-c("X","Y","Z")])) #IC in Fig4B is at CessT (first row is t=0, so +1)
names(ic)<-nms
dput(ic)  # c(X = 771.745489104474, Y = 6.26564295379477, Z = 103.729321625289)
