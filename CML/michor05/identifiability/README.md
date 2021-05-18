#  FME Identifiability Analysis 
The goal here is to see how well parameters can be recovered from bi-exponential model simulations. 

First, log(BCRABL/BCR) is normal, as
nlmer fit residuals show a flat variance-mean relationship 
(the main characteristic of normality).   
```
rm(list=ls())
library(tidyverse)
library(myelo)
library(lme4)
head(d<-michor)
names(d)=c("TIME","DV","ID")  #change names to NONMEM style
head(d<-d%>%mutate(LNDV=log(DV)))
(startvec <- c(lA = log(0.995), alpha = -0.05, lB=log(0.005), beta = -0.005))
nform <- ~ log(exp(lA+alpha*input) + exp(lB+beta*input))
nfun <- deriv(nform, namevec=c("lA", "alpha", "lB", "beta"),
              function.arg=c("input","lA", "alpha", "lB", "beta"))
(M<-nlmer(LNDV ~ nfun(TIME, lA, alpha, lB, beta) ~ (lA|ID) + (alpha|ID) +(lB|ID) + (beta|ID),
          data=d,start = startvec,
          nAGQ = 0L,
          control = nlmerControl(tolPwrss = 1e-4)) )
head(d<-d%>%mutate(EY=predict(M),res=LNDV-EY,var=res^2,sd=sqrt(var)))
d%>%ggplot(aes(x=TIME,y=LNDV))+facet_wrap(ID~.,ncol=4)+geom_point(size=1)+geom_line(aes(y=EY)) 
# ggsave("outs/nlmerFits.pdf",width=6,height=30) #same as in myelo/CML/michor05/alphaBeta 
d%>%ggplot(aes(x=EY,y=var))+geom_point() #+scale_y_log10() #no trend, so log(Y) makes it normal
ggsave("outs/nlmerResids.png",width=6,height=30)
sqrt(mean(d$var)) # sd=0.6
```
![](outs/nlmerResids.png)


The deSolve C coded bi-exponential model (see adjacent compile folder) is loaded and simulated as follows. 
```
library(deSolve)
library(FME)

pars=c(Ke=0.05, Kpc=0.005,Ac0=0.995,Ap0=0.005)
(f=file.path(system.file(paste("libs",
    Sys.getenv("R_ARCH"),sep=""),package = "myelo"),
    paste("myelo",.Platform$dynlib.ext,sep="")))
dyn.load(f)
biExp <- function (pars) {
  y0=with(as.list(pars),c(Ac=Ac0,Ap=Ap0))
  out=ode(y=y0,times=seq(0,360,1),func="derivsBiExp",
          dllname = "myelo",initfunc = "parmsBiExp",
          parms=pars[1:2],
          nout = 1, outnames = c("ratio"))
  D=as.data.frame(out)
  D%>%mutate(lrat=log(ratio))
}

head(D<-biExp(pars))
head(D<-D%>%mutate(lrat=rnorm(dim(D)[1],mean=lrat,sd=0.6),sd=0.6))

tc=function(sz) theme_classic(base_size=sz)
gx=xlab("Days")
D%>%ggplot(aes(x=time,y=lrat))+geom_point()+tc(14)+gx
ggsave("outs/simData.png",width=4,height=4)
```
![](outs/simData.png)


Now define cost as the penalty of a lack of fit, assuming normality.
```
head(DataLogRat<-D%>%select(time, lrat,sd))
# we narrowed D above to make lrat the only choice 
biCost <- function (pars) {    #for the output 
  out = biExp(pars)            #picked out from out 
  return(modCost(model = out, obs = DataLogRat, err = "sd")) #here
}
x=biCost(pars)
x$model
gy=ylab("Weighted Residuals")
x$residuals%>%ggplot(aes(x,res))+geom_point()+xlab("Days")+gy
ggsave("outs/resids.png",width=4,height=4)
```
![](outs/resids.png)

The following chunk shows the impact of a 30% increase in Ke 
```
par(fig=c(0,1,0,1))
ref  <- biExp(pars)
pert <- biExp(pars*c(1.3,1,1,1))
plot(ref$time,ref$lrat,type="l",lwd=2,xlab="Days",ylab="lrat",
     main="local sensitivity, parameter Ke")
lines(pert$time,pert$lrat,lwd=2,lty=2)

arrseq <- seq(50,150,20)#c(10,30,50,70,90)
arrows(ref$time[arrseq],ref$lrat[arrseq],
       ref$time[arrseq],pert$lrat[arrseq], length= 0.075)
legend("bottomleft",c("Ke=0.05","Ke=0.065"),lty=c(1,2))
par(new=TRUE)
par(fig=c(0.5,0.99,0.5,0.95))

# note, denom below is negative, as is num, so this
ss <- (pert$lrat-ref$lrat)/pert$lrat # is positive
plot(ref$time,ss,type="l",lwd=2,
     xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
points(ref$time[arrseq],ss[arrseq])
title("Sensitivity functions ",line=0.5,cex.main=1)
par(fig=c(0,1,0,1))
# save as KeSensTimeCrs.png by hand from gui
```
![](outs/KeSensTimeCrs.png)
Arrows in the main plot match circles in the upper right. 

Sensitivities for all parameters are created below. 
```
Sfun <- sensFun(biCost, pars)
summary(Sfun)
plot(Sfun, which = c("lrat"), xlab="time", lwd = 2,legpos="topright")
# save as sensTimeCourses.png via gui
```
![](outs/sensTimeCourses.png)
Defined as output changes/values, sensitivity numerator and denominator minus signs cancel. 


The following creates pairwise phase-plane-like sensitivity plots. 
```
pairs(Sfun, which = c("lrat"), col = c("blue"))  
```
![](outs/pairsSfun.png)
For example, in row 1 col 2, starting at the origin, Ke sensitivity jumps first and Kpc sensitivity  later.
Lower diags are correlation coefficients. 

Next we look at collinearity
```
ident <- collin(Sfun)
head(ident, n = 20)
#    Ke Kpc Ac0 Ap0 N collinearity
# 1   1   1   0   0 2          1.1
# 2   1   0   1   0 2          1.9
# 3   1   0   0   1 2          1.2
# 4   0   1   1   0 2          1.0
# 5   0   1   0   1 2          4.5
# 6   0   0   1   1 2          1.1
# 7   1   1   1   0 3          2.0  #pick this
# 8   1   1   0   1 3          5.3
# 9   1   0   1   1 3          2.0  #or this
# 10  0   1   1   1 3          4.6
# 11  1   1   1   1 4          5.5
plot(ident, log = "y") #collinearity.png
```
![](outs/collinearity.png)
This suggests fixing the slope or intercept of the beta exponential.


Perturbing initial estimates 2-fold, we get back true values
```
biCost2 <- function(lpars) #nest to insure positive
  biCost(c(exp(lpars), n = 900)) # parameters
Pars <- pars[1:4] * 2  #double to make it work to find the opt
Fit <- modFit(f = biCost2, p = log(Pars))
exp(coef(Fit))
deviance(Fit)

ini   <- biExp(pars = c(Pars, n = 900))
final <- biExp(pars = c(exp(coef(Fit)), n = 900))
plot(D$time,D$lrat, xlab = "time", ylab = "lrat")
lines(ini$time, ini$lrat, lty = 2)
lines(final$time, final$lrat)
legend("topright", c("data", "initial", "fitted"),
       lty = c(NA,2,1), pch = c(1, NA, NA)) #=> initFitted.png
```
![](outs/initFitted.png)
which shows how much initial estimates lie away from final estimates. 


The following MCMC code takes time to run, so we save the result. 
```    
# MCMC <- modMCMC(f = biCost2, p = Fit$par,niter=1e5)
# save(MCMC, file="data/mcmc.Rdata")
```

Now check convergence
```
load("data/mcmc.Rdata")
options(width = 50)
MCMC$pars <- exp(MCMC$pars)
summary(MCMC)
par(mar=c(4, 4, 3, 1) + .1)
plot(MCMC, Full = TRUE) #converge.png
```
![](outs/converge.png)
Ac0 struggles due to 
saturation of Y/(Y+2) with Y=Ac+Ap. 

Now look at the correlations
```
pairs(MCMC,nsample=1000,cex.labels=1.4,cex=0.7)
```
![](outs/pairsMCMC.png)
Ac0 multimodality is consistent with convergence plateaus. Strong slope-intercept correlation
of the beta exponential is also apparent. 

Beta intercept burial in the alpha intercept is inferred from this
```
sR=sensRange(func=biExp,parms=pars,parInput=MCMC$par)
plot(summary(sR),xlab="time")
```
![](outs/sensRangeTimeCrs.png)

The next plot shows what happens when parameters vary by 25%. 
```
parRange=cbind(min=0.75*pars, max=1.25*pars)
crlfun=function(pars) return(m=mean(biExp(pars)$lrat))
CRL=modCRL(fun=crlfun,parRange=parRange,num=500)
cor(CRL)[5,]
plot(CRL,ylab="lrat",cex=0.5,trace=TRUE)
```
![](outs/globalSens.png)
Weak dependence of lrat on Ac0 goes with slow Ac0 convergence.



