library(myelo)   # load definition of function fokas91
# model of Michor et al (2005) Nature, 435, 1267-70 
# coded by Tom Radivoyevitch on 8/6/2005
d0=0.003;d1=0.008;d2=0.05;d3=1
ax=0.8
bx=5
cx=100
ry=0.008
ay=2*ax
byy=2*bx
cy=cx
rz=0.023
az=ay
bz=byy
cz=cy

X0o=1e7
X0o=2e6 # needed to get ratio at t=0 about right
X1o=ax*X0o/d1  # at steady state RHS = 0 = ax*X0o - d1*X1o 
X2o=bx*X1o/d2
X3o=cx*X2o/d3
Y0o=2.5e5
Y1o=ay*Y0o/d1
Y2o=byy*Y1o/d2
Y3o=cy*Y2o/d3

"baoverbRatio"
Y3o/(Y3o+X3o)  

(y0<-c(X0=X0o,X1=X1o,X2=X2o,X3=X3o,Y0=Y0o,Y1=Y1o,Y2=Y2o,Y3=Y3o,Z0=0,Z1=0,Z2=0,Z3=0))
m=.5  # This slope tuning paramater of the controller lambda is immaterial 
#  since X's are decoupled from other state variables. 

# first show without mutations and stopping therapy after 400 days (i.e. b in figure 4)
y0["Z0"]=0
y0["Z1"]=0
y0["Z2"]=0
y0["Z3"]=0
u=0   
library(myelo)
out1=lsoda(y=y0,times=seq(0,365,1),michor05, parms=c(trt=1), rtol=1e-4, atol= rep(1e-4,12))
ny0=out1[nrow(out1),2:13]
out2=lsoda(y=ny0,times=seq(365,500,1),michor05, parms=c(trt=0), rtol=1e-4, atol= rep(1e-4,12))
outs=data.frame(rbind(out1,out2))
attach(outs)

par(mfcol=c(5,2))
plot(time,Y0,type="l",log="y",ylab="SC",xlab="Time (days)",main="Treatment halted at 1 year")
plot(time,Y1,type="l",log="y",ylab="PC",xlab="Time (days)")
plot(time,Y2,type="l",log="y",ylab="DC",xlab="Time (days)")
plot(time,Y3,type="l",log="y",ylab="TC",xlab="Time (days)")
plot(time,ratio*100,type="l",log="y",ylab="BCR-ABL/BCR",xlab="Time (days)")
detach(outs)


# now keep therapy on but include resistance growth

y0["Z0"]=10
y0["Z1"]=az*10/d1
y0["Z2"]=bz*y0["Z1"]/d2
y0["Z3"]=cz*y0["Z2"]/d3

u=4e-7
u=0  # starting off with 10 resistant cells 5 lines up, so no mutation was used
out1=lsoda(y=y0,times=seq(0,500,1),michor05, parms=c(trt=1), rtol=1e-4, atol= rep(1e-4,12))
outs=as.data.frame(out1)
attach(outs)
plot(time,Y0,type="l",log="y",ylab="SC",xlab="Time (days)",main="Emerging Resistance",ylim=range(c(Y0,Z0),finite=T,na.rm=T))
lines(time,Z0)
plot(time,Y1,type="l",log="y",ylab="PC",xlab="Time (days)",ylim=range(c(Y1,Z1),finite=T,na.rm=T) )
lines(time,Z1)
plot(time,Y2,type="l",log="y",ylab="DC",xlab="Time (days)",ylim=range(c(Y2,Z2),finite=T,na.rm=T))
lines(time,Z2)
plot(time,Y3,type="l",log="y",ylab="TC",xlab="Time (days)",ylim=range(c(Y3,Z3),finite=T,na.rm=T) )
lines(time,Z3)
plot(time,ratio*100,type="l",log="y",ylab="BCR-ABL/BCR",xlab="Time (days)")
par(mfrow=c(1,1))

detach(outs)



