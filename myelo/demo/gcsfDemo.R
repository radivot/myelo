library(myelo)

scholzPars

#Z<-function(x,min,max,nor,b){
#	if((nor<max)&(nor>min))
#		res=max-(max-min)*exp(-log((max-min)/(max-nor))*x^b)elseres=nor
#	res
#}
#this equivalent form shows better that
#x=1=>Z=nor   #x=inf=>Z=max   #x=0=>Z=min
Z<-function(x,min,max,nor,b)res=max-(max-min)*((max-nor)/(max-min))^(x^b)


(x=10^c(-12:5))
attach(as.list(scholzPars))
c(AminPGBF,AmaxPGBF,AnorPGBF,AbPGBF)
c(AminPGBP,AmaxPGBP,AnorPGBP,AbPGBP)
yF=Z(x,AminPGBF,AmaxPGBF,AnorPGBF,AbPGBF)
plot(x,yF,log="xy",ylab="Z",xlab="G-CSF",ylim=c(.01,1e3))
yP=Z(x,AminPGBP,AmaxPGBP,AnorPGBP,AbPGBP)
lines(x,yP)
#Figure5 claims an X-axis scale factor of 10^10, but I get 10^5 from this plot!
#If x=0 on a logscale,why rescale x? 
data.frame(x,yF,yP)



parsF=c(ksc=kFsc,km=kFm,vmax=vFmax,ku=kFu,kcp=kFcp,kpc=kFpc,VD=VFD,vGRAmax=vGRAFmax,kGRAm=kGRAFm)
parsP=c(ksc=kPsc,km=kPm,vmax=vPmax,ku=kPu,kcp=kPcp,kpc=kPpc,VD=VPD,vGRAmax=vGRAPmax,kGRAm=kGRAPm)





detach(as.list(scholzPars))





