library(myelo)

Z<-function(x,min,max,nor,b){
 if((nor<max)&(nor>min))
		res=max-(max-min)*exp(-log((max-min)/(max-nor))*x^b) else res=nor
	res
}
#this equivalent form shows better that
#x=1=>Z=nor   #x=inf=>Z=max   #x=0=>Z=min
Z2<-function(x,min,max,nor,b) res = max - (max-min)*((max-nor)/(max-min))^(x^b)

(x=10^c(-12:5))
length(scholzPars)
attach(as.list(scholzPars))
y=Z(x,AminPGBF,AmaxPGBF,AnorPGBF,AbPGBF)
y2=Z2(x,AminPGBF,AmaxPGBF,AnorPGBF,AbPGBF)
plot(x,y,log="xy",ylab="Z",xlab="G-CSF",ylim=c(.01,1e3))
lines(x,y2,col="green")
detach(as.list(scholzPars))

library(myelo)
X0=rep(1,24)
G6nor=GRAnor=1
PexoIV=0
names(X0)<-c("g","g1","g2","g3","g4","S","CG","PGB",
             "G4a","G4b","G4c","G4d","G4e",
             "G5a","G5b","G5c","G5d","G5e",
             "G6a","G6b","G6c","G6d","G6e","GRA")
times <- seq(1, 30, by = 1)
out   <- ode(X0, times, scholz12, scholzPars)
tail(out)


