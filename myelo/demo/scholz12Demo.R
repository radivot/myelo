########## !!!!!!!!!!! SCHOLZ ET AL BELOW NEEDS FIXING SO IGNORE FOR NOW ########
# The following code pertains to the G-CSF model of M Scholz TBMM 2012
library(myelo)
# first we show that the Z function form they use is more complicated than need be
Z<-function(x,Amin,Amax,Anor,b){
	if((Anor<Amax)&(Anor>Amin))
		res=Amax-(Amax-Amin)*exp(-log((Amax-Amin)/(Amax-Anor))*x^b) else res=Anor
	res }  # since this 
#equivalent form shows better that x=1=>Z=Anor,  x=inf=>Z=Amax, and  x=0=>Z=Amin
Z2<-function(x,Amin,Amax,Anor,b) res = Amax - (Amax-Amin)*((Amax-Anor)/(Amax-Amin))^(x^b)
(x=10^c(-12:5)) #length(scholzPars) # 128 parameters
with(as.list(scholzPars),{ #prove equal by plotting both: y as points, y2 as lines
			y=Z(x,AminPGBF,AmaxPGBF,AnorPGBF,AbPGBF)
			y2=Z2(x,AminPGBF,AmaxPGBF,AnorPGBF,AbPGBF)
			plot(x,y,log="xy",ylab="Z",xlab="G-CSF",ylim=c(.01,1e3))
			lines(x,y2,col="green")})
# now run a simulation to show that RHS code works at least syntactically
X0=rep(1,24)
names(X0)<-c("g","g1","g2","g3","g4","S","CG","PGB","G4a","G4b","G4c","G4d","G4e",
		"G5a","G5b","G5c","G5d","G5e","G6a","G6b","G6c","G6d","G6e","GRA")
# WARNING: the following is needed because the model is not completely specified! 
APGBout=APGBin=ACGout=ACGin=GSS=G6nor=GRAnor=1
out   <- ode(X0, c(0,10^(1:4)), scholz12, scholzPars)
tail(out) # it hit steady state
(X0=out[dim(out)[1],2:25])
(eventdat <- data.frame(var = "g", time = 0, value = 2, method = "mult"))
out   <- ode(X0, -20:600, scholz12, scholzPars,events = list(data = eventdat))
plot(out)
