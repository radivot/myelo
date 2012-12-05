library(myelo)

scholzPars

Z<-function(x,min,max,nor,b) {
	if ((nor<max)&(nor>min)) 
		res= max- (max-min)*exp(-log((max-min)/(max-nor))*x^b) else res=nor
	res
} 


(X=10^c(-5:10))


