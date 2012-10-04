library(myelo)   # load definition of function fokas91
library(deSolve)
# first do CML state since up higher in plots
pars=c(	fb=0.8, # blast fraction remaining active per division (nb = 3 is assumed)
		fp=0.8, # progenitor fraction remaining active per division (np=2 is assumed)
		fm=0.7, # myelocyte fraction remaining active per division  (nm=2 is assumed)
		T= 20, # hours per stage  I'M NOT SURE IF THIS AND THE NEXT STAY AS IS
		q=72e6/16  # flux into the blast pool per hour per kg of marrow
)  # 

X0=rep(0,14)
names(X0)<-c(paste("A",1:7,sep=""),paste("M",1:7,sep=""))
times <- seq(1, 700, by = 1)
out   <- ode(X0, times, fokas91, pars)
tail(out)  # Mismatch with Table 3: Nm/Np = 1.9 < 2.3, Nm = 2.1e9 < 2.6e9,=one problem, c(2.3/1.9,2.6/2.1)

X0=out[700,2:15] # start at steady state, then hit with chemo dropping all cells to 1%
(eventdat <- data.frame(var = names(X0), time = 0, value = 0.01, method = "mult"))
out   <- ode(X0, -20:300, fokas91, pars,events = list(data = eventdat))
tail(out) 

plot(out[ , 1], out[ , "Ntot"], type = "l", xlab = "time", ylab = "Marrow Cells per kg of body mass",
		main = "Model of Fokas et al. 1991", log="y", lwd = 2,col=1)

# now do normal state
pars=c(	fb=0.8, # blast fraction remaining active per division (nb = 3 is assumed)
		fp=0.6, # progenitor fraction remaining active per division (np=2 is assumed)
		fm=0.5, # myelocyte fraction remaining active per division  (nm=2 is assumed)
		T= 20, # hours per stage
		q=72e6/16  # flux into the blast pool per hour per kg of marrow
)  # 

X0=rep(0,14)
names(X0)<-c(paste("A",1:7,sep=""),paste("M",1:7,sep=""))
times <- seq(1, 700, by = 1)
out   <- ode(X0, times, fokas91, pars)
tail(out)  # Mismatch with Table 3: Nm/Np = 1.9 < 2.3, Nm = 2.1e9 < 2.6e9,=one problem, c(2.3/1.9,2.6/2.1)

X0=out[700,2:15] # start at steady state, then hit with chemo dropping all cells to 1%
(eventdat <- data.frame(var = names(X0), time = 0, value = 0.01, method = "mult"))
out   <- ode(X0, -20:300, fokas91, pars,events = list(data = eventdat))
tail(out) 

lines(out[ , 1], out[ , "Ntot"],lwd = 2,col=2)
legend(150,1e8, c("CML","Normal"), col = 1:2, bty="n",lwd=2)
# recoveries from chemo in 7 days as shown here seem reasonable

