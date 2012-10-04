fokas91<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				dA1 = fb*q -  A1/T
				dA2 = 2*fb*A1/T - A2/T
				dA3 = 2*fb*A2/T - A3/T
				dA4 = 2*fp*A3/T - A4/T
				dA5 = 2*fp*A4/T - A5/T
				dA6 = 2*fm*A5/T - A6/T
				dA7 = 2*fm*A6/T - A7/T
				dM1 = (1-fb)*q -  M1/T 
				dM2 = 2*(1-fb)*A1/T - M2/T  + M1/T
				dM3 = 2*(1-fb)*A2/T - M3/T  + M2/T
				dM4 = 2*(1-fp)*A3/T - M4/T  + M3/T
				dM5 = 2*(1-fp)*A4/T - M5/T  + M4/T
				dM6 = 2*(1-fm)*A5/T - M6/T  + M5/T
				dM7 = 2*(1-fm)*A6/T - M7/T  + M6/T
				Nb=A1+A2+A3+M1+M2+M3
				Np=A4+A5+M4+M5
				Nm=A6+A7+M6+M7
				Ntot=Nb+Np+Nm
				Q=(2*A7+M7)/T
				return(list(c(dA1,dA2,dA3,dA4,dA5,dA6,dA7,dM1,dM2,dM3,dM4,dM5,dM6,dM7),
							c(Nb=Nb,Np=Np,Nm=Nm,povb=Np/Nb,movb=Nm/Np,Qovm=Q/Nm,Qovq=Q/q,Ntot=Ntot,Q=Q)))
			})
}
