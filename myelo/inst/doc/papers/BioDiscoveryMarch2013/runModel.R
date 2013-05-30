#library(deSolve) 
#library(devtools)
# install.packages("Rcpp")
#install_github("myelo",subdir="myelo",username="radivot")

library(myelo)   # load definition of function rcj12
################################ the following is from the rcj12 help page
pars=parsBase=c(
		Vx.p=1.7e4, Kil3=1, 
#		Vx.p=4.25e3, Kil3=0.1, # in the Appendix 
		Vp.p=.0203, Kgm=2, Kcn=100, 
		Qf=.0373, Vq.p=.53778, Kp=7e4,
		Vp.n=.026, Kg=1, 
#		kp.m=1.309e-3, # progen to monocytes
		Vn.tn=.781, Kil8=396, Ksb=250, 
		Vn.x=2e-3, 	Kmcl1=2.438, 
		Vtn.x=.038, 
		Vtn.an=.001, Knap2=.1, 
		Kmo=1,
		Van.tn=.09,  Vx.gm=2.007, Kn=5e5, Vgm.x=.05,
		Vx.nap2=.002, Kabl=1, Kdas=5, Vnap2.x=.2, Vx.rac2=1.5, Kan=100, Vrac2.x=.1, 
		Vx.mcl1=1, 
		Kabl2=1, Vmcl1.x=.5, Vx.ros=980, Krac2=250,
		Vgsh.gssg=1680, Vgssg.gsh=150,Kgssg=30,Kros=30, Kgsh=1e5, Vn.m=1.309e-3,Kros2=.026,
		Kseli=.3, Kifna=1, Vbact.x=5.756e-2, nHr=2.5,kkmax=320,Vx.abl=.22,Ktnfa=1,Kimab=1, Vabl.x=.0145,Vil8.x=5e-4,
		Vtnfa.x=5e-6,  
		g=1, il3=1, gsh=1e6,                          
		 volP=2900, volM=1400, volT=65700,
		bu=0,mo=0,das=0,sb=0,cn=0,seli=0,imab=0,ifna=0,antiox=7)      

# Starting values for state variables
X0=X0Base=c(P=6.468063e6, Q=5.656005e5, N=4.556524e6, TN=7.524667e4, AN=76.00674,
		GM=3.969130,NAP2=.01, RAC2=26.40101, Mcl1=2.0, ROS=1.950843, bact=0, abl=0, il8=4, tnfa=0, GSSG=49.297)  

#6468063 565600.5 4556524 75246.67 76.00674 3.969130 0.01000 26.40101 2.000 1.950843    0   0   4    0 49.29748 
#this steady state is the same for both rcj12 and rcj12conIL8 below 

times <- seq(1,80, by = 1)
out   <- ode(X0, times, rcj12, pars)
matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
		main = "Simulation 1: Uninhibited control", log="y", lwd = 2)
legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")
################################## the rest extends the rcj12 help page, running simulations
tail(out)   # steady state checks out with appendix

###### END Simulation 1

## Simulation 2  IL8 increase to 50 nM constant => need to change function

rcj12conIL8 <- function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				vX.P =Vx.p*il3/(il3+Kil3)
				vP.P =P*Vp.p*GM/(GM+Kgm)/(1+cn/Kcn)*(1-Qf)
				vP.N =P*Vp.n*g/(g+Kg)
				vN.TN =N*Vn.tn*il8/(il8+Kil8)/(1+sb/Ksb)
				vTN.X =TN*Vtn.x/(1+Mcl1/Kmcl1)
				vTN.AN =TN*Vtn.an*NAP2/(NAP2+Knap2)
				vAN.TN =AN*Van.tn
				vP.Q =P*Vp.p*GM/(GM+Kgm)/(1+cn/Kcn)*Qf
				vQ.P =Q*Vq.p/(1+P/Kp)
				vX.GM =Vx.gm/(1+N/Kn)/(1+mo/Kmo)
				vGM.X=GM*Vgm.x
				vX.NAP2=Vx.nap2*(1+abl/(Kabl*(1+das/Kdas)))
				vNAP2.X=NAP2*Vnap2.x
				vX.RAC2=Vx.rac2*(1+AN/Kan)
				vRAC2.X=RAC2*Vrac2.x
				vN.X=N*Vn.x/(1+Mcl1/Kmcl1)
				vX.MCL1=Vx.mcl1*(1+abl/Kabl2)/(1+seli/Kseli)/(1+ifna/Kifna)
				vMCL1.X=Mcl1*Vmcl1.x
				vX.ROS=RAC2*Vx.ros/(RAC2+Krac2)
				vGSH.GSSG=ROS*Vgsh.gssg/(ROS+Kros)*(gsh/Kgsh/(1+gsh/Kgsh))
				vN.M=P*Vn.m
				vGSSG.GSH=GSSG*Vgssg.gsh/(GSSG+Kgssg)
				vROS.X=ROS*antiox*Kros2
				vBACT.X=bact*Vbact.x*ROS^nHr/(ROS^nHr+kkmax)
				vX.ABL=tnfa/(1+imab/Kimab)*Vx.abl/(tnfa+Ktnfa)
				vABL.X=abl*Vabl.x
				vIL8.X=Vil8.x*bact*Vbact.x*ROS^nHr/(ROS^nHr+kkmax)
				vTNFA.X=Vtnfa.x*bact*Vbact.x*ROS^nHr/(ROS^nHr+kkmax)
				
				dP =vX.P+vP.P-vP.N-vP.Q+vQ.P-vN.M;
				dQ =vP.Q-vQ.P;
				dN =vP.N*volM/volP-vN.TN-vN.X
				dTN=vN.TN*volP/volT-vTN.X-vTN.AN+vAN.TN
				dAN=vTN.AN-vAN.TN
				dGM =vX.GM-vGM.X
				dNAP2   =vX.NAP2-vNAP2.X
				dRAC2   =vX.RAC2-vRAC2.X
				dMcl1   =vX.MCL1-vMCL1.X
				dROS    =vX.ROS-vGSH.GSSG-vROS.X
				dBACT   =-vBACT.X
				dabl    =vX.ABL-vABL.X
				dil8    =0
				dtnfa   =-vTNFA.X
				dGSSG   =vGSH.GSSG-vGSSG.GSH
				return(list(c(dP,dQ,dN,dTN,dAN,dGM,dNAP2,dRAC2,dMcl1,dROS,dBACT,dabl,dil8,dtnfa,dGSSG),
								c(vP.P=vP.P,vP.N=vP.N,vN.TN=vN.TN,vTN.X=vTN.X,vTN.AN=vTN.AN,vAN.TN=vAN.TN,vP.Q=vP.Q,vQ.P=vQ.P,
										vX.GM=vX.GM,vGM.X=vGM.X,vX.NAP2=vX.NAP2,vNAP2.X=vNAP2.X,vX.RAC2=vX.RAC2,vN.X=vN.X,
										vX.MCL1=vX.MCL1,vMCL1.X=vMCL1.X,tvN.M=vN.M,vBACT.X=vBACT.X,vIL8.X=vIL8.X,
										vTNFA.X=vTNFA.X,v14=vRAC2.X,v18=vX.ROS,
										v19=vGSH.GSSG,v21=vGSSG.GSH,v22=vROS.X,v23=vX.ABL,v24 =vABL.X))
				)
			})
}

# Starting values for state variables
X0=c(P=6.466e6, Q=5.655e5, N=4.558e6, TN=7.564e4, AN=76.84,
		GM=3.97,NAP2=.01004, RAC2=26.59, Mcl1=2.021, ROS=1.967, bact=0, abl=0, il8=4, tnfa=0, GSSG=0)  
times <- seq(-10,100, by = 1)
(IL8add <- data.frame(var = "il8", time = 0, value = 12.5,method = "mult"))
out   <- ode(X0, times, rcj12conIL8, pars,events = list(data = IL8add))
matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
		main = "Simulation 2. IL8 change: 4 nM -> 50 nM", log="y", lwd = 2)
legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")
tail(out)   # steady state checks out with appendix
############################################# End Sim 2


############# Start Sim 3
pars["sb"]=5000  # = IL8 antagonist
out   <- ode(X0, 0:300, rcj12conIL8, pars)
matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
		main = "Simulation 3. SB: 0 uM -> 5 uM", log="y", lwd = 2)
legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")


############# Start Sim 4
X0["il8"]=50
out   <- ode(X0, 0:300, rcj12conIL8, pars)
matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
		main = "Simulation 4. SB: 0 uM -> 5 uM + il8 4 -> 50", log="y", lwd = 2)
legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")

# set back to baseline
X0["il8"]=4
pars["sb"]=0


pars=parsBase
X0=X0Base

############# Start Sim 5
pars["mo"]=10
out   <- ode(X0, seq(0,480,40), rcj12, pars)
matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
		main = "Simulation 5. mo: 0 uM -> 10 nM", log="y", lwd = 2)
legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")
tail(out)

pars["mo"]=0
############# Start Sim 6
pars["g"]=5
out   <- ode(X0, seq(0,480,40), rcj12conIL8, pars)
matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
		main = "Simulation 6. G-CSF 1 uM -> 5 nM", log="y", lwd = 2)
legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")


pars["g"]=1
pars=parsBase
X0=X0Base

############# Start Sim 7
(GMmult <- data.frame(var = "GM", time = 0, value = 12.5,method = "mult"))
out   <- ode(X0, seq(-100,300,10), rcj12conIL8, pars, events = list(data = GMmult))
matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
		main = "Simulation 7. GM-CSF: 4 uM -> 50 nM", log="y", lwd = 2)
legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")

# next three are bateria simulations
pars=parsBase
X0=X0Base
############# Start Sim 8
(bactadd <- data.frame(var = "bact", time = 0, value = 0e6,method = "add"))
#out   <- ode(X0, seq(-50,5e5,1e4), rcj12conIL8, pars, events = list(data = bactadd))
out   <- ode(X0, seq(-50,5e5,1e4), rcj12, pars, events = list(data = bactadd))
#out   <- ode(X0, seq(-50,250,20), rcj12conIL8, pars, events = list(data = bactadd))
matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
		main = "Simulation 8. bacteria: 1e6", log="y", lwd = 2)
legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")
rbind(head(out),tail(out))
############# Start Sim 9
(bactadd <- data.frame(var = "bact", time = 0, value = 1e5,method = "add"))
out   <- ode(X0, seq(-50,250,20), rcj12conIL8, pars, events = list(data = bactadd))
matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
		main = "Simulation 8. bacteria: 1e5", log="y", lwd = 2)
legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")
out
############# Start Sim 10
(bactadd <- data.frame(var = "bact", time = 0, value = 1e7,method = "add"))
out   <- ode(X0, seq(-50,250,20), rcj12conIL8, pars, events = list(data = bactadd))
matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
		main = "Simulation 8. bacteria: 1e7", log="y", lwd = 2)
legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")
out

tail(out)

