 # myeloLV.R   PD Model of Bcr-Abl+ cell lines                1 March 2014
 # This version models in vitro drug exposure.
 #
 library(deSolve)
 
 # P   = progenitors (old CFU-GM)
 # Q   = quiescent progenitors
 # N   = Neutrophils
 #TN   = Tissue Neutrophils
 #AN   = Activated Neutrophils
 # D   = dead cell
 #GM   = GM-CSF   
 #stat5 = nuclear transcription factor activated by Bcr-Abl
 #stat3 = mitochondrial transcription factor associated with increased ROS production
 #Mcl1 = mantle cell lymphoma 1 (antiapoptosis factor)
 #ROS  = reactive oxygen species
 
 # Drugs used to manipulate the system
 #das= dasatinib                   
 #seli = seliciclib = cdk9 inhibitor, blocks transcription of Mcl-1
 #antiox = ascorbic acid, antoxidant 
 #peitc, inducer of oxidative stress  
 
  
 fmyelo <- function(Time, State, Pars) {
           with(as.list(c(State, Pars)), {
           redox=gsh/(gsh+2*GSSG) 				                     
 
 	 vP.P =P*Vp.p*GM/(GM+Kgm)*(1-Qf)*(1+MYC/Kmyc)
 	 vP.Q =P*Vp.p*GM/(GM+Kgm)*Qf*(1+MYC/Kmyc)
 	 vQ.P =Q*Vq.p/(1+P/Kp)
 	 vX.STAT5=Vx.stat5*(1+BcrAbl*redox/(Kabl*(1+das/Kdas)))
	 vSTAT5.X=stat5*Vstat5.x
 	 vX.STAT3=Vx.stat3*(1+BcrAbl*redox/(Kabl3*(1+das/Kdas)))
 	 vSTAT3.X=stat3*Vstat3.x
 	 vP.X=P*Vp.x/(1+Mcl1/Kmcl1)
 	 vX.MCL1=Vx.mcl1*stat5/(stat5+Kstat5m)/(1+seli/Kseli)
 	 vMCL1.X=Mcl1*Vmcl1.x
 	 vX.ROS=stat3*Vx.ros/(stat3+Kstat3)
 	 vGSH.GSSG=ROS*Vgsh.gssg/(ROS+Kros)*(gsh/Kgsh/(1+gsh/Kgsh))
 	 
 	 vGSSG.GSH=GSSG*Vgssg.gsh/(GSSG+Kgssg)
 	 vROS.X=ROS*antiox*Kros2
 	 vGSH.X=peitc*Vgsh.x/(peitc+Kgsh.x)

 
 dP     = vP.P-vP.Q+vQ.P-vP.X
 dQ     =vP.Q-vQ.P;
 dstat5 =vX.STAT5-vSTAT5.X
 dstat3 =vX.STAT3-vSTAT3.X
 dMcl1  =vX.MCL1-vMCL1.X
 dROS   =vX.ROS-vGSH.GSSG-vROS.X
 dgsh   =2.0*vGSSG.GSH-vGSH.GSSG-vGSH.X
 dGSSG  =(vGSH.GSSG+vGSH.X)/2.0-vGSSG.GSH

 return(list(c(dP,dQ,dstat5,dstat3,dMcl1,dROS,dgsh,dGSSG),
        c(vP.P=vP.P,vP.Q=vP.Q,vQ.P=vQ.P,vX.STAT5=vX.STAT5,vSTAT5.X=vSTAT5.X,vX.STAT3=vX.STAT3,vP.X=vP.X,
 vX.MCL1=vX.MCL1,vMCL1.X=vMCL1.X,vSTAT3.X=vSTAT3.X,vX.ROS=vX.ROS,v19=vGSH.GSSG,v21=vGSSG.GSH,v22=vROS.X,v25=vGSH.X)) )
  })
 }
 
 
 pars=c(Kabl=1, Kabl3=1, Kdas=85, Kgm=2, Kgsh=30, Kgsh.x=5.0, Kgssg=30, Kmcl1=2.438, Kmyc=1.0, Kp=7e4, Kros=30, Kros2=.104, Kseli=.3, Kstat3=250, Kstat5m=1,  
      Vgsh.x=6.0e4, Vgsh.gssg=382.8, Vgssg.gsh=77.0, Vmcl1.x=.5, Vp.x=2e-3, Vp.p=.00203, Vq.p=.53778, Vstat3.x=.1, Vstat5.x=.2, 
      Vx.mcl1=106.8, Vx.ros=980, Vx.stat3=6.67, Vx.stat5=.002,   
      
      Qf=.0373, MYC=9.01, GM=1.5, BcrAbl=13, MYC=10.2,                                                  # these are boundary  conditions
 
      das=0,seli=0,antiox=7,peitc=0)                                                             # Drug doses, nM except antiox which is mM
 
      X0=c(P=9.72e6, Q=1.00e6, stat5=.09901, stat3=661.3, Mcl1=19.25, ROS=482.0, gsh=7e5, GSSG=1.5e5)	      # Starting values for variables                 						 # Enter dose of drug being modelled in mg/m^2

 times <- seq(1,80, by = 1)
 out   <- ode(X0, times, fmyelo, pars)

 head(out)
 tail(out)

 X0=out[80,2:19]
 
 matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
 main = "myeloLV.R: Uninhibited control", log="y", lwd = 2)
 legend(40,1e5, c("Proliferating Cells", "Quiescent Cells","STAT5","STAT3","Mcl-1"), col = 1:5, lty = 1:5,bty="n")