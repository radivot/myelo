 # myeloL8B.R PK/PD Model of CML 1 March 2014
 #
 library(deSolve)
 
 # P   = progenitors (old CFU-GM)
 # Q   = quiescent progenitors
 # N   = Neutrophils
 #TN   = Tissue Neutrophils
 #AN   = Activated Neutrophils
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
 
 # PK parameters for current drug (imatinib mesylate)
 # ke = elimination rate constant, hr-1
 # ka = rate constant for intestinal absorption, hr-1
 # Fpc = percent oral bioavailability, F%
 # Vd = volume of distribution, L/kg
 # Vg = gut volume, L/kg
 # MW = molecular weight
 # Kdas = IC50, nM
 
 fmyelo <- function(Time, State, Pars) {
           with(as.list(c(State, Pars)), {
           redox=gsh/(gsh+2*GSSG)
           g=gb 
           das=Dp				                     #drug being modelled is TKI
 
 	 vX.P =Vx.p*il3/(il3+Kil3)
	 vP.P =P*Vp.p*GM/(GM+Kgm)*(1-Qf)*(1+MYC/Kmyc)
 	 vP.N =P*Vp.n*g/(g+Kg)
 	 vN.TN =N*Vn.tn*il8/(il8+Kil8)	
	 vTN.X =TN*Vtn.x/(1+Mcl1/Kmcl1)
 	 vTN.AN =TN*Vtn.an*stat5/(stat5+Kstat5)
 	 vAN.TN =AN*Van.tn
	 vP.Q =P*Vp.p*GM/(GM+Kgm)*Qf*(1+MYC/Kmyc)
 	 vQ.P =Q*Vq.p/(1+P/Kp)
 	 vX.GM =Vx.gm*(1+stat5/Kastat5)
 	 vGM.X=GM*Vgm.x*(1+N/Kn)
 	 vX.STAT5=Vx.stat5*(1+BcrAbl*redox/(Kabl*(1+das/Kdas)))
	 vSTAT5.X=stat5*Vstat5.x
 	 vX.STAT3=Vx.stat3*(1+BcrAbl*redox/(Kabl3*(1+das/Kdas)))
 	 vSTAT3.X=stat3*Vstat3.x
 	 vN.X=N*Vn.x/(1+Mcl1/Kmcl1)
 	 vX.MCL1=Vx.mcl1*stat5/(stat5+Kstat5m)/(1+seli/Kseli)
 	 vMCL1.X=Mcl1*Vmcl1.x
 	 vX.ROS=stat3*Vx.ros/(stat3+Kstat3)
 	 vGSH.GSSG=ROS*Vgsh.gssg/(ROS+Kros)*(gsh/Kgsh/(1+gsh/Kgsh))
 	 vP.M=P*Vp.m
 	 vGSSG.GSH=GSSG*Vgssg.gsh/(GSSG+Kgssg)
 	 vROS.X=ROS*antiox*Kros2
  	 vX.ABL=0
 	 vABL.X=0
 	 vGSH.X=peitc*Vgsh.x/(peitc+Kgsh.x)
 	 vDg.Dp=Dg/Vg*ka
 	 vDp.x=Dp*ke
 
 dP =vX.P+vP.P-vP.N-vP.Q+vQ.P-vP.M
 dQ =vP.Q-vQ.P;
 dN =vP.N*volM/volP-vN.TN-vN.X
 dTN=vN.TN*volP/volT-vTN.X-vTN.AN+vAN.TN
 dAN=vTN.AN-vAN.TN
 dGM =vX.GM-vGM.X
 dstat5   =vX.STAT5-vSTAT5.X
 dstat3  =vX.STAT3-vSTAT3.X
 dMcl1   =vX.MCL1-vMCL1.X
 dROS    =vX.ROS-vGSH.GSSG-vROS.X
 dabl    =vX.ABL-vABL.X
 dg=(vX.ABL-vABL.X)/(Ka2*(1+das/Kdas))
 dil3    =(vX.STAT5-vSTAT5.X)/Km0
 dMYC    =(vX.STAT5-vSTAT5.X)/Kmyc
 dgsh=2.0*vGSSG.GSH-vGSH.GSSG-vGSH.X
 dGSSG   =(vGSH.GSSG+vGSH.X)/2.0-vGSSG.GSH
 dDg     =-vDg.Dp
 dDp=vDg.Dp*Fpc/100/Vd-vDp.x

 return(list(c(dP,dQ,dN,dTN,dAN,dGM,dstat5,dstat3,dMcl1,dROS,dabl,dg,dil3,dMYC,dgsh,dGSSG,dDg,dDp),
        c(v0=vX.P,vP.P=vP.P,vP.N=vP.N,vN.TN=vN.TN,vTN.X=vTN.X,vTN.AN=vTN.AN,vAN.TN=vAN.TN,vP.Q=vP.Q,vQ.P=vQ.P,
 vX.GM=vX.GM,vGM.X=vGM.X,vX.STAT5=vX.STAT5,vSTAT5.X=vSTAT5.X,vX.STAT3=vX.STAT3,vN.X=vN.X,
 vX.MCL1=vX.MCL1,vMCL1.X=vMCL1.X,vP.M=vP.M,v14=vSTAT3.X,v18=vX.ROS,
 v19=vGSH.GSSG,v20=vP.M,v21=vGSSG.GSH,v22=vROS.X,v25=vGSH.X))
    )
  })
 }
 
 
 pars=c(Ka2=.5, Kabl=1, Kabl2=1, Kabl3=1, Kan=100, Kastat5=.1, Kdas=85, Kg=1, Kgm=2, Kgsh=30, Kgsh.x=5.0, Kgssg=30, Kil3=100, Kil8=396,
      Km0=1.23e-3, Kmcl1=2.438, Kmyc=1.0, Kn=5e5, Kp=7e4, Kros=30, Kros2=.104, Kseli=.3, Kstat3=250, Kstat5=.1, Kstat5m=1,  
      Vabl.x=.0145, Van.tn=.09, Vbact.x=5.756e-2, Vgm.x=4.943e-3, Vgsh.x=6.0e4, Vgsh.gssg=382.8, Vgssg.gsh=77.0, Vmcl1.x=.5, Vn.tn=.781, Vn.x=2e-3, Vp.m=1.309e-3, 
      Vp.n=.026, Vp.p=.00203, Vq.p=.53778, Vstat3.x=.1, Vstat5.x=.2, Vtn.an=.001, Vtn.x=.038, 
      Vx.abl=.22, Vx.gm=.1372, Vx.mcl1=106.8, Vx.p=5.97e5,Vx.ros=980, Vx.stat3=6.67, Vx.stat5=.002,   
      
      il8=4,Qf=.0373, volP=2900, volM=1400, volT=65700,                    
      MYC=9.01, Vg=.045,                                                                         # these are boundary  conditions
 
      ke=.0173, ka=.61, Fpc=98.0, Vd=4.89, MW=589.7,                                             # PK/PD parameters for imatinib
#     ke=.161, ka=.98, Fpc=80, Vd=2.9, MW=354,                                                   # PK/PD parameters for seliciclib
 
      das=0,seli=0,antiox=7,peitc=0)                                                             # Drug doses, mg/m2
 
      X0=c(P=9.72e6, Q=1.00e6, N=1.50e7, TN=1.21e6, AN=6711, GM=1.79, stat5=.09901, stat3=661.3, Mcl1=19.25, ROS=482.0,  
      BcrAbl=13, gb=51.2, il3=36, MYC=10.2, gsh=7e5, GSSG=1.5e5, 				 # Starting values for variables
      
      Dg=0 *22857.1/589.7, Dp=0)                 						 # Enter dose of drug being modelled in mg/m^2

#  initial drug dose in nanomole/L = mg/m^2 * 1.6 * 1e6 /70/MW

 times <- seq(1,80, by = 1)
 out   <- ode(X0, times, fmyelo, pars)

 head(out)
 tail(out)

 X0=out[80,2:19]
 
 matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
 main = "myeloL8B.R: Uninhibited control", log="y", lwd = 2)
 legend(20,2e5, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
 "activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")
