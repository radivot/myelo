# myeloN3.R Model of Myeloid Cell Maturation and trafficking 14 Oct 2012
# state variables are capitalized
# P   = progenitors (old CFU-GM)
# Q   = quiescent progenitors
# N   = Neutrophils
#TN   = Tissue Neutrophils
#AN   = Activated Neutrophils
# D   = dead cell
#GM   = GM-CSF   
#NAP2 = neutrophil activating peptide-2 
#RAC2 = G protein associated with increased ROS production
#Mcl1 = mantle cell lymphoma 1 (antiapoptosis factor)
#ROS  = reactive oxygen species

# Drugs used to manipulate the system
#bu = Busulfan = DNA crosslinker  (1)  drop state variable P by 10^(2*bu/20)
#mo = MOR103 = Ab blocker of GM   (2)
#das= dasatinib                   (3)
#sb = SB272844 = IL8 inhibitor    (4)
#cn = CNDAC = sapacitabine becomes this = SSB producer = DSBs in S via HR  (5)
#seli = seliciclib = cdk9 inhibitor, blocks transcription of Mcl-1 (6)
#imab = imatinib = anti-TNF monoclonal antibody (7)
#ifna = interferon alpha (8)
#antiox = ascorbic acid, antoxidant (9)

rcj12 <- function(Time, State, Pars) {
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
				dil8    =-vIL8.X
				dtnfa   =-vTNFA.X
				dGSSG   =vGSH.GSSG-vGSSG.GSH
				return(list(c(dP,dQ,dN,dTN,dAN,dGM,dNAP2,dRAC2,dMcl1,dROS,dBACT,dabl,dil8,dtnfa,dGSSG),
								c(vP.P=vP.P,vP.N=vP.N,vN.TN=vN.TN,vTN.X=vTN.X,vTN.AN=vTN.AN,vAN.TN=vAN.TN,vP.Q=vP.Q,vQ.P=vQ.P,
										vX.GM=vX.GM,vGM.X=vGM.X,vX.NAP2=vX.NAP2,vNAP2.X=vNAP2.X,vX.RAC2=vX.RAC2,vN.X=vN.X,
										vX.MCL1=vX.MCL1,vMCL1.X=vMCL1.X,tvN.M=vN.M,vBACT.X=vBACT.X,vIL8.X=vIL8.X,vTNFA.X=vTNFA.X,v14=vRAC2.X,v18=vX.ROS,
										v19=vGSH.GSSG,v21=vGSSG.GSH,v22=vROS.X,v23=vX.ABL,v24 =vABL.X))
				)
			})
}

