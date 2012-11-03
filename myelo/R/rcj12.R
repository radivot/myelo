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



#'Granulopoiesis model of Robert C. Jackson (2012)
#'
#'This function returns the right hand side of the ODEs of a model that Dr.
#'Robert C. Jackson is actively developing. The intent is not to use this
#'function directly, but rather to pass it as an argument to \code{ode()} of
#'the \code{deSolve} package.
#'
#'This model is a work in progress. Its current form is depicted below as a
#'graph.  
#'\if{html}{\figure{rcj.svg}} \if{latex}{\figure{rcj.png}{options: width=5in}}

#'@param Time The parameters of this model do not depend on time: this argument
#'exists here as a dummy only because deSolve expects it to exist.
#'@param State Vector of current states. The elements of this vector must be
#'named because within this function's definition is
#'\code{with(as.list(State,Pars), {XXX} )} and code in XXX presupposes that it
#'can call named elements as variables with those names.
#'@param Pars Vector of parameter values with \emph{named} (see above)
#'elements.
#'@return A list of length 2 (as expected by deSolve) where the first list
#'element is a vector of the derivatives (i.e. the right hand side of the
#'ODEs), and the second element of the list is a vector of auxiliary variables
#'that one may wish to track over time.
#'@author Robert C. Jackson (\email{bob@@aol.com})
#'@seealso \code{\link{myelo-package}, \link{fokas91}}
#'@keywords IO
#' @export 
#'@examples
#'
#'\dontrun{
#'library(myelo)
#'pars=c(Vx.p=1.7e4, Kil3=1, Vp.p=.0203, Kgm=2, Kcn=100, Vp.n=.026, Kg=1, Kmo=1,
#'		Vn.tn=.781, Kil8=396, Ksb=250, Vtn.x=.038, Kmcl1=2.438, Vtn.an=.001, Knap2=.1, 
#'		Van.tn=.09, Vq.p=.53778, Kp=7e4, Vx.gm=2.007, Kn=5e5, Vgm.x=.05,
#'		Vx.nap2=.002, Kabl=1, Kdas=5, Vnap2.x=.2, Vx.rac2=1.5, Kan=100, Vrac2.x=.1, 
#'		Vn.x=2e-3, Vx.mcl1=1, Kabl2=1, Vmcl1.x=.5, Vx.ros=980, Krac2=250,
#'		Vgsh.gssg=1680, Vgssg.gsh=150,Kgssg=30,Kros=30, Kgsh=1e5, Vn.m=1.309e-3,Kros2=.026,
#'		Kseli=.3, Kifna=1, Vbact.x=5.756e-2, nHr=2.5,kkmax=320,Vx.abl=.22,Ktnfa=1,Kimab=1, Vabl.x=.0145,Vil8.x=5e-4,
#'		Vtnfa.x=5e-6,  
#'		g=1, il3=1, gsh=1e6,                          
#'		Qf=.0373, volP=2900, volM=1400, volT=65700,
#'		bu=0,mo=0,das=0,sb=0,cn=0,seli=0,imab=0,ifna=0,antiox=7)      
#'
#'# Starting values for state variables
#'X0=c(P=6.466e6, Q=5.655e5, N=4.558e6, TN=7.564e4, AN=76.84,
#'		GM=3.97,NAP2=.01004, RAC2=26.59, Mcl1=2.021, ROS=1.967, bact=0, abl=0, il8=4, tnfa=0, GSSG=0)  
#'times <- seq(1,80, by = 1)
#'out   <- ode(X0, times, rcj12, pars)
#'matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
#'		main = "Simulation 1: Uninhibited control", log="y", lwd = 2)
#'legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
#'				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")
#'
#'# DNA cross-linking
#'X0=out[80,2:16]
#'bu=20 # set to MTD
#'(eventdat <- data.frame(var = "P", time = 0, value = 1/10^(2*bu/20), method = "mult"))
#'out   <- ode(X0, -20:600, rcj12, pars,events = list(data = eventdat))
#'
#'matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
#'		main = "Simulation 2: busulfan treatment", log="y", lwd = 2)
#'legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
#'				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")
#'} 
#'
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

