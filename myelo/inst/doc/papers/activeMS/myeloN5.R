# myeloN5.R Model of the innate immune response             3 June 2017
# 
setwd("C:/R")
library(deSolve)
# state variables are capitalized
# P   = progenitors (old CFU-GM)
# Q   = quiescent progenitors
# N   = Neutrophils
#TN   = Tissue Neutrophils
#AN   = Activated Neutrophils
# D   = dead cell
#GM   = GM-CSF   
#g    = G-CSF
#NAP2 = neutrophil activating peptide-2, driven by STAT5 (not currently modelled)
#STAT3 = mitochondrial transcription factor associated with increased ROS production.
#STAT5 = transcription factor activated by abl
#RAC2 = G protein associated with increased ROS production, driven by STAT3 (not currently modelled).
#Mcl1 = mantle cell lymphoma 1 (antiapoptosis factor)
#ROS  = reactive oxygen species
#Mb = monoblasts
#M   = monocytes
#Mphage = macrophages
#AM = activated macrophages
#My = myeloblasts
#M_CSF = total (free+bound) CSF. M_CSFf = free (unbound) M-CSF
#LPS = lipopolysaccharide, bacterial cell wall component
#TNF = tumour necrosis factor
#CCL2 = cytokine (C-C motif) ligand 2

# Drugs used to manipulate the system
#bu      = Busulfan = DNA crosslinker  (1)  drop state variable P by 10^(2*bu/20)
#mo     = MOR103 = Ab blocker of GM   (2)
#das     = dasatinib                   (3)
#sb       = SB272844 = IL8 inhibitor    (4)
#cn       = CNDAC = sapacitabine becomes this = SSB producer = DSBs in S via HR  (5)
#seli     = seliciclib = cdk9 inhibitor, blocks transcription of Mcl-1 (6)
#imab   = imatinib = anti-TNF monoclonal antibody (7)
#ifna     = interferon alpha (8)
#antiox = ascorbic acid, antoxidant (9)
#leva   = levamisole, macrophage activator
#DAC = 2'-deoxy-5-azacytidine, hypomethylating agent

bacteria=1.0e+3; bact=bacteria; Kbact=1.0e-3; Flps=0.5; Rcount=4e-6; KiM_CSF=8.5e-7; QI=5.655e5; Q=QI; Tet2=1.0; flag99=0; if(bacteria)flag99=1;t99=1e4; 
DAC=0;stat5I=.01004; stat5=stat5I; stat3I=26.59; stat3=stat3I
PI=6.78e6; P=PI; MI=3.0e5; M=MI; NI=4.558e6; N=NI; TNI=7.564e4; TN=TNI; ANI=76.84; AN=ANI; eryth=1.52e8; MKP=1.52e4; MphageI=2.64e4; Mphage=MphageI; AMI=43.7; 
AM=AMI; GMI=3.73; GM=GMI; Mcl1I=2.021; Mcl1=Mcl1I; ROSI=1.967; ROS=ROSI; g=1; il3=1; gshI=1.0e6; gsh=gshI; GSSGI=0; GSSG=GSSGI;
MyI=2.23e6; My=MyI; MbI=6.83e5; Mb=MbI; CCL2I=1.27e3; CCL2=CCL2I; M_CSFI=7.5; M_CSF=M_CSFI; abl=bacteria*5e-6; il8=4+bacteria*Kbact*Flps; monotot=M+Mphage+AM; 
KaHMA=1.0; M_CSFR=monotot*Rcount*(Tet2+(1-Tet2)*DAC/KaHMA)
M_CSFf=M_CSF-M_CSFR; if(M_CSFf<0)M_CSFf=0;
weightP=1.0; weightE=0.1; weightMK=1.5; 
CFUlimit=(PI+MyI+MbI)*weightP+eryth*weightE+MKP*weightMK
Vbact.x = 1.5
ix = 5 # inner loop index
nlines=24; dim=(nlines+1)*16; timex=0; timeW=9999; flag1=1
Cmatrix<-matrix(1:dim, ncol=16)

colnames(Cmatrix)<-c("Time","P","M","N","TN","AN","Mphage","AM","bacteria","Q","GM","stat5","stat3","Mcl1","MCSF","MCSFf")


fmyelo <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    vX.P =Vx.p*(1+il3/Ka0)
    vP.P =P*Vp.p*GM/(GM+Kgm)/(1+cn/Kcn)*(1-Qf)
    vP.My =P*Vp.my*g/(g+Kg)
    vN.TN =N*Vn.tn*il8/(il8+Kil8)/(1+sb/Ksb)
    vTN.X =TN*Vtn.x/(1+Mcl1/Kmcl1)
    vTN.AN =TN*Vtn.an*stat5/(stat5+Kstat5)
    vAN.TN =AN*Van.tn
    vP.Q =P*Vp.p*GM/(GM+Kgm)/(1+cn/Kcn)*Qf
    vQ.P =Q*Vq.p/(1+P/Kp)
    vX.GM =Vx.gm/(1+N/Kn)/(1+mo/Kmo)
    vGM.X=GM*Vgm.x
    vX.STAT5=Vx.stat5 *(1+abl/(Kabl*(1+das/Kdas)))
    vSTAT5.X=stat5*Vstat5.x
    vX.STAT3=Vx.stat3*(1+AN/Kan)
    vSTAT3.X=stat3*Vstat3.x
    vN.X=N*Vn.x/(1+Mcl1/Kmcl1)
    vX.MCL1=Vx.mcl1*(1+abl/(Kabl2*(1+das/Kdas)))/(1+seli/Kseli)/(1+ifna/Kifna)
    vMCL1.X=Mcl1*Vmcl1.x
    vX.ROS=stat3*Vx.ros/(stat3+Kstat3)
    vGSH.GSSG=ROS*Vgsh.gssg/(ROS+Kros)*(gsh/Kgsh/(1+gsh/Kgsh))
    vP.Mb=P*Vp.mb/(P+Kp.mb)*M_CSFf/KaM_CSF
    vGSSG.GSH=GSSG*Vgssg.gsh/(GSSG+Kgssg)
    vROS.X=ROS*antiox*Kros2
    vBACT.BACT=bact*Kbac
    vBACT.X=bact*Vbact.x*ROS^nHr/(ROS^nHr+kkmax)
    vX.ABL=tnfa/(1+imab/Kimab)*Vx.abl/(tnfa+Ktnfa)
    vABL.X=abl*Vabl.x
    vIL8.X=Vil8.x*bact*Vbact.x*ROS^nHr/(ROS^nHr+kkmax)
    vTNFA.X=Vtnfa.x*bact*Vbact.x*ROS^nHr/(ROS^nHr+kkmax)
    vMy.N=My*k32
    vMy.My=My*Vmy.my*GM/(GM+Kgm)/(1+cn/Kcn)
    vMb.Mb=Mb*Vmb.mb*GM/(GM+Kgm)/(1+cn/Kcn)
    vMb.M=Mb*k30
    vM.Mphage=M*k26
    vMphage.X=Mphage*k27
    vMphage.AM = Vmphage.am *Mphage/Kmacro /(1+Mphage/Kmacro)* (LPS/Klps+ leva/Kleva)/ (1+LPS/Klps+leva/Kleva)*((CCL2/Kaccl2)^nHccl2)/(1+(CCL2/Kaccl2)^nHccl2)
    vX.AM = Vx.am
    vAM.X = AM*k29
    vX.LPS=vBACT.BACT*Kbact
    vLPS.X=vBACT.X*Kbact
    
    dP       = vX.P+vP.P-vP.My-vP.Q+vQ.P-vP.Mb;
    dM        = vMb.M*volM/volP-vM.Mphage;
    dN       = vMy.N*volM/volP-vN.TN-vN.X
    dTN     = vN.TN*volP/volT-vTN.X-vTN.AN+vAN.TN
    dAN     = vTN.AN-vAN.TN
    dMphage = vM.Mphage-vMphage.X
    dAM   = vMphage.AM+vX.AM-vAM.X
    dQ       = vP.Q-vQ.P;
    dGM    = vX.GM-vGM.X
    dstat5   = vX.STAT5-vSTAT5.X
    dstat3   = vX.STAT3-vSTAT3.X
    dMcl1  = vX.MCL1-vMCL1.X
    dROS    = vX.ROS-vGSH.GSSG-vROS.X
    dbact     = vBACT.BACT-vBACT.X; if(bact<.1)dbact=0
    dabl       = vX.ABL-vABL.X
    dil8        = dbact * Kbact * Flps
    dtnfa      = dbact * 5e-6
    dg	= (vX.ABL-vABL.X)/(Ka2*(1+das/Kdas))
    dil3 = (vX.ABL-vABL.X)/(Km0*(1+das/Kdas));if(il3<.1)dil3=0
    dgsh	= vGSSG.GSH-vGSH.GSSG
    dGSSG = vGSH.GSSG-vGSSG.GSH
    dMy      = vP.My+vMy.My-vMy.N
    dMb      = vP.Mb+vMb.Mb-vMb.M
    dLPS     = vX.LPS-vLPS.X;if(LPS<0){LPS=0;dLPS=0}
    dCCL2    = dAM*Kccl2
    dM_CSF = dAM*KiM_CSF
    
    return(list(c(dP,dM,dN,dTN,dAN,dMphage,dAM,dbact,dQ,dGM,dstat5,dstat3,dMcl1,dROS,dabl,dil8,dtnfa,dg,dil3, dgsh,dGSSG,dMy,dMb,dLPS,dCCL2,dM_CSF),
                c(vP.P=vP.P,vP.My=vP.My,vN.TN=vN.TN,vTN.X=vTN.X,vTN.AN=vTN.AN,vAN.TN=vAN.TN,vP.Q=vP.Q,vQ.P=vQ.P,vX.GM=vX.GM,vGM.X=vGM.X,vX.STAT5=vX.STAT5,vSTAT5.X=vSTAT5.X,vX.STAT3=vX.STAT3,vN.X=vN.X,vX.MCL1=vX.MCL1,vMCL1.X=vMCL1.X,vB.B=vBACT.BACT,vBACT.X=vBACT.X,vIL8.X=vIL8.X,vTNFA.X=vTNFA.X,v14=vSTAT3.X,v18=vX.ROS,v19=vGSH.GSSG,v21=vGSSG.GSH,v22=vROS.X,v23=vX.ABL,vX.P =vX.P, v26=vM.Mphage, v27=vMphage.X, v28=vMphage.AM, v29=vAM.X,v30=vMb.M, v31=vMb.Mb, v32=vMy.N, v33=vMy.My, M_CSF=M_CSF))
    )
  }	)
}

pars=c(Vx.p=4.168e4, Ka0=1, Vp.p=7.4e-3, Kgm=2, Kcn=100, Vp.my=.026, Kg=1, Kmo=1,  Vn.tn=.781, Kil8=396, Ksb=250, Vtn.x=.038, Kmcl1=2.438, Vtn.an=.001, Kstat5=.1, Van.tn=.09, Vq.p=.53778, Kp=7e4, Vx.gm=2.007, Kn=5e5, Vgm.x=.05, Vx.stat5=.002, Kabl=1, Kdas=5, Vstat5.x=.2, Vx.stat3=1.5, Kan=100, Vstat3.x=.1, Vn.x=2e-3, Vx.mcl1=1, Kabl2=1, Vmcl1.x=.5, Vx.ros=980, Kstat3=250, Vgsh.gssg=382.8, Vgssg.gsh=155.8,Kgssg=30,Kros=30, Kgsh=9.5e5, Vp.mb=2.0e5, Kp.mb=5.0e6, KaM_CSF=29.0, Kros2=.104, Kseli=.3, Kifna=1,  nHr=2.5,kkmax=1000,Vx.abl=.22,Ktnfa=1,Kimab=1, Vabl.x=.0145,Vil8.x=5e-4, Vtnfa.x=5e-6, Ka2=.5, Km0=.1, Klps=1.0e3, Vmy.my=2.406e-4, k32=.04, Vmb.mb=2.406e-4, k30=.04, k26=.04, k27=0.5, Kbac=1.134, Vmphage.am=5.0e7, Kmacro=3.5e5,k29=.08,Kccl2=1.0,Kaccl2=1e4, nHccl2=3.0, Kleva=200, Vx.am=3.5,
       Qf=.0373, volP=2900, volM=1400, volT=65700, M_CSFf=M_CSF-M_CSFR,     
       bu=0,mo=0,das=0,sb=0,cn=0,seli=0,imab=0,ifna=0,antiox=7,leva=0) # Drug concentrations

X0=c(P=PI, M=MI, N=NI, TN=TNI, AN=ANI, Mphage=MphageI, AM=AMI, bact=bacteria, Q=QI, GM=GMI,stat5=stat5I, stat3=stat3I, Mcl1=Mcl1I, ROS=ROSI, abl=bacteria*5e-6, il8=4+bacteria*Kbact*Flps, tnfa=bacteria*5e-6, g=1, il3=1, gsh=gshI, GSSG=GSSGI, My=MyI, Mb=MbI, LPS=bacteria*Kbact, CCL2=CCL2I, M_CSF=M_CSFI)                  # Starting values for variables

Cmatrix[1,1]=0; Cmatrix[1,2]=P;Cmatrix[1,3]=M; Cmatrix[1,4]=N; Cmatrix[1,5]=TN; Cmatrix[1,6]=AN; Cmatrix[1,7]=Mphage; Cmatrix[1,8]=AM; Cmatrix[1,9]=bact; Cmatrix[1,10]=Q; Cmatrix[1,11]=GM; Cmatrix[1,12]=stat5; Cmatrix[1,13]=stat3; Cmatrix[1,14]=Mcl1; Cmatrix[1,15]=M_CSF; Cmatrix[1,16]=M_CSFf 

times <- seq(1,ix, by = 1)
out   <- ode(X0, times, fmyelo, pars)

P=out[ix,2]

for(i in seq(from=2, to=nlines, by=1)) {
  M=out[ix,3];N=out[ix,4];TN=out[ix,5];AN=out[ix,6];Mphage=out[ix,7];AM=out[ix,8];bact=out[ix,9];Q=out[ix,10]; GM=out[ix,11];stat5=out[ix,12];stat3=out[ix,13];
  Mcl1=out[ix,14];ROS=out[ix,15];abl=out[ix,16];il8=out[ix,17]; tnfa=out[ix,18];g=out[ix,19];il3=out[ix,20];gsh=out[ix,21];GSSG=out[ix,22];My=out[ix,23];Mb=out[ix,24]; 
  LPS=out[ix,25];CCL2=out[ix,26];M_CSF=out[ix,27];monotot=M+Mphage+AM;
  M_CSFR=monotot*Rcount*(Tet2+(1-Tet2)*DAC/KaHMA);M_CSFf=M_CSF-M_CSFR;if(M_CSFf<0)M_CSFf=0;
  if(flag99==1 && bact<1){flag99=0; t99=out[ix,1]}
  
  X0=c(P,M,N,TN,AN,Mphage,AM,bact,Q,GM,stat5,stat3,Mcl1,ROS,abl,il8,tnfa,g,il3,gsh,GSSG,My,Mb,LPS,CCL2, M_CSF)
  
  Cmatrix[i,1] =  out[ix,1];  Cmatrix[i,2] = out[ix,2] ; Cmatrix[i,3 ]=out[ix,3];    Cmatrix[i,4]=out[ix,4];
  Cmatrix[i,5] =  out[ix,5];  Cmatrix[i,6] = out[ix,6];  Cmatrix[i,7 ]=out[ix,7];    Cmatrix[i,8]=out[ix,8]; 
  Cmatrix[i,9] =  out[ix,9];  Cmatrix[i,10]=out[ix,10]; Cmatrix[i,11]=out[ix,11]; Cmatrix[i,12]=out[ix,12];
  Cmatrix[i,13]=out[ix,13]; Cmatrix[i,14]=out[ix,14];Cmatrix[i,15] =M_CSF; Cmatrix[i,16]=M_CSFf 
  
  times<-seq(ix*(i-1)+1,ix*i, by=1)
  out <-ode(X0, times, fmyelo, pars) 
  
  #Competition function
  progenitors=(P+My+Mb)*weightP+eryth*weightE+MKP*weightMK;
  if(progenitors>CFUlimit){ ratioP=P*weightP/progenitors; ratioMy=My*weightP/progenitors;
  ratioMb=Mb*weightP/progenitors; ratioE=eryth*weightE/progenitors;
  ratioMK=MKP*weightMK/progenitors;
  P=ratioP*CFUlimit/weightP;My=ratioMy*CFUlimit/weightP;Mb=ratioMb*CFUlimit/weightP;
  eryth=ratioE*CFUlimit/weightE; MKP=ratioMK*CFUlimit/weightMK; }
}

Cmatrix[i+1,1] =  out[ix,1];  Cmatrix[i+1,2] = out[ix,2] ; Cmatrix[i+1,3 ]=out[ix,3];    Cmatrix[i+1,4]=out[ix,4];
Cmatrix[i+1,5] =  out[ix,5];  Cmatrix[i+1,6] = out[ix,6];  Cmatrix[i+1,7 ]=out[ix,7];    Cmatrix[i+1,8]=out[ix,8]; 
Cmatrix[i+1,9] =  out[ix,9];  Cmatrix[i+1,10]=out[ix,10]; Cmatrix[i+1,11]=out[ix,11]; Cmatrix[i+1,12]=out[ix,12];
Cmatrix[i+1,13]=out[ix,13]; Cmatrix[i+1,14]=out[ix,14];Cmatrix[i+1,15]=M_CSF; Cmatrix[i+1,16]=M_CSFf

write.table(Cmatrix, file="tableM.txt", append=FALSE, quote=FALSE, sep=" ", row.names=TRUE, col.names=TRUE);
write(t99, file="tableM.txt", append=TRUE)

Cmatrix
#if(flag99==0 && Cmatrix[i+1,9]<1 && t99<1e4){print("bacteria eliminated at ");print(t99)}

matplot(Cmatrix[,1],Cmatrix[,2:9],type="l",xlab="time",ylab="count",main="myeloN5.R: bacteria=1.0e+3",log="y",lwd=2)
legend(8,1e4,c("progenitors","monocytes","blood neutro","tissue neutro","activated neutro","macrophages","activated macros","bacteria"),col=1:8,lty=1:8,bty="n")






