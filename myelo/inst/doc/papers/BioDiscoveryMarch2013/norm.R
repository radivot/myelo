fnorm <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    # let these two ride as auxilliary variable relative to ABL
    g      =ABL/(Ka2*(1+das/Kdas))
    il3    =ABL/(Km0*(1+das/Kdas))
    il8    =BACT * 5e-4 +4
    tnfa   =BACT * 5e-6 +0.01   
    
    vX.P =Vx.p*il3/(il3+Kil3)
    vP.P =P*Vp.p*GM/(GM+Kgm)/(1+cn/Kcn)*(1-Qf)
    vP.N =P*Vp.n*g/(g+Kg)
    vN.TN =N*Vn.tn*il8/(il8+Kil8)/(1+sb/Ksb)
    vTN.X =TN*Vtn.x/(1+MCL1/Kmcl1)
    vTN.AN =TN*Vtn.an*NAP2/(NAP2+Knap2)
    vAN.TN =AN*Van.tn
    vP.Q =P*Vp.p*GM/(GM+Kgm)/(1+cn/Kcn)*Qf
    vQ.P =Q*Vq.p/(1+P/Kp)
    vX.GM =Vx.gm/(1+N/Kn)/(1+mo/Kmo)
    vGM.X=GM*Vgm.x
    vX.NAP2=Vx.nap2*(1+ABL/(Kabl*(1+das/Kdas)))
    vNAP2.X=NAP2*Vnap2.x
    vX.RAC2=Vx.rac2*(1+AN/Kan)
    vRAC2.X=RAC2*Vrac2.x
    vN.X=N*Vn.x/(1+MCL1/Kmcl1)
    vX.MCL1=Vx.mcl1*(1+ABL/(Kabl2*(1+das/Kdas)))/(1+seli/Kseli)/(1+ifna/Kifna)
    vMCL1.X=MCL1*Vmcl1.x
    vX.ROS=RAC2*Vx.ros/(RAC2+Krac2)
    vGSH.GSSG=ROS*Vgsh.gssg/(ROS+Kros)*(GSH/Kgsh/(1+GSH/Kgsh))
    vP.M=P*Vp.m
    vGSSG.GSH=GSSG*Vgssg.gsh/(GSSG+Kgssg)
    vROS.X=ROS*antiox*Kros2
    vBACT.X=BACT*Vbact.x*ROS^nHr/(ROS^nHr+kkmax)
    vX.ABL=tnfa/(1+imab/Kimab)*Vx.abl/(tnfa+Ktnfa)
    vABL.X=ABL*Vabl.x

    
    dP =vX.P+vP.P-vP.N-vP.Q+vQ.P-vP.M;
    dQ =vP.Q-vQ.P;
    dN =vP.N*volM/volP-vN.TN-vN.X
    dTN=vN.TN*volP/volT-vTN.X-vTN.AN+vAN.TN
    dAN=vTN.AN-vAN.TN
    dGM =vX.GM-vGM.X
    dNAP2   =vX.NAP2-vNAP2.X
    dRAC2   =vX.RAC2-vRAC2.X
    dMCL1   =vX.MCL1-vMCL1.X
    dROS    =vX.ROS-vGSH.GSSG-vROS.X
    dBACT   =-vBACT.X
    dABL    =vX.ABL-vABL.X
    dGSH     =vGSSG.GSH-vGSH.GSSG
    dGSSG   =vGSH.GSSG-vGSSG.GSH
    return(list(c(dP,dQ,dN,dTN,dAN,dGM,dNAP2,dRAC2,dMCL1,dROS,dBACT,dABL,dGSH,dGSSG),
                c(g=g,il3=il3,
                  vP.P=vP.P,vP.N=vP.N,vN.TN=vN.TN,vTN.X=vTN.X,vTN.AN=vTN.AN,vAN.TN=vAN.TN,vP.Q=vP.Q,vQ.P=vQ.P,
                  vX.GM=vX.GM,vGM.X=vGM.X,vX.NAP2=vX.NAP2,vNAP2.X=vNAP2.X,vX.RAC2=vX.RAC2,vN.X=vN.X,
                  vX.MCL1=vX.MCL1,vMCL1.X=vMCL1.X,vP.M=vP.M,vBACT.X=vBACT.X,v14=vRAC2.X,v18=vX.ROS,
                  v19=vGSH.GSSG,v21=vGSSG.GSH,v22=vROS.X,v23=vX.ABL,v24 =vABL.X))
    )
  })
}



pars=c(
#        Vx.p=1.7e4, Kil3=1, 
        Vx.p=4.2e4, Kil3=0.1, # values in Appendix. Didn't help much
       Vp.p=.0203, Kgm=2, Kcn=100, Vp.n=.026, Kg=1, Kmo=1,
       Vn.tn=.781, Kil8=396, Ksb=250, Vtn.x=.038, Kmcl1=2.438, Vtn.an=.001, Knap2=.1, 
       Van.tn=.09, Vq.p=.53778, Kp=7e4, Vx.gm=2.007, Kn=5e5, Vgm.x=.05,
       Vx.nap2=.002, Kabl=1, Kdas=5, Vnap2.x=.2, Vx.rac2=1.5, Kan=100, Vrac2.x=.1, 
       Vn.x=2e-3, Vx.mcl1=1, Kabl2=1, Vmcl1.x=.5, Vx.ros=980, Krac2=250,
       Vgsh.gssg=382.8, Vgssg.gsh=155.8,Kgssg=30,Kros=30, Kgsh=9.5e5, Vp.m=1.309e-3,Kros2=.104,
       Kseli=.3, Kifna=1, Vbact.x=5.756e-2, nHr=2.5,kkmax=320,Vx.abl=.22,Ktnfa=1,Kimab=1, Vabl.x=.0145,Vil8.x=5e-4,
       Vtnfa.x=5e-6, Ka2=.5, Km0=.1,                                 
       Qf=.0373, volP=2900, volM=1400, volT=65700,        # these are boundary  conditions
       bu=0,mo=0,das=0,sb=0,cn=0,seli=0,imab=0,ifna=0,antiox=7)      # Drug concentrations

# X0=c(P=6.466e6, Q=5.655e5, N=4.558e6, TN=7.564e4, AN=76.84,
#      GM=3.97,NAP2=.01004, RAC2=26.59, MCL1=2.021, ROS=1.967, BACT=0, ABL=0, GSH=1e6, GSSG=0)                                                                            # Starting values for variables
# When bact>0 initialise il8 to bact*5e-4 +4
X0=structure(c(51856311.7168859, 19729860.5277474, 17012655.6091217, 
            299967.506077166, 343.818967027626, 1.14602836074423, 0.011502219187436, 
            66.5728450541439, 2.3004438374872, 84.2664337973667, 0, 0.150221918743598, 
            999607.375711163, 392.624288836607), .Names = c("P", "Q", "N", 
            "TN", "AN", "GM", "NAP2", "RAC2", "MCL1", "ROS", "BACT", "ABL", 
            "GSH", "GSSG"))

