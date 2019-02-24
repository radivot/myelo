#'iron metabolism model of Parmar and Mendes 2019
#'
#'This function returns the right hand side of an ordinary differential equation
#'model of iron metabolism in mice published by Parmar and Mendes
#'in Plos Comp Biol in 2019.  The intended use of this function is as
#'an argument to \code{ode()} of the \code{deSolve} package. To emulate hemachromatosis set ksHepci to 0.
#'Our goal here is to reproduce some of the Table 1 values at 365 days.
#'
#'
#'@param Time The parameters of this model do not depend on time: this argument
#'  exists here as a dummy only because deSolve expects it to exist.
#'@param State Vector of current states. The elements of this vector must be
#'  named because within this function's definition is
#'  \code{with(as.list(State,Pars), {XXX} )} and code in XXX presupposes that it
#'  can call named elements as variables with those names.
#'@param Pars Vector of parameter values with \emph{named} (see above) elements.
#'@return A list of length 2 (as expected by deSolve) where the first list
#'  element is the vector of derivatives, i.e. the right hand side of the ODEs,
#'  and the second element of the list is a vector of auxiliary variables that
#'  one may wish to track over time.
#'@author Tom Radivoyevitch
#'@seealso \code{\link{myelo-package}, \link{parmar17}}
#'@references Jignesh H. Parmar and Pedro Mendes,
#'  A computational model to understand mouse iron physiology and disease,  \
#'  emph{PLoS Comp Biol}, \bold{15}, 1 (2019) 
#'@keywords IO
#'@name parmar19
#'@export
#'@import deSolve

parmar19<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)),{
    
    # Fe1Tf=TfTot-Tf-Fe2Tf #monos = total - free - duos (all in Moles)
    # FeTotalN=(FeDuo+FeLiver+FeSplee+FeRBC+FeRest+2*Fe2Tf+Fe1Tf+NTBI+FeBM) * quantit # Total_Fe (particle)
    # aveConc=FeTotalN/(6.02214e+23*(Duodenu+Liver+Plasma+RBC+RestOfB+Spleen)) #Total_Fe (conc)
    # FePlasmN=(2*Fe2Tf+Fe1Tf+NTBI)*quantit # FePlasma (particle)
    # gTot=FeTotalN*55.845/6.02214e+23 # Total_Fe (g); atomic weight of Fe is 55.845 gm/Mole
    # TfTotN=(Tf+Fe1Tf+Fe2Tf)*quantit # Total_Tf (particle)
    # TfTotC=TfTotN/(6.02214e+23*Plasma) # Total_Tf (conc) M/L
    # TfTotDen=TfTotC*80 #Tf is 80 kDa so 80,000 gm/Mole so TfTotC (M/L)x(80000 gm/M)x(L/1000 mL), so this is gm/ml, not mg/ml  
    # FePlas=FePlasmN/(Plasma*6.02214e+23) # FePlasma(conc)
    # TfSatur=100*(2*Fe2Tf +Fe1Tf)/(2*(Tf+Fe2Tf+Fe1Tf)) #TfSaturation
    # Tf_c=Tf/Plasma
    # FeDuo_c=FeDuo/Duodenu
    # FeRest_c=FeRest/RestOfB
    # FeBM_c=FeBM/BoneMar
    # NTBI_c=NTBI/Plasma
    # FeLiver_c=FeLiver/Liver
    # FeSplee_c=FeSplee/Spleen
    # Hepcidi_c=Hepcidi/Plasma
    # Fe2Tf_c=Fe2Tf/Plasma
    # FeRBC_c=FeRBC/RBC
    # Fe1Tf_c=Fe1Tf/Plasma
    # FeOutsi_c=FeOutsi/RestOfB
    
    
    # Assignment Model Entities:
    #model entity 'TfSaturation':assignment
    TfSatur=100*(2*Fe2Tf * 6.02214e+23 +Fe1Tf * 6.02214e+23 )/(2*(Tf * 6.02214e+23 +Fe2Tf * 6.02214e+23 +Fe1Tf * 6.02214e+23 ))
    #model entity 'Total_Fe':assignment
    Total_F=FeDuo * 6.02214e+23 +FeLiver * 6.02214e+23 +FeSplee * 6.02214e+23 +FeRBC * 6.02214e+23 +FeRest * 6.02214e+23 +2*Fe2Tf * 6.02214e+23 +Fe1Tf * 6.02214e+23 +NTBI * 6.02214e+23 +FeBM * 6.02214e+23 
    #model entity 'Total_Fe_concentration':assignment
    Total=Total_F/(6.02214e+23*(Duodenu+Liver+Plasma+RBC+RestOfB+Spleen))
    #model entity 'FePlasma':assignment
    FePlasm=2*Fe2Tf * 6.02214e+23 +Fe1Tf * 6.02214e+23 +NTBI * 6.02214e+23 
    #model entity 'Total_Fe_(g)':assignment
    Total_1=Total_F*55.845/6.02214e+23
    #model entity 'Total_Tf':assignment
    Total_T=Tf * 6.02214e+23 +Fe1Tf * 6.02214e+23 +Fe2Tf * 6.02214e+23 
    #model entity 'Total_Tf_conc':assignment
    Total_2=Total_T/(6.02214e+23*Plasma)
    #model entity 'Total_Tf_mg/ml':assignment
    Total_3=Total_2*76000
    #model entity 'FePlasma(conc)':assignment
    FePlas=FePlasm/(Plasma*6.02214e+23)
    
    
    NTBI_c=NTBI/Plasma
    Tf_c=Tf/Plasma
    FeDuo_c=FeDuo/Duodenu
    FeBM_c=FeBM/BoneMar
    FeRest_c=FeRest/RestOfB
    FeLiver_c=FeLiver/Liver
    Fe1Tf_c=Fe1Tf/Plasma
    FeSplee_c=FeSplee/Spleen
    EPO_c=EPO/Plasma
    Hepcidi_c=Hepcidi/Plasma
    FeRBC_c=FeRBC/RBC
    Fe2Tf_c=Fe2Tf/Plasma
    FeOutsi_c=FeOutsi/RestOfB
    
    #model entity 'FeTftotal':assignment
    FeTftot=2*Fe2Tf_c+Fe1Tf_c
    #model entity 'vDiet':assignment
    vDiet=1*0.00346965
    #model entity 'rel([EPO])':assignment
    relEPO=EPO_c/2.4267e-12
    #model entity 'rel([FeBM])':assignment
    relFeBM=FeBM_c/0.00208051
    #model entity 'rel([FeDuo])':assignment
    relFeDu=FeDuo_c/0.0100222
    #model entity 'rel([FeLiver])':assignment
    relFeLi=FeLiver_c/0.0013
    #model entity 'rel([FeRBC])':assignment
    relFeRB=FeRBC_c/0.023
    #model entity 'rel([FeRest])':assignment
    relFeRe=FeRest_c/0.00040176
    #model entity 'rel([FeSpleen])':assignment
    relFeSp=FeSplee_c/0.00226838
    #model entity 'rel([Hepcidin])':assignment
    relHepc=Hepcidi_c/2.30092e-08
    #model entity 'rel([NTBI])':assignment
    relNTBI=NTBI_c/4e-08
    #model entity 'rel([FeTotal])':assignment
    relFeTo=Total_F/1.72456e+19
    #model entity 'FeTftotal':assignment
    FeTftot=2*Fe2Tf_c+Fe1Tf_c
    

    # VhepDeg=k1*Hepcidi_c
    # Vduo=VmaxDuo*Duodenu*FeDuo_c/((Km+FeDuo_c)*(1+Hepcidi_c/Ki))
    # massAc=kInBM*Fe2Tf_c*Plasma
    # massAc1=vRBCSpl*FeRBC_c*RBC
    # Vspl=VmaxSpl*Spleen*FeSplee_c/((Km+FeSplee_c)*(1+Hepcidi_c/Ki))
    # Vliv=VmaxLiv*Liver*FeLiver_c/((Km+FeLiver_c)*(1+Hepcidi_c/Ki))
    # MassAc2=kNTBI_F*NTBI_c*Tf_c
    # massAc3=kDuoLos*FeDuo_c*Duodenu
    # massAc4=kInLive*Fe2Tf_c*Plasma
    # massAc5=kInRest*Fe2Tf_c*Plasma
    # MassAc6=kRestLo*FeRest_c
    # massAc7=kInDuo*Fe2Tf_c*Plasma
    # Vrest=VmaxRest*RestOfB*FeRest_c/((Km+FeRest_c)*(1+Hepcidi_c/Ki))
    # MassAc8=kFe1Tf_*Fe1Tf_c*NTBI_c
    # massA=kInLive*Fe1Tf_c*Plasma
    # massA1=kInBM*Fe1Tf_c*Plasma
    # massA2=kInRest*Fe1Tf_c*Plasma
    # VinDuo=kInDuo*Fe1Tf_c*Plasma
    # massA4=kInRBC*FeBM_c*BoneMar
    # massA5=kBMSple*FeBM_c*BoneMar
    
    
    Constan=vDiet
    LinearA=ksHepci*FeTftot*KEPOHep^hEPOHep/(KEPOHep^hEPOHep+EPO_c^hEPOHep)
    MassAct=kdHepci*Hepcidi_c
    MMWComp=vDuoNTB*Duodenu*FeDuo_c/((KmFeFPN+FeDuo_c)*(1+Hepcidi_c/KiHepci))
    MassAc=kInBM*EPO_c*Fe2Tf_c*Plasma
    MassAc1=kRBCSpl*FeRBC_c*RBC
    MMWCom=vSpleen*Spleen*FeSplee_c/((KmFeFPN+FeSplee_c)*(1+Hepcidi_c/KiHepci))
    MMWCom1=vLiverN*Liver*FeLiver_c/((KmFeFPN+FeLiver_c)*(1+Hepcidi_c/KiHepci))
    MassAc2=kNTBI_F*NTBI_c*Tf_c
    MassAc3=kDuoLos*FeDuo_c*Duodenu
    MassAc4=kInLive*Fe2Tf_c*Plasma
    MassAc5=kInRest*Fe2Tf_c*Plasma
    MassAc6=kInDuo*Fe2Tf_c*Plasma
    MMWCom2=vRestNT*RestOfB*FeRest_c/((KmFeFPN+FeRest_c)*(1+Hepcidi_c/KiHepci))
    MassAc7=kFe1Tf_*Fe1Tf_c*NTBI_c
    MassAc8=kInLive*Fe1Tf_c*Plasma
    MassA=kInBM*EPO_c*Fe1Tf_c*Plasma
    MassA1=kInRest*Fe1Tf_c*Plasma
    MassA2=kInDuo*Fe1Tf_c*Plasma
    MassA3=kInRBC*EPO_c*FeBM_c*BoneMar
    MassA4=kBMSple*FeBM_c*BoneMar
    UptakeW=vNTBILi*Plasma*(NTBI_c/(KmNTBI+NTBI_c))*(NTBI_c/(KaNTBI+NTBI_c))
    BasalFl=vEPO*KiEPORB^hEPO/(KiEPORB^hEPO+FeRBC_c^hEPO)
    MassA5=kdEPO*EPO_c
    MassA6=kRestLo*FeRest_c
    MassA7=kdTf*Tf_c
    Consta=vTf
    MassA8=kdTf*Fe2Tf_c
    MassA9=kdTf*Fe1Tf_c 
    
    # Equations:
    dNTBI=MMWComp+MMWCom+MMWCom1-MassAc2*Plasma+MMWCom2-MassAc7*Plasma-UptakeW+2*MassA8*Plasma+MassA9*Plasma
    dTf=MassAc-MassAc2*Plasma+MassAc4+MassAc5+MassAc6+MassAc8+MassA+MassA1+MassA2-MassA7*Plasma+Consta*Plasma
    dFeDuo=Constan*Duodenu-MMWComp-MassAc3+2*MassAc6+MassA2
    dFeBM=2*MassAc+MassA-MassA3-MassA4
    dFeRest=2*MassAc5-MMWCom2+MassA1-MassA6*RestOfB
    dFeLiver=-MMWCom1+2*MassAc4+MassAc8+UptakeW
    dFe1Tf=MassAc2*Plasma-MassAc7*Plasma-MassAc8-MassA-MassA1-MassA2-MassA9*Plasma
    dFeSplee=MassAc1-MMWCom+MassA4
    dEPO=BasalFl*Plasma-MassA5*Plasma
    dHepcidi=LinearA*Plasma-MassAct*Plasma
    dFeRBC=-MassAc1+MassA3
    dFe2Tf=-MassAc-MassAc4-MassAc5-MassAc6+MassAc7*Plasma-MassA8*Plasma
    
    list(c(dNTBI, dTf, dFeDuo, dFeBM, dFeRest, dFeLiver,dFe1Tf, dFeSplee,dEPO,dHepcidi,dFeRBC,dFe2Tf),
         c(NTBI_c=NTBI_c,Tf_c=Tf_c,FeDuo_c=FeDuo_c,FeBM_c=FeBM_c,FeRBC_c=FeRBC_c,
           FeLiver_c=FeLiver_c,FeSplee_c=FeSplee_c,FeRest_c=FeRest_c,Hepcidi_c=Hepcidi_c)
         
    # list(c(dTf, dFeDuo, dFeRest, dFeBM, dNTBI, dFeLiver, dFeSplee,dHepcidi,dFe2Tf,dFeRBC),
    #      c(Fe1Tf_c=Fe1Tf_c,Fe2Tf_c=Fe2Tf_c,VinBM=massA1,VinRBC=massA4,VinSpl=massA5,
    #        FePlas=FePlas,FeDuo_c=FeDuo_c,FeBM_c=FeBM_c,FeRBC_c=FeRBC_c,
    #        FeLiver_c=FeLiver_c,FeSplee_c=FeSplee_c,FeRest_c=FeRest_c,Hepcidi_c=Hepcidi_c,
    #        gTot=gTot)
    )
  }) # end with.
}
