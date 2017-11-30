#'iron metabolism model of Parmar et al. 2017
#'
#'This function returns the right hand side of an ordinary differential equation
#'model of iron metabolism in mouse fed a normal iron diet, published by Parmar
#'et al in BMC Systems Biology in 2017.  The intended use of this function is as
#'an argument to \code{ode()} of the \code{deSolve} package.
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
#'@seealso \code{\link{myelo-package}, \link{fokas91}}
#'@references Jignesh H. Parmar, Grey Davis, Hope Shevchuk and Pedro Mendes,
#'  Modeling the dynamics of mouse iron body distribution: hepcidin is necessary
#'  but not sufficient,  \emph{BMC Systems Biology}, \bold{11}, 57 (2017)
#'@keywords IO
#'@export
#'@import deSolve

parmar17<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)),{
    Fe1Tf=TfTot-Tf-Fe2Tf #monos = total - free - duos (all in Moles)
    FeTotalN=(FeDuo+FeLiver+FeSplee+FeRBC+FeRest+2*Fe2Tf+Fe1Tf+NTBI+FeBM) * quantit # Total_Fe (particle)
    aveConc=FeTotalN/(6.02214e+23*(Duodenu+Liver+Plasma+RBC+RestOfB+Spleen)) #Total_Fe (conc)
    FePlasmN=(2*Fe2Tf+Fe1Tf+NTBI)*quantit # FePlasma (particle)
    gTot=FeTotalN*55.845/6.02214e+23 # Total_Fe (g); atomic weight of Fe is 55.845 gm/Mole
    TfTotN=(Tf+Fe1Tf+Fe2Tf)*quantit # Total_Tf (particle)
    TfTotC=TfTotN/(6.02214e+23*Plasma) # Total_Tf (conc) M/L
    TfTotDen=TfTotC*80 #Tf is 80 kDa so 80,000 gm/Mole so TfTotC (M/L)x(80000 gm/M)x(L/1000 mL), so this is gm/ml, not mg/ml  
    FePlas=FePlasmN/(Plasma*6.02214e+23) # FePlasma(conc)
    TfSatur=100*(2*Fe2Tf +Fe1Tf)/(2*(Tf+Fe2Tf+Fe1Tf)) #TfSaturation
    Tf_c=Tf/Plasma
    FeDuo_c=FeDuo/Duodenu
    FeRest_c=FeRest/RestOfB
    FeBM_c=FeBM/BoneMar
    NTBI_c=NTBI/Plasma
    FeLiver_c=FeLiver/Liver
    FeSplee_c=FeSplee/Spleen
    Hepcidi_c=Hepcidi/Plasma
    Fe2Tf_c=Fe2Tf/Plasma
    FeRBC_c=FeRBC/RBC
    Fe1Tf_c=Fe1Tf/Plasma
    FeOutsi_c=FeOutsi/RestOfB
    
    VhepDeg=k1*Hepcidi_c
    Vduo=VmaxDuo*Duodenu*FeDuo_c/((Km+FeDuo_c)*(1+Hepcidi_c/Ki))
    massAc=kInBM*Fe2Tf_c*Plasma
    massAc1=vRBCSpl*FeRBC_c*RBC
    Vspl=VmaxSpl*Spleen*FeSplee_c/((Km+FeSplee_c)*(1+Hepcidi_c/Ki))
    Vliv=VmaxLiv*Liver*FeLiver_c/((Km+FeLiver_c)*(1+Hepcidi_c/Ki))
    MassAc2=kNTBI_F*NTBI_c*Tf_c
    massAc3=kDuoLos*FeDuo_c*Duodenu
    massAc4=kInLive*Fe2Tf_c*Plasma
    massAc5=kInRest*Fe2Tf_c*Plasma
    MassAc6=kRestLo*FeRest_c
    massAc7=kInDuo*Fe2Tf_c*Plasma
    Vrest=VmaxRest*RestOfB*FeRest_c/((Km+FeRest_c)*(1+Hepcidi_c/Ki))
    MassAc8=kFe1Tf_*Fe1Tf_c*NTBI_c
    massA=kInLive*Fe1Tf_c*Plasma
    massA1=kInBM*Fe1Tf_c*Plasma
    massA2=kInRest*Fe1Tf_c*Plasma
    VinDuo=kInDuo*Fe1Tf_c*Plasma
    massA4=kInRBC*FeBM_c*BoneMar
    massA5=kBMSple*FeBM_c*BoneMar
    
    # Equations:
    dTf=massAc-MassAc2*Plasma+massAc4+massAc5+massAc7+massA+massA1+massA2+VinDuo
    dFeDuo=vDiet*Duodenu-Vduo-massAc3+2*massAc7+VinDuo
    dFeRest=2*massAc5-MassAc6*RestOfB-Vrest+massA2
    dFeBM=2*massAc+massA1-massA4-massA5
    dNTBI=Vduo+Vspl+Vliv-MassAc2*Plasma+Vrest-MassAc8*Plasma
    dFeLiver=-Vliv+2*massAc4+massA
    dFeSplee=massAc1-Vspl+massA5
    dHepcidi=vHepSyn*Plasma-VhepDeg*Plasma
    dFe2Tf=-massAc-massAc4-massAc5-massAc7+MassAc8*Plasma
    dFeRBC=-massAc1+massA4
    list(c(dTf, dFeDuo, dFeRest, dFeBM, dNTBI, dFeLiver, dFeSplee,dHepcidi,dFe2Tf,dFeRBC),
         c(Fe1Tf_c=Fe1Tf_c,Fe2Tf_c=Fe2Tf_c,FeRBC_c=FeRBC_c,VinBM=massA1,VinRBC=massA4,VinSpl=massA5,gTot=gTot)
    )
  }) # end with.
}
