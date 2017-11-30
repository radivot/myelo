# Jignesh H. Parmar, Grey Davis, Hope Shevchuk and Pedro Mendes,
#  Modeling the dynamics of mouse iron body distribution: hepcidin is necessary
#  but not sufficient,  BMC Systems Biology, 11:57 (2017)
# these estimates all correspond to column #6 in Table 2
parmarPars=c(FeOutsi=0, Duodenu=3.84615e-05, RBC=0.00079, Spleen=6.73077e-05, Liver=0.0011619, 
             Plasma=0.0013, RestOfB=0.0196948, BoneMar=0.000214286, 
             # organ volumes in Liters (Table 1) mouse= ~25gm = ~.023L
             kNTBI_F=1.08432e+09,  #kNTBI_Fe1Tf
             kInDuo=0.0689984, kInLive=2.97791, kInRBC=1.08448, kInRest=6.16356, 
             Km=0.0159421, Ki=1e-09, kFe1Tf_=1.08432e+09,   #kFe1Tf_Fe2Tf
             VmaxDuo=0.200242, VmaxLiv=0.0261148, VmaxSpl=1.3422, VmaxRest=0.0109451, # all 4 NTBI, in Moles/Day
             kDuoLos=0.0270113,  #kDuoLoss
             vRBCSpl=0.0235286, kRestLo=0.0235332, kInBM=15.7691, kBMSple=0.061903, 
             vDiet=0.00377422,vDietH=0.00415624, vDietL=0,           #Hi vs. low Fe diet
             quantit=6.02214e+23, # Avogadr=6.02214e+23, 
             #HepcidinSynthesis Hi vs. Low Fe diet (adaptation?)
             vHepSyn=1.7393e-08, vHepSynH=2.30942e-08,  vHepSynL=8.54927e-09, 
             k1=0.75616  #reaction 'HepcidinDecay':  kinetic ,eter 'k1'
)
sum(parmarPars[1:8]) #23.2 gm mouse
save(parmarPars,file="parmarPars.rda")

