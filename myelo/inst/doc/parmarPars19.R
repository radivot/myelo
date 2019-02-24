# Jignesh H. Parmar and Pedro Mendes, Plos Comp Biol in 2019

# Fixed Model Entities:
#metabolite 'FeOutside': fixed
FeOutsi=0
#compartment 'Duodenum':fixed
Duodenu=3.84615e-05
#compartment 'RBC':fixed
RBC=0.00079
#compartment 'Spleen':fixed
Spleen=6.73077e-05
#compartment 'Liver':fixed
Liver=0.0011619
#compartment 'Plasma':fixed
Plasma=0.0013
#compartment 'RestOfBody':fixed
RestOfB=0.0196948
#compartment 'BoneMarrow':fixed
BoneMar=0.000214286
#global quantity 'kNTBI_Fe1Tf':fixed
kNTBI_F=1.004e+09
#global quantity 'kInDuo':fixed
kInDuo=0.0405971
#global quantity 'kInLiver':fixed
kInLive=2.11666
#global quantity 'kInRBC':fixed
kInRBC=5.03844e+11
#global quantity 'kInRest':fixed
kInRest=4.78121
#global quantity 'KmFeFPN':fixed
KmFeFPN=0.112511
#global quantity 'KiHepcidinFPN':fixed
KiHepci=6.3e-09
#global quantity 'kFe1Tf_Fe2Tf':fixed
kFe1Tf_=1.004e+09
#global quantity 'vDuoNTBI':fixed
vDuoNTB=0.200907
#global quantity 'vLiverNTBI':fixed
vLiverN=0.0444795
#global quantity 'vSpleenNTBI':fixed
vSpleen=2.06738
#global quantity 'vRestNTBI':fixed
vRestNT=0.0101453
#global quantity 'kDuoLoss':fixed
kDuoLos=6.80738e-05
#global quantity 'kRBCSpleen':fixed
kRBCSpl=0.03
#global quantity 'kInBM':fixed
kInBM=4.06878e+12
#global quantity 'kBMSpleen':fixed
kBMSple=0.103218
#global quantity 'vDiet_basal':fixed
vDiet_b=0.00346965
#global quantity 'KaNTBI':fixed
KaNTBI=0.000255016
#global quantity 'KmNTBI':fixed
KmNTBI=0.000679291
#global quantity 'vNTBILiver':fixed
vNTBILi=14.1511
#global quantity 'fDiet':fixed
fDiet=1
#global quantity 'ksHepcidin':fixed
ksHepci=0.000398766
#global quantity 'vEPO':fixed
vEPO=2.62675e-09
#global quantity 'kdEPO':fixed
kdEPO=4.8
#global quantity 'hEPO':fixed
hEPO=6.5
#global quantity 'KiEPORBC':fixed
KiEPORB=0.01
#global quantity 'KEPOHepcidin':fixed
KEPOHep=5e-12
#global quantity 'hEPOHepcidin':fixed
hEPOHep=4
#global quantity 'hNTBI':fixed
hNTBI=1
#global quantity 'kRestLoss':fixed
kRestLo=0.016862
#global quantity 'vTf':fixed
vTf=1.548e-05
#global quantity 'kdTf':fixed
kdTf=0.4
#global quantity 'kdHepcidin':fixed
kdHepci=0.75616

(pars=dput(ls()))
names(pars)=pars
dput(pars)

parmarPars19=c(BoneMar = BoneMar, Duodenu = Duodenu, fDiet = fDiet, 
  FeOutsi = FeOutsi, hEPO = hEPO, hEPOHep = hEPOHep, hNTBI = hNTBI, 
  KaNTBI = KaNTBI, kBMSple = kBMSple, kdEPO = kdEPO, kdHepci = kdHepci, 
  kdTf = kdTf, kDuoLos = kDuoLos, KEPOHep = KEPOHep, kFe1Tf_ = kFe1Tf_, 
  KiEPORB = KiEPORB, KiHepci = KiHepci, kInBM = kInBM, kInDuo = kInDuo, 
  kInLive = kInLive, kInRBC = kInRBC, kInRest = kInRest, 
  KmFeFPN = KmFeFPN, KmNTBI = KmNTBI, kNTBI_F = kNTBI_F, 
  kRBCSpl = kRBCSpl, kRestLo = kRestLo, ksHepci = ksHepci, 
  Liver = Liver, Plasma = Plasma, RBC = RBC, RestOfB = RestOfB, 
  Spleen = Spleen, vDiet_b = vDiet_b, vDuoNTB = vDuoNTB, 
  vEPO = vEPO, vLiverN = vLiverN, vNTBILi = vNTBILi, vRestNT = vRestNT, 
  vSpleen = vSpleen, vTf = vTf)

# parmarPars=c(FeOutsi=0, Duodenu=3.84615e-05, RBC=0.00079, Spleen=6.73077e-05, Liver=0.0011619, 
#              Plasma=0.0013, RestOfB=0.0196948, BoneMar=0.000214286, 
#              # organ volumes in Liters (Table 1) mouse= ~25gm = ~.023L
#              kNTBI_F=1.08432e+09,  #kNTBI_Fe1Tf
#              kInDuo=0.0689984, kInLive=2.97791, kInRBC=1.08448, kInRest=6.16356, 
#              Km=0.0159421, Ki=1e-09, kFe1Tf_=1.08432e+09,   #kFe1Tf_Fe2Tf
#              VmaxDuo=0.200242, VmaxLiv=0.0261148, VmaxSpl=1.3422, VmaxRest=0.0109451, # all 4 NTBI, in Moles/Day
#              kDuoLos=0.0270113,  #kDuoLoss
#              vRBCSpl=0.0235286, kRestLo=0.0235332, kInBM=15.7691, kBMSple=0.061903, 
#              vDiet=0.00377422,vDietH=0.00415624, vDietL=0,           #Hi vs. low Fe diet
#              quantit=6.02214e+23, # Avogadr=6.02214e+23, 
#              #HepcidinSynthesis Hi vs. Low Fe diet (adaptation?)
#              vHepSyn=1.7393e-08, vHepSynH=2.30942e-08,  vHepSynL=8.54927e-09, 
#              k1=0.75616  #reaction 'HepcidinDecay':  kinetic ,eter 'k1'
# )
# sum(parmarPars[1:8]) #23.2 gm mouse
save(parmarPars19,file="parmarPars19.rda")

