# Changjing Zhuge, Jinzhi Lei, and Michael C. Mackey
# Neutrophil dynamics in response to chemotherapy and G-CSF
# Journal of Theoretical Biology 293 (2012) 111–120

zhugePars=c(
		Qss=1.1e6, #cells/kg
		gamS=0.07, # 1/day
		gamMinS=0.03,
		gamMaxS=0.20, 
		tauS=2.8, #days1,2
		k0=8.0, # days
		the2=0.3e6, # cells/kg
		s2=4,
		Nss=6.3e8, # cells/kg
		gamN=2.4, # 1/day
		tauNP=5, # days
		tauNM=6, # days
		tauNMgcsf=2, # days
		tauN=11, # days
		etaNP=2.5420, # 1/days
		etaMinNP=2.0420, # 1/days
		etaMaxNP=3.0552, # 1/days
		gam0=0.27, # 1/day
		gamMin0=0.12, 
		f0=0.40, # 1/days
		the1=0.36e8, # cells/kg 
		s1=1, 
		kdel=0.01, # 1/days
		T=21, # days
      T1=4 # day
)
save(zhugePars,file="zhugePars.rda")

