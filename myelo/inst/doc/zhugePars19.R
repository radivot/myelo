# Changjing Zhuge, Michael C. Mackey, Jinzhi Lei
# Origins of oscillation patterns in cyclical thrombocytopenia
# Journal of Theoretical Biology 462 (2019) 432–445

zhugePars19=c(
		Sss=1.1, # 1e6 cells/kg
		gamS=0.1, # 1/day
		tauS=2.8, #days1,2
		k0=3.0, # days
		the2=0.5, # 1e6 cells/kg
		s2=4,
		Nss=6.9, # 1e8 cells/kg
		gamN=2.4, # 1/day
		tauNM=3.5, # days
		Anmax=70635, # none
		etaNbar=0.27, # 1/days
		f0=0.40, # 1/days
		the1=0.36, # 1e8 cells/kg 
		s1=1, 
		the1script=79.7, # 1e8 cells/kg 
		nu1=1, # 1/days
		
		kR=0.00469,  # 1/days
		
		Pss=3.1071, #1e10 cell/kg
		gamP=0.05, # 1/day
		tauPM=7, # days
		tauPS=9.5, # days
		Kpbar=.372, # 1/day
		Apmax=524288, #=2^19 # none
		etaPbar=0.55, # 1/days
		Kp=11.66, # 1e10 cells/kg 
		s4=1.29, 
		the4script=35.9, # 1e10 cells/kg 
		nu4=1 # 1/days
)
save(zhugePars19,file="zhugePars19.rda")

