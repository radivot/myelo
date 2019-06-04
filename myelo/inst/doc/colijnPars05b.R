# Caroline Colijn and Michael C. Mackey
# A mathematical model of hematopoiesis: II. Cyclical neutropenia
# Journal of Theoretical Biology 237 (2005) 133-146

colijnPars05b=c(
		Qss=1.1, # 1e6  cells/kg
		gamS=0.07, #0.1, # 1/day
		tauS=2.8, #days1,2
		k0=8.0,  #3.0, # days
		the2=0.3, #0.5e6, # 1e6  cells/kg
		s=4,
		
		Nss=7, #  1e8 cells/kg or 3.55   6.9e8 in Z19
		gamN=2.4, # 1/day
		tauNM=3.5, # days
		# An=75200, # none  # 70635 in Z19, 
		An=752, # none  # 70635 in Z19, 
		# etaNbar=0.27, # 1/days
		f0=0.40, # 1/days
		the1=0.36, # 1e8  cells/kg 
		n=1, # n=1
		# the1script=79.7e8, # cells/kg 
		# nu1=1, # 1/days
		
		# kR=0.00469,  # 1/days
		Rss=3.5,  # 1e11 cell/kg
		gamR=0.001, # 1/day
		tauRM=6, # days
		tauSum=120, # days
		tauRet=2.8, # days
		# Ar=56300,  #563000, # none  #lost factor of 10 since 2005a !!!!!!!
		Ar=5.63,  #563000, # none  #lost factor of 10 since 2005a !!!!!!!
		Krbar=1.17,  #1.1, # 1/day
		Kr=0.0382, #  1e11 cells/kg. Deal with units raised to -m ??
    m=6.96,		
		
		Pss=2.94,  #2.14e10, #3.1071e10, # 1e10  cell/kg
		gamP=0.15,   #0.05, # 1/day
		tauPM=7, # days
		tauPS=9.5, # days
		Kpbar=.372, # 1/day
		# Ap=28200,  #282000, #524288, #=2^19 #### Down Factor of 10 again!!!!
		Ap=28.2,  #282000, #524288, #=2^19 #### Down Factor of 10 again!!!!
		etaPbar=0.55, # 1/days
		Kp=11.66, #   units of (1e10 cells/kg)-s4
		r=1.29 
		# the4script=35.9e10, # cells/kg 
		# nu4=1 # 1/days
)
save(colijnPars05b,file="colijnPars05b.rda")

