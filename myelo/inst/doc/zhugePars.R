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
		tauN=11, # days
		etaNP=2.5420, # 1/days
		etaMinNP=2.0420, # 1/days
		etaMaxNP=3.0552, # 1/days
		gam0=0.27, # 1/day
		gamMin0=0.12, 
		f0=0.40, # 1/days
		the1=0.36e8, # cells/kg 
		s1=1, 
		kdel=0.01 # 1/days
)
save(zhugePars,file="zhugePars.rda")


brooksPars=c(
		Qss=1.12e6, #cells/kg
		gamS=0.1043, # 1/day
		gamMinS=0.03,
		gamMaxS=0.40, 
		tauS=2.83, #days1,2
		k0=8.0, # days
		the2=0.0826e6, # cells/kg
		s2=2,
		Nss=5.59e8, # cells/kg
		gamN=2.4, # 1/day
bn=0.05, # mg=ml 9
		tauNP=5.9, # days
		tauNM=3.8, # days
		tauN=9.7, # days
bv=0.001, # mg=ml 9
Vmax=3.8,
		etaNP=2.1995, # 1/days
		etaMinNP=0.4, # 1/days
		etaMaxNP=2.5444, # 1/days
AN=1549.58e2,
		gam0=0.27, # 1/day
		gamMin0=0.12, 
		f0=0.154605, # 1/days
		the1=0.0154848e8, # cells/kg 
		s1=.5, 
		kdel=0.0134, # 1/days
Xss=0.1, # ug/kg
Gss=0, # ug/ml 
VB=76, # ml/kg
Gprod=0, # ug/(mlxday)
kT=1.68, # 1/days
kB=6.4,# 1/days
sig=0.72,# kg/day
gG=4.36,# days
k=10, #(ug/ml)^2   #Chemotherapy
del=100, # 1/days
phi=32.07, #1/days
hS=0.0702, #kg/(mgxday)
hNP=0.4275 #kg/(mgxday)
)
