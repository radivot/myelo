#Grace Brooks, Gabriel Provencher, Jinzhi Lei, Michael C. Mackey
#Neutrophil dynamics after chemotherapy and G-CSF: The role
#of pharmacokinetics in shaping the response
#Journal of TheoreticalBiology 315(2012) 97–109

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
bn=0.05, # mg/ml    # not indented = param not in Zhuge et al
		tauNP=5.9, # days       Note: many values changed
		tauNM=3.8, # days
		tauN=9.7, # days
bv=0.001, # ug/ml 
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
gamG=4.36,# days
k=10, #(ug/ml)^2   #Chemotherapy
del=100, # 1/days
phi=32.07, #1/days
hS=0.0702, #kg/(mgxday)
hNP=0.4275 #kg/(mgxday)
)

save(brooksPars,file="brooksPars.rda")
