#Grace Brooks, Gabriel Provencher, Jinzhi Lei, Michael C. Mackey
#Neutrophil dynamics after chemotherapy and G-CSF: The role
#of pharmacokinetics in shaping the response
#Journal of TheoreticalBiology 315(2012) 97–109

brooksPars=c(
# Stem cell params
		Qss=1120693, #cells/kg
		gamS=0.1043, # 1/day
		gamMinS=0.03,
		gamMaxS=0.40,   # NOT USED  ###############
		tauS=2.83, #days1,2
		k0=8.0, # days
		the2=0.0826e6, # cells/kg  the=theta
		s2=2,
# Neutrophil params
		Nss=559290277, # cells/kg
		gamN=2.4, # 1/day
		bn=0.05, # mg/ml NOTE: This is b0 in Eq. 16  
		tauNP=5.9, # days       
		tauNMmax=3.8, # days
		tauN=9.7, # days
		bv=0.001, # ug/ml 
		Vmax=3.8,
		etaNP=2.1995, # 1/days
		etaMinNP=0.4, # 1/days
		etaMaxNP=2.5444, # 1/days
		AN=1549.58e2,   # NOT USED, COMPUTED FROM OTHER PARAMS
		gam0=0.27, # 1/day
		gamMin0=0.12, 
		f0=0.154605, # 1/days
		the1=0.0154848e8, # cells/kg 
		s1=.5, 
		kdel=0.0134, # 1/days   increased >10-fold in Figure 5 to destabilize
# G-CSF params
		Xss=0.1, # ug/kg   SHOULDN'T THIS BE ZERO?????
		Gss=0, # ug/ml 
		VB=76, # ml/kg
		Gprod=0, # ug/(mlxday)
		kT=1.68, # 1/days
		kB=6.4,# 1/days
		sig=0.72,# kg/day
		gamG=4.36,# days
		kG=10, #(ug/ml)^2   
		T1=5, # figure 9 shows this is the best time for GCSF in days after chemo starts
		Delg=1, # spread GCSF injection of 1 day to help numerics (i.e. need samples at higher res
		Dg=0,  # default is no GCSF, as in Figures 3 and 5
#Chemotherapy params
		del=100, # 1/days
		phi=32.07, #1/days
		hS=0.0702, #kg/(mgxday)
		hNP=0.4275, #kg/(mgxday)
		cycles=1, # stop after one cycle of therapy
		Delc=1, # spread chemo injection over 1 day to help numerics
		Dc=135, # chemo/taxol dose of 135 mg/kg 
		T=21, # default cycle length = 21 days
		bS=0.01, # not in pdf, but in C++ code
		cn=0.085 # same here
)

save(brooksPars,file="brooksPars.rda")

