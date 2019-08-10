#Lena E. Friberg, Anja Henningsson, Hugo Maas, Laurent Nguyen, and Mats O. Karlsson
#Model of Chemotherapy-Induced Myelosuppression With Parameter Consistency Across Drugs, 
#\emph{Journal of Clinical Oncology},  \bold{20},  4713-4721 (2002)

fribergPars02=c(
		Circ0=5.05, # 1e9 cells/L
		# MTT=88.7, # hours
		ktr=24*4/88.7, # (n+1)/MTT where n is number of transit compartments
		gam=0.161, # unitless power
		slope=8.58, # (1/uM)
		# Emax=83.9,
		# EC50=7.17,
		### below are from Bruno et al J.  Pharmacokinetics and Biopharmaceutics, Vol.24, No.2, 1996
		k12 = 24*1.06,# 1/hr  put in 1/days
		k21 = 24*1.51, 
		k13 = 24*1.26, 
		k31 = 24*0.084, 
		k10 = 24*5.2,# CL=38.5 L/hr , V1=7.4L  38.5/7.4=5.2
		V1=7.4,
		mw=0.808 
)
save(fribergPars02,file="fribergPars02.rda")

