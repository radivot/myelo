# These values were all obtained from the TBMM 2012 supplement 

cells=c(Snor= 1, # normal value of stem cells set
		TS=8, # duration of cell cycle  WL*
		p=0.1, # self-renewal probability  WL
		aminS=0.01, # proliferative fraction under minimal stimulation  WL
		anorS=0.15, # proliferative fraction under normal stimulation  WL
		aintS=0.45, # proliferative fraction under intensified stimulation  WL
		amaxS=1, # proliferative fraction under maximal stimulation WL
		wG=0.4, # weighting parameter G for regulation of a  WL
		wS=1, # weighting parameter S for regulation of a  WL
		thetaG=-10, # weighting parameter G for regulation of p  WL
		aminCG=0.1205, # proliferative fraction under minimal stimulation  fitted
		anorCG=0.1252, # proliferative fraction under normal stimulation  fitted
		aintCG=0.8340, # proliferative fraction under intensified stimulation  fitted
		amaxCG=1, # proliferative fraction under maximal stimulation  WL*
		NG4=5,# number of subcompartments in G4  set
		NG5=5,# number of subcompartments in G5  set
		NG6=5,# number of subcompartments in G6 set
		TnorGRA=5.576, # transition time of granulocytes  fitted, [4]
		TpredGRA=0.466)# prolongation of TnorGRA under Prednisone  fitted, [5]

gProd=c(
		Pendomax=257, # maximal G-CSF production  fitted
		Pendonor=1, #normal G-CSF production set
		Pendomin=0.318, # minimal G-CSF production  fitted
		Pendob=0.022, # sensitivity parameter of G-CSF production  fitted
		wGRA= 1, # influence o	f GRA on G-CSF production set, [1]
		wG6=0.2) # influence of G6 on G-CSF production  set, [1]


gPK=c(
		kFsc=0.161, # subcutaneous absorption  fitted
		kFm=34.7, #Michaelis-Menten constant of subcutaneous elimination  fitted
		vFmax=67.3, #Maximum of subcutaneous elimination  fitted
		kFu=0.441, # unspecific elimination  fitted
		kFcp=0.000, # transition central to peripheral  fitted
		kFpc=NA, # transition peripheral to central - not determinable
		VFD=1.156, # distribution volume  fitted
		vGRAFmax=4.77, # Maximum of specific elimination  fitted
		kGRAFm=22.4 #Michaelis-Menten constant of specific elimination  fitted
) # steady state natural g-csf = 0.02 ug/l = (.02ug/L)(nmol/20ug)= 1 pM

pegPK=c(
		kPsc=0.107, # subcutaneous absorption  fitted
		kPm=5.5,# Michaelis-Menten constant of subcutaneous elimination  fitted
		vPmax=16.5, # Maximum of subcutaneous elimination  fitted
		kPu=0.087, # unspecific elimination  fitted
		kPcp=0.075, # transition central to peripheral  fitted
		kPpc=0.548, # transition peripheral to central  fitted
		VPD=4.091, # distribution volume  fitted
		vGRAPmax=5.16, # Maximum of specific elimination  fitted
		kGRAPm=30.8) # Michaelis-Menten constant of specific elimination  fitted
pegPK

gPD=c(
		AminCGF=0.910, # amplification in CG under minimal stimulation  fitted
		AnorCGF=105, # amplification in CG under normal stimulation  fitted
		AmaxCGF=206,# amplification in CG under maximal stimulation  fitted
		AbCGF=0.024, # sensitivity of amplification in CG  fitted
		TminCGF=47.3, #CGF transition time in CG under minimal stimulation  fitted
		TnorCGF=78.2, # transition time in CG under normal stimulation  fitted
		TmaxCGF=286, #transition time in CG under maximal stimulation  fitted
		TbCGF=0.590, #sensitivity of transition time in CG  fitted
		AminPGBF=1.31, # amplification in PGB under minimal stimulation  fitted
		AnorPGBF=61.2, # amplification in PGB under normal stimulation  fitted
		AmaxPGBF=815, # amplification in PGB under maximal stimulation  fitted
		AbPGBF=0.721, # sensitivity of amplification in PGB  fitted
		TminPGBF=4.64, # transition time in PGB under minimal stimulation  fitted
		TnorPGBF=40.9, # transition time in PGB under normal stimulation  fitted
		TmaxPGBF=217, # transition time in PGB under maximal stimulation  fitted
		TbPGBF=0.104, # sensitivity of transition time in PGB  fitted
		AnorG4F=1, # postmitotic amplification in G4 set
		TminG4F= 119, # transition time in G4 under minimal stimulation fitted
		TnorG4F=11.4, # transition time in G4 under normal stimulation  fitted
		TmaxG4F=3.93, # transition time in G4 under maximal stimulation  fitted
		TbG4F=0.366, # sensitivity of transition time in G4  fitted
		AnorG5F=1, # postmitotic amplification in G5  set
		TminG5F=48.3, #transition time in G5 under minimal stimulation  fitted
		TnorG5F=37.0, # transition time in G5 under normal stimulation  fitted
		TmaxG5F=4.64, # transition time in G5 under maximal stimulation  fitted
		TbG5F=0.459, # sensitivity of transition time in G5  fitted
		AminG6F=0.201, #postmitotic amplification in G6 under minimal stimulation  fitted
		AnorG6F=0.249, # postmitotic amplification in G6 under normal stimulation  fitted
		AmaxG6F=0.850, # postmitotic amplification in G6 under maximal stimulation  fitted
		AbG6F=0.503, # sensitivity of postmitotic amplification in G6  fitted
		TminG6F=141, # transition time in G6 under minimal stimulation  fitted
		TnorG6F=82.0, # transition time in G6 under normal stimulation  fitted
		TmaxG6F=41.4, # transition time in G6 under maximal stimulation  fitted
		TbG6F=0.526, # sensitivity of transition time in G6  fitted
		DFGCSF=1.08, # Delay of Filgrastim action  fitted
		NFGCSF=4) # number of delay compartments  set

pegPD=c(
		AminCGP=4.21, # amplification in CG under minimal stimulation  fitted
		AnorCGP=357, #amplification in CG under normal stimulation  fitted
		AmaxCGP=669, #amplification in CG under maximal stimulation  fitted
		AbCGP=0.103, # sensitivity of amplification in CG  fitted
		TminCGP=0.848, #transition time in CG under minimal stimulation  fitted
		TnorCGP=1.69, # transition time in CG under normal stimulation  fitted
		TmaxCGP=57.6, # transition time in CG under maximal stimulation  fitted
		TbCGP=0.040, # sensitivity of transition time in CG  fitted
		
		AminPGBP=0.067, # amplification in PGB under minimal stimulation  fitted
		AnorPGBP=24.5, # amplification in PGB under normal stimulation  fitted
		AmaxPGBP=47.0, # amplification in PGB under maximal stimulation  fitted
		AbPGBP=0.327, # sensitivity of amplification in PGB  fitted
		TminPGBP=0.756, #transition time in PGB under minimal stimulation  fitted
		TnorPGBP=24.7, # transition time in PGB under normal stimulation  fitted
		TmaxPGBP=24.7, # transition time in PGB under maximal stimulation  fitted
		TbPGBP=0.118, # sensitivity of transition time in PGB  fitted
		
		AnorG4P=1, # postmitotic amplification in G4 set
		TminG4P=153,# transition time in G4 under minimal stimulation  fitted
		TnorG4P=14.2,# transition time in G4 under normal stimulation  fitted
		TmaxG4P=4.17,# transition time in G4 under maximal stimulation  fitted
		TbG4P=0.225, # sensitivity of transition time in G4  fitted
		
		AnorG5P=1, # postmitotic amplification in G5 set
		TminG5P=153,# transition time in G5 under minimal stimulation  fitted
		TnorG5P=14.2,# transition time in G5 under normal stimulation  fitted
		TmaxG5P=4.17,# transition time in G5 under maximal stimulation  fitted
		TbG5P=0.225, # sensitivity of transition time in G5  fitted
		
		AminG6P=0.139, # postmitotic amplification in G6 under minimal stimulation  fitted
		AnorG6P=0.368, #postmitotic amplification in G6 under normal stimulation fitted
		AmaxG6P=1, # postmitotic amplification in G6 under maximal stimulation set
		AbG6P=0.321, # sensitivity of postmitotic amplification in G6  fitted
		TminG6P=153,# transition time in G6 under minimal stimulation  fitted
		TnorG6P=14.2,# transition time in G6 under normal stimulation  fitted
		TmaxG6P=4.17,# transition time in G6 under maximal stimulation  fitted
		TbG6P=0.225, # sensitivity of transition time in G6  fitted
		
		dPgcsf=0.984, #Delay of Pegfilgrastim action  fitted
		nPgcsf=4, # number of delay compartments  set
		wPmin=0, # minimum of weighting function set
		wPnor=0.499, # value of weighting function for 1ug Pegfilgrastim  fitted
		wPmax=1, # maximum of weighting function set
		wPb=0.068) #sensitivity parameter of weighting function fitted
pegPD


chop=c(
		ffc 		=1.09, # first cycle effect fitted
		DCX 		=0.0486, # delay of toxicity fitted 
		kS   		=0.212,  #toxicity in S compartment  fitted
		kCG  		=0.464,  #toxicity in CG compartment fitted
		kPGB 		=0.168, #toxicity in PGB compartment  fitted
		kMGB 		=0.000148, #toxicity in MGB compartment  fitted
		DCXLYM 	=0.0223, #delay of lymphocyte toxicity  fitted
		kLYM 		=16.4)  #lymphocyte toxicity 
chop


scholzPars=c(cells,gProd,gPK,gPD,pegPK,pegPD,chop)
scholzPars
length(scholzPars)
length(unique(names(scholzPars)))
save(scholzPars,file="scholzPars.rda")

