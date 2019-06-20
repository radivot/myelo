# Morgan Craig ... Michael C. Mackey
# M. Craig A. R. Humphries M. C. Mackey
#'Granulopoiesis Incorporating the Negative Feedback Dynamics and Kinetics of 
#'G-CSF/Neutrophil Binding and Internalization, 
#'Bull Math Biol (2016) 78:2304-2357  

craigPars16=c(
  # table 4
  gamSss=0.1, # 1/day
  tauS=2.8, #days
  AQss = 1.5116, # = 2*exp(-0.1*2.8), text 1.512 
  fQ = 8,      # 1/day
  s2 = 2,
  the2 = 0.0809, # 1e6 cells/kg
  
  kapDel = 0.0146, # 1/day
  kapss = 0.0073325, # 1/day
  kapmin = 0.0073325, # 1/day
  
  s1 = 1.5,  ######## was 2,
  etaNPss = 1.6647, # 1/day         exp(1.6647*7.3074-0.0064*2.7-0.1577*3.9)=101958.3
  bNP = 0.022868, ####### was 11.2679, # ng/mL
  
  etaNPmin = 1.4060, ####### was 0.4,    # 1/day
  # etaNPmax = 2.544,  # 1/day (not used?)
  tauNP = 7.3074,    # days
  
  Vmax = 7.8670, #was 10,         # max maturation velocity
  bV = 0.24611, ######  was 3.5,       # ng/mL

  aNM = 3.9,      # days
  gamNMss = 0.1577,# 1/day
  
  phiNrss = 0.3640,  # 1/day (text 0.387)   # ftransss = 0.3640,  # 1/day (text 0.387)
  phiNrmax=4.1335,                 ########was transmax = 1.456,  # 1/day
  
  bG = 1.8924e-4, ###### was 11.2679  # ng/mL
  gamNr = 0.0063661, # 0.0064, # 1/day
  gamNss= 2.1875, # 1/day  (= 35/16)
  G1ss=0.025, # ng/mL # was Gss = 0.0246,  # ng/mL
  GBFss=1.5823e-5, 
  Gprod = 0.014161, #### was 0.2535,  # ng/mL/day
  V=0.525, # bound GCSF conversion factor
 
  # all new vals in this chunk match value 4 (0.8 at top) in Table 1
  kren = 0.16139, ########  10.3,  # 1/day  
  kint = 462.42,  ######## 114.48, # 1/day
  k12g = 2.2423, #(ng/mL)-Pow days-1
  k21g = 184.87, #1/day
  Pow = 1.4608,

  # table 5  # not used explicitle but used to make other parameters
  Qss=1.1, # 1e6 cells/kg
  betaQss = 0.043, # 1/day
  Nss = 0.3761, # =0.22/0.585  1e9 cells/kg  
  Ncircss = 0.22, # 1e9 cells/kg  
  Nrss = 2.26, # 1e9 cells/kg 
  Npss = 0.93, # 1e9 cells/kg 
  Nmss = 4.51, # 1e9 cells/kg 
  G2ss = 2.1899e-5, # ng/mL
  tauNr = 2.7,    # days was 2.5
  tauNcircss = 0.4571429,    # days =16/35
  tauhalf = 7.6,    # hrs 
  ANss = 1.0378e5, # 103780,  #, text 103777     exp(etaNPss*tauNP-gamNr*tauNr-gamNM*aNM) =>
  bbarV=0.031283, #matches value 4 (0.8 at top) in table 2
  phiRatio=11.356,  # max/ss
  phiMin=0.020056, 
  theta=0.15096, # ratio in GCSF KO mice (uses other cytokines)
  Cko=0.25,  # N fraction in KO (relative to WT? or other WBC?)
  mu=0.84458, # min/ss prolif rates  (value 4 in table 2, 0.8 at top matches table 3)

  ## Table 6 PK params
  Vd300=4754.7,
  F300=0.64466,
  ka300=8.0236,  
  Vd375=2322.9,
  F375=0.49964,
  ka375=6.6133,  
  Vd750=2178.0,
  F750=0.75,
  ka750=5.143,  
  # set defaults to 750 (override if needed)
  Vd=2178.0,
  F=0.75,
  ka=5.143,  
  # 2015 JTB had
  # Vd = 1788,   #mL
  # F = 0.6020,
  # ka = 13.5,  # 1/day
  k21 = 18.2222,# 1/day  map 1=p, 2=f, 3=sl1, 4=sl2  (don't bother rewriting) 
  k31 = 0.6990,# 1/day
  k12 = 90.2752,# 1/day
  k13 = 8.2936,# 1/day
  kelC = 132.0734, # 1/day
  k42 = 62.5607,# 1/day
  k24 = 9.2296,# 1/day
  BSA = 1.723, # m^2
  
  hQ=0.0079657,  # was hS = 0.1, (new is Value 2 in Table 3)    
  EC50 = 0.75390, #was 2.3056, # ng/mL (new is not in Table 3, but close to Value 2)
  EM50 = 0.75390*32.7, #map to mass in ug in V1 to match units of chemo ODE states
  # see table 3 of Cancer Chemother Pharmacol (2012) 69:15–24 to find 32.7
  Sc = 0.89816, # was h=3 in (Cp/EC50)^h in 2015, now (Cp/EC50)^Sc in 2016
  etaNPInf=0  # chemo at Inf shuts down P completely 
  # # below was in 2015 but not in 2016 
  # bS = 11.2679, # ng/mL
  # Kd = 1.44,   # ng/mL
  # chi = 0.0654,  # (ng/mL)/(1e9 cells/kg)  # 0.0246/0.3761
  # kapNss = 0.0073, # 1/day
  # fN = 0.0088,    # 1/day
  # the1 = 0.8409,  # 1e9 cells/kg  
  # gamSss = 0.1,  # 1/day
  # gamSmin = 0.1, # 1/day
  # gamSmax = 0.4, # 1/day  #not used????
  # gamNMmin = 0.12,   # 1/day
  # gamNMmax = 0.67,   # 1/day
  # bNM = 11.2679, # ng/mL
)
save(craigPars16,file="craigPars16.rda")
library(tools)
showNonASCIIfile("R/myelo-package.R")
# library(tools)
# showNonASCIIfile("R/craig16.R")

