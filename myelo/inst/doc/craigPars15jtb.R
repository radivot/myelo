# Morgan Craig ... Michael C. Mackey
# Neutrophil dynamics during concurrent chemotherapy and G-CSF administration: 
# Mathematical modelling guides dose optimisation to minimise neutropenia
# Journal of Theoretical Biology 385 (2015) 77-89

craigPars15jtb=c(
  Qss=1.1, # 1e6 cells/kg
  gamS=0.1, # 1/day
  tauS=2.8, #days
  AQss = 1.5116,
  kapDelta = 0.0140, # 1/day
  betaQss = 0.043, # 1/day
  fQ = 8,      # 1/day
  s2 = 2,
  the2 = 0.0809, # 1e6 cells/kg
  #
  Nrhomeo = 2.26, # 1e9 cells/kg 
  Nhomeo = 0.3761, # 1e9 cells/kg 
  Nhomeocirc = 0.22, # 1e9 cells/kg  
  
  gamN=2.1875, # 1/day
  tauNM=3.5, # days
  An=752, # none  # 70635 in Z19, 
  f0=0.40, # 1/days
  the1=0.36, # 1e8 cells/kg 
  n=1, # n=1

  tauNP = 7.3074,    # days
  aNM = 3.9,      # days
  tauNr = 2.7,    # days
  gamNr = 0.0064, # 1/day
  gamNM = 0.1577,# 1/day
  kapNss = 0.0073, # 1/day
  ANss = 103780,
  etaNPss = 1.6647, # 1/day
  fN = 0.0088,    # 1/day
  s1 = 2,
  the1 = 0.8409,  # 1e9 cells/kg  
  fhomeotrans = 0.3640,  # 1/day
  
  kelC = 132.0734, # 1/day
  k12 = 90.2752,# 1/day
  k21 = 18.2222,# 1/day
  k13 = 8.2936,# 1/day
  k31 = 0.6990,# 1/day
  k24 = 9.2296,# 1/day
  k42 = 62.5607,# 1/day
  BSA = 1.723, # m^2
  
  Gss = 0.0246,  # ng/mL
  Gprod = 0.2535,  # ng/mL/day
  kren = 10.3,  # 1/day
  chi = 0.0654,  # (ng/mL)/(1e9 cells/kg)
  kint = 114.48, # 1/day
  Kd = 1.44,   # ng/mL
  ka = 13.5,  # 1/day
  F = 0.6020,
  Vd = 1788,   #mL
  
  gamSss = 0.1,  # 1/day
  gamSmin = 0.1, # 1/day
  gamSmax = 0.4, # 1/day
  hS = 0.1,      
  bS = 11.2679, # ng/mL
  h = 3,
  EC50 = 2.3056, # ng/mL
  
  etaMaxNP = 2.544,  # 1/day
  etaMinNP = 0.4,    # 1/day
  Vmax = 10,         # max maturation velocity
  gamMinNM = 0.12,   # 1/day
  gamMaxNM = 0.67,   # 1/day
  transmax = 1.456,  # 1/day
  bV = 3.5,       # ng/mL
  bNP = 11.2679, # ng/mL
  bNM = 11.2679, # ng/mL
  bG = 11.2679  # ng/mL
)
save(craigPars15jtb,file="craigPars15jtb.rda")
# library(tools)
# showNonASCIIfile("R/myelo-package.R")
