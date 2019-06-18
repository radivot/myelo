# Morgan Craig ... Michael C. Mackey
# Neutrophil dynamics during concurrent chemotherapy and G-CSF administration: 
# Mathematical modelling guides dose optimisation to minimise neutropenia
# Journal of Theoretical Biology 385 (2015) 77-89

craigPars15=c(
  Qss=1.1, # 1e6 cells/kg
  gamSss=0.1, # 1/day
  tauS=2.8, #days
  AQss = 1.5116, # = 2*exp(-0.1*2.8), text 1.512 
  kapDel = 0.0140, # 1/day
  betaQss = 0.043, # 1/day
  fQ = 8,      # 1/day
  s2 = 2,
  the2 = 0.0809, # 1e6 cells/kg
  #
  Nrss = 2.26, # 1e9 cells/kg 
  Nss = 0.3761, # 1e9 cells/kg 
  Ncircss = 0.22, # 1e9 cells/kg  
  
  gamNss=2.1875, # 1/day
  tauNP = 7.3074,    # days
  aNM = 3.9,      # days
  tauNr = 2.7,    # days 2.5
  gamNr = 0.0064, # 1/day
  gamNMss = 0.1577,# 1/day
  kapNss = 0.0073, # 1/day
  ANss = 103780,  #, text 103777     exp(etaNPss*tauNP-gamNr*tauNr-gamNM*aNM) =>
  etaNPss = 1.6647, # 1/day         exp(1.6647*7.3074-0.0064*2.7-0.1577*3.9)=101958.3
  fN = 0.0088,    # 1/day
  s1 = 2,
  the1 = 0.8409,  # 1e9 cells/kg  
  ftransss = 0.3640,  # 1/day (text 0.387)
  
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
  chi = 0.0654,  # (ng/mL)/(1e9 cells/kg)  # 0.0246/0.3761
  kint = 114.48, # 1/day
  Kd = 1.44,   # ng/mL
  ka = 13.5,  # 1/day
  F = 0.6020,
  Vd = 1788,   #mL
  
  gamSss = 0.1,  # 1/day
  gamSmin = 0.1, # 1/day
  gamSmax = 0.4, # 1/day  #not used????
  hS = 0.1,      
  bS = 11.2679, # ng/mL
  h = 3,
  EC50 = 2.3056, # ng/mL
  
  etaNPmax = 2.544,  # 1/day
  etaNPmin = 0.4,    # 1/day
  Vmax = 10,         # max maturation velocity
  gamNMmin = 0.12,   # 1/day
  gamNMmax = 0.67,   # 1/day
  transmax = 1.456,  # 1/day
  bV = 3.5,       # ng/mL
  bNP = 11.2679, # ng/mL
  bNM = 11.2679, # ng/mL
  bG = 11.2679  # ng/mL
)
save(craigPars15,file="craigPars15.rda")
# library(tools)
# showNonASCIIfile("R/myelo-package.R")
