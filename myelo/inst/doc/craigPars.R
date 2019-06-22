# M. Craig A. R. Humphries M. C. Mackey
#'Granulopoiesis Incorporating the Negative Feedback Dynamics and Kinetics of 
#'G-CSF/Neutrophil Binding and Internalization, 'Bull Math Biol (2016) 78:2304-2357  
# short verion. 44 actually used, out of 75 in craigPars16 
craigPars=c(
  Qss=1.1, # 1e6 cells/kg   #define Qss parms[0]
  Nrss = 2.26, # 1e9 cells/kg 
  Nss = 0.3761, # =0.22/0.585  1e9 cells/kg  
  G1ss=0.025, # ng/mL # was Gss = 0.0246,  # ng/mL
  G2ss = 2.1899e-5, # ng/mL
  Aqss = 1.5116, # = 2*exp(-0.1*2.8), text 1.512 
  Anss = 1.0378e5, # 103780,  #, text 103777     
  gamNMss = 0.1577,# 1/day
  gamNrss = 0.0063661, # 0.0064, # 1/day
  gamNss= 2.1875, # 1/day  (= 35/16)
  tauS=2.8, #days
  tauNP = 7.3074,    # days

  fQ = 8,      # 1/day  #define fQ   parms[12]
  the2 = 0.0809, # 1e6 cells/kg
  s2 = 2,
  kapss = 0.0073325, # 1/day
  kapDel = 0.0146, # 1/day  #define kapDel parms[16]
  
  etaNPss = 1.6647, # 1/day  parms[17]  //growth rates
  etaNPmin = 1.4060, ####### was 0.4,    # 1/day
  bNP = 0.022868, ####### was 11.2679, # ng/mL #define bNP parms[19]
  
  ka=5.143,  #define ka parms[20]
  kren = 0.16139, ########  10.3,  # 1/day  
  kint = 462.42,  ######## 114.48, # 1/day
  Gprod = 0.014161, #### was 0.2535,  # ng/mL/day
  k12g = 2.2423, #(ng/mL)-Pow days-1
  k21g = 184.87, #1/day
  Pow = 1.4608,  #define Pow parms[26]
  
  Vmax = 7.8670, #was 10,         # max maturation velocity
  bV = 0.24611, ######  was 3.5,       # ng/mL  #define bV parms[28]
  
  phiNrss = 0.3640,  # 1/day (text 0.387)   # ftransss = 0.3640,  # 1/day (text 0.387)
  phiNrmax=4.1335,                 ########was transmax = 1.456,  # 1/day
  bG = 1.8924e-4, ###### was 11.2679  # ng/mL
  V=0.525, # bound GCSF conversion factor  #define V parms[32]
  k12 = 90.2752,# 1/day
  k21 = 18.2222,# 1/day  map 1=p, 2=f, 3=sl1, 4=sl2  (don't bother rewriting) 
  k24 = 9.2296,# 1/day
  k42 = 62.5607,# 1/day
  k13 = 8.2936,# 1/day
  k31 = 0.6990,# 1/day
  kelC = 132.0734, # 1/day  #define kelC parms[39]
  
  hQ=0.0079657,  # was hS = 0.1, (new is Value 2 in Table 3)      #define hQ parms[40]
  EM50 = 0.75390*32.7, #map to mass in ug in V1 to match units of chemo ODE states
  # see table 3 of Cancer Chemother Pharmacol (2012) 69:15–24 to find 32.7
  Sc = 0.89816 # was h=3 in (Cp/EC50)^h in 2015, now (Cp/EC50)^Sc in 2016  #define Sc parms[42]
)
save(craigPars,file="craigPars.rda")
library(tools)
showNonASCIIfile("R/myelo-package.R")
# library(tools)
# showNonASCIIfile("R/craig16.R")

