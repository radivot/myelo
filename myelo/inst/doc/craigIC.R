# M. Craig A. R. Humphries M. C. Mackey
#'Granulopoiesis Incorporating the Negative Feedback Dynamics and Kinetics of 
#'G-CSF/Neutrophil Binding and Internalization, Bull Math Biol (2016) 78:2304-2357  
# Initial Conditions. Note that G1ss is both an IC and a parameter in the model.
craigIC=c(     #Similary tauNP is both a parameter and a component of  Tn IC = tauNP + aNM 
  Q=1.1, # 1e6 cells/kg   #define Qss parms[0]
  Nr = 2.26, # 1e9 cells/kg 
  N = 0.3761, # =0.22/0.585  1e9 cells/kg  
  G1= 0.025, # ng/mL # was Gss = 0.0246,  # ng/mL
  G2 = 2.1899e-5, # ng/mL
  Tn = 7.3074 + 3.9,      # tauNP + aNM = days
  An = 1.0378e5, # 103780,  #, text 103777     
  Aq = 1.5116, # = 2*exp(-0.1*2.8), text 1.512 
  Cp = 0,   
  Cf = 0,    
  Cs1= 0,    
  Cs2= 0,    
  Gs = 0,
  Ic = 0,
  Ig = 0) 
save(craigIC,file="craigIC.rda")
library(tools)
showNonASCIIfile("R/myelo-package.R")
# library(tools)
# showNonASCIIfile("R/craig16.R")

