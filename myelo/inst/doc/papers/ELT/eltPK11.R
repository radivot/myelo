library(myelo)  
library(deSolve)
library(rodeoExt)
# Pharmacokinetics and Hepatic Uptake of Eltrombopag, a Novel
# Platelet-Increasing Agent,   S Kazuya Takeuchi, Tomoko Sugiura, Saki Umeda,
# Kazuki Matsubara, Masato Horikawa, Noritaka Nakamichi, David L. Silver,
# Norihisa Ishiwata, and Yukio Kato, DRUG METABOLISM AND DISPOSITION Copyright
# 2011 by The American Society for Pharmacology and Experimental Therapeutics
# DMD 39:1088–1096, 2011    # ELT has a COO- group, and rifampicin blocks anion transport into liver

#IV in rats, 40% out in bile unchanged at 72 hours

# table 1
CLtot=33.3  #ml/(kg*h)  
CLbile=13.6  #ml/(kg*h)  [0 to 72 hours]
CLurine=0.07  #ml/(kg*h) [at most]
V0=41.1 #ml/kg
Vdss=287 #ml/kg
#table 2
CLbile=2.38  #ml/(kg*h) [0 to 8 hours]
CLtot=34.3  #ml/(kg*h)
CLuptake=44.8  #ml/(kg*h)

