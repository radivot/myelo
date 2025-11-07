## Fig5A_1mkAges.R        Step 1 of 3 before running Fig5.R  
# Use single age US mortality data to get mean ages of age groups 0 to 19 in SEER mortality
setwd("~/ccf/CMLR/seer25")
# SEERaBomb::mkMrtLocal() #converts HMD files in ~/data/hmd_countries/USA into ~/data/mrt/mrtUSA.RData
load("~/data/mrt/mrtUSA.RData") #Human Mortality Database (HMD) goes 1933 to 2023 (91 years)
# mrt loaded above is a list of two sexes, each with a matrix of 111 age rows (0 to 110+) and 91 year columns
x=mrt[[1]]; dim(x) 
x[c(1:3,106:111),80:91]
years=as.character(1975:2023)
ages=as.character(0:19) # new age group of 85-89 added in April 2025 release, so now up to 19 instead of 18 
aN=length(ages)
yN=length(years)
(Z=matrix(rep(0,yN*aN),nrow=aN,dimnames=list(ages,years)))
Ages=NULL
for (sex in c("Male","Female")){
  (M=mrt[[sex]][,years])
  (A=Z) # replace zeros here with Age midpoints below
  for (y in years) {
    (m=M[,y])
    age=rep(0,19)
    (age[1]=exp(-m[1])/2)  # first age group is 1 year
    (age[2]=1+4*exp(-sum(m[2:5]))/2) # next is 4 years
    for (i in 3:19) {                # all rest are 5 years
      strt=5*(i-2) # starting point
      mul=1
      tmp=0
      for (j in 1:5) {
        mul=mul*exp(-m[(5*(i-2)+j)]) #prob of alive at end of 1 Y
        tmp=tmp+mul/2
      }
      age[i]=strt+tmp
    }
    strt=90 # starting point
    mul=1
    tmp=0
    for (j in 1:15) { # use HMD data out to ages of 105 years
      mul=mul*exp(-m[(90+j)]) #prob of alive at end of one more year
      tmp=tmp+mul/2
    }
    age[20]=strt+tmp
    A[1:20,y]=age
  }
  Ages[[sex]]=A
}
save(Ages,file="~/data/mrt/Ages.RData")  ####### used by step 2 file Fig5_2mkMrts.R
