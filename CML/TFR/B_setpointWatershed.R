graphics.off();rm(list=ls())
library(pracma)
F <- function(p) {
  pz=p[1]; Kz=p[2]
  az=2;rz=200;m=0.0001; py=1.658
  a=-(rz - az*(py/m))
  b= -pz*(py/m)
  c = -Kz^2*(rz - az*(py/m))
  r <- rep(NA, length(p))
  r[1] <- -SP + (-b - sqrt(b^2-4*a*c))/(2*a) 
  r[2] <- -WS + (-b + sqrt(b^2-4*a*c))/(2*a) 
  return(r)
}
SP=100*2 ####  hold setpoint at MR4
WS=1000*2 ## and watershed at MR3
p0 <- c(pz=3000, Kz=500) #initial guess
(f=fsolve(F, p0)) 
cat("pz=",f$x[1],"Kz=",f$x[2]) #pz= 4373.462 Kz= 632.4555
