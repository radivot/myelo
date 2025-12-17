# Graphical abstrac plots
library(phaseR)
hahnYZp<-function(Time, State, Pars) {
  with(as.list(c(Time, State, Pars)),{
    dy =  py*y*(1-y/Ky)   -  m*y*z 
    dz =  rz    -   a*z                   + pz*y*z/(Kz^2+y^2)  
    list(c(dy,dz))
  })
}
(pars=c(pz=4e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 400, m = 1e-04))
equi=findEquilibrium(hahnYZp,parameters=pars, state.names=c("y","z"),y0 = c(200, 2e3))
(y0=matrix(c(1, pars["rz"]/pars["a"],equi$ystar[1], 0),ncol=2,byrow=T))

hahnYZ_trajectory <- trajectory(hahnYZp, add=F,
                                parameters=pars,
                                xlim = c(0, 5e3),
                                ylim = c(0, 5e3),
                                state.names=c("y","z"),
                                y0   = y0,
                                tlim = c(0, 100),log="",ylab="z",xlab="y")
dev.copy2pdf(file="~/Results/twoCities/hiroGAsimp.pdf",width=3,height=2.5) 

(pars=c(pz=4.1e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 400, m = 1e-04))
equi=findEquilibrium(hahnYZp,parameters=pars, state.names=c("y","z"),y0 = c(200, 2e3))
(y0=matrix(c(1, pars["rz"]/pars["a"],equi$ystar[1], 0),ncol=2,byrow=T))

hahnYZ_trajectory <- trajectory(hahnYZp, add=F,
                                parameters=pars,
                                xlim = c(0, 5e3),
                                ylim = c(0, 5e3),
                                state.names=c("y","z"),
                                y0   = y0,
                                tlim = c(0, 100),log="",ylab="z",xlab="y")
dev.copy2pdf(file="~/Results/twoCities/nagaGAsim.pdf",width=3,height=2.5) 

