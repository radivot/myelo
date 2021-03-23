library(phaseR)
hahnYZp<-function(Time, State, Pars) {
  with(as.list(c(Time, State, Pars)),{
    dy =  py*y*(1-y/Ky)   -  m*y*z 
    dz =  rz    -   a*z                   + pz*y*z/(Kz^2+y^2)  
    list(c(dy,dz))
  })
}


(pars=c(pz=4.1e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 200, m = 1e-04))
equi=findEquilibrium(hahnYZp,parameters=pars, state.names=c("y","z"),y0 = c(200, 2e3))
# (y0=matrix(c(1, pars["rz"]/pars["a"],equi$ystar[1], 0),ncol=2,byrow=T))
(y0=matrix(c(1, pars["rz"]/pars["a"]),ncol=2,byrow=T))
hahnYZ_flowField  <- flowField(hahnYZp, 
                               parameters=pars,
                               state.names=c("y","z"),
                               xlim = c(0, 5e3), 
                               ylim = c(0, 1e4),
                               add  = FALSE)
grid()
hahnYZ_nullclines <- nullclines(hahnYZp,add.legend=FALSE,
                                parameters=pars,
                                state.names=c("y","z"),
                                xlim = c(0, 5e3), 
                                ylim = c(0, 1e4))

hahnYZ_trajectory <- trajectory(hahnYZp,
                                parameters=pars,
                                state.names=c("y","z"),
                                y0   = y0,
                                tlim = c(0, 100))
title("4.1")
dev.copy2pdf(file="~/Results/twoCities/pzrz4_1.pdf",width=4,height=3.5) 


##################
(pars=c(pz=4.1e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 400, m = 1e-04))
equi=findEquilibrium(hahnYZp,parameters=pars, state.names=c("y","z"),y0 = c(200, 2e3))
(y0=matrix(c(1, pars["rz"]/pars["a"],equi$ystar[1], 0),ncol=2,byrow=T))
hahnYZ_flowField  <- flowField(hahnYZp, 
                               parameters=pars,
                               state.names=c("y","z"),
                               xlim = c(0, 5e3), 
                               ylim = c(0, 1e4),
                               add  = FALSE)
grid()
hahnYZ_nullclines <- nullclines(hahnYZp,add.legend=FALSE,
                                parameters=pars,
                                state.names=c("y","z"),
                                xlim = c(0, 5e3), 
                                ylim = c(0, 1e4))

hahnYZ_trajectory <- trajectory(hahnYZp,
                                parameters=pars,
                                state.names=c("y","z"),
                                y0   = y0,
                                tlim = c(0, 100))
title("4.1x2")
dev.copy2pdf(file="~/Results/twoCities/pzrz4_1x2.pdf",width=4,height=3.5) 


(pars=c(pz=4e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 200, m = 1e-04))
equi=findEquilibrium(hahnYZp,parameters=pars, state.names=c("y","z"),y0 = c(200, 2e3))
# (y0=matrix(c(1, pars["rz"]/pars["a"],equi$ystar[1], 0),ncol=2,byrow=T))
(y0=matrix(c(1, pars["rz"]/pars["a"]),ncol=2,byrow=T))
hahnYZ_flowField  <- flowField(hahnYZp, 
                               parameters=pars,
                               state.names=c("y","z"),
                               xlim = c(0, 5e3), 
                               ylim = c(0, 1e4),
                               add  = FALSE)
grid()
hahnYZ_nullclines <- nullclines(hahnYZp,add.legend=FALSE,
                                parameters=pars,
                                state.names=c("y","z"),
                                xlim = c(0, 5e3), 
                                ylim = c(0, 1e4))

hahnYZ_trajectory <- trajectory(hahnYZp,
                                parameters=pars,
                                state.names=c("y","z"),
                                y0   = y0,
                                tlim = c(0, 100))
title("4")
dev.copy2pdf(file="~/Results/twoCities/pzrz4.pdf",width=4,height=3.5) 


##################
(pars=c(pz=4e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 400, m = 1e-04))
equi=findEquilibrium(hahnYZp,parameters=pars, state.names=c("y","z"),y0 = c(200, 2e3))
(y0=matrix(c(1, pars["rz"]/pars["a"],equi$ystar[1], 0),ncol=2,byrow=T))
hahnYZ_flowField  <- flowField(hahnYZp, 
                               parameters=pars,
                               state.names=c("y","z"),
                               xlim = c(0, 5e3), 
                               ylim = c(0, 1e4),
                               add  = FALSE)
grid()
hahnYZ_nullclines <- nullclines(hahnYZp,add.legend=FALSE,
                                parameters=pars,
                                state.names=c("y","z"),
                                xlim = c(0, 5e3), 
                                ylim = c(0, 1e4))

hahnYZ_trajectory <- trajectory(hahnYZp,
                                parameters=pars,
                                state.names=c("y","z"),
                                y0   = y0,
                                tlim = c(0, 100))
title("4x2")
dev.copy2pdf(file="~/Results/twoCities/pzrz4x2.pdf",width=4,height=3.5) 



