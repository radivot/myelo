# http://tbb.bio.uu.nl/rdb/grindR/tutorial.pdf
source("~/ccf/hobbs/methods/bifurc/grind/grind.R")
model <- function(t, state, parms) {
  with(as.list(c(state, parms)),{
    dy =  py*y*(1-y/Ky)   -  m*y*z 
    dz =  rz    -   a*z                   + pz*y*z/(Kz^2+y^2)  
    list(c(dy,dz))
  })
}  
s <- c(y=1,z=100)     

p=c(pz=4.1e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 200, m = 1e-04)
plane(xmin=1,ymin=1,xmax=1e7, ymax=1e5,log="xy",tstep=0.5,portrait=T,main="4.1_200",grid=9)
dev.copy2pdf(file="~/Results/twoCities/phaseDeBoer4p1rz200.pdf",width=5,height=4) 

par(mar=c(4.1,4.1,1.5,0.6))
p=c(pz=4.1e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 400, m = 1e-04)
plane(xmin=1,ymin=1,xmax=1e7, ymax=1e5,log="xy",tstep=0.5,portrait=T,main="pz = 4100   rz = 400",grid=9)
dev.copy2pdf(file="~/Results/twoCities/phaseDeBoer4p1rz400.pdf",width=4,height=3) 

p=c(pz=4.0e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 200, m = 1e-04)
plane(xmin=1,ymin=1,xmax=1e7, ymax=1e5,log="xy",tstep=0.5,portrait=T,main="pz = 4000   rz = 200",grid=9)
dev.copy2pdf(file="~/Results/twoCities/phaseDeBoer4p0rz200.pdf",width=4,height=3) 

p=c(pz=4.0e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 400, m = 1e-04)
plane(xmin=1,ymin=1,xmax=1e7, ymax=1e5,log="xy",tstep=0.5,portrait=T,main="4.0_400",grid=9)
dev.copy2pdf(file="~/Results/twoCities/phaseDeBoer4p0rz400.pdf",width=5,height=4) 



p=c(pz=4.0e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 200, m = 1e-04)
f=newton(c(y=615.8385,z=2498.4604))
continue(f,x="pz",xmax=1e4,y="y",step=0.1,ymax=5e3,ymin=0,log="",main="rz = 200") 
dev.copy2pdf(file="~/Results/twoCities/bifurc.pdf",width=4,height=3) 
# continue(f,x="Kz",xmax=1e4,y="y",step=0.1,ymin=100,ymax=2e7,log="y") 
# continue(f,x="Kz",xmax=1e4,y="y",step=0.1,ymax=3e3) 
# continue(f,x="pz",xmax=1e4,y="Z",step=0.1,ymax=3e3) 
# continue(f,x="Kz",xmax=1e4,y="Z",step=0.1,ymax=3e3) 


p=c(pz=3.84e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 200, m = 1e-04)
plane(xmin=1,ymin=1,xmax=1e6,ymax=1e5,log="xy",main="pz = 3840   rz = 200 vs 400")
p=c(pz=3.84e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 400, m = 1e-04)
plane(xmin=1,ymin=1,xmax=1e6,ymax=1e5,log="xy",add=T)
dev.copy2pdf(file="~/Results/twoCities/phaseDeBoer3p8rzUp.pdf",width=4,height=3) 


p=c(pz=3.84e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 200, m = 1e-04)
plane(xmin=1,ymin=1,xmax=1e6,ymax=1e5,log="xy",main="pz = 3840 vs 4000   rz = 200")
p=c(pz=4.0e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 200, m = 1e-04)
plane(xmin=1,ymin=1,xmax=1e6,ymax=1e5,log="xy",add=T)
dev.copy2pdf(file="~/Results/twoCities/phaseDeBoer3p8to4.pdf",width=4,height=3) 

p=c(pz=3.84e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 200, m = 1e-04)
plane(xmin=1,ymin=1,xmax=1e6,ymax=1e5,log="xy",main="pz = 3840 vs 4000 vs 5000   rz = 200")
p=c(pz=4e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 200, m = 1e-04)
plane(xmin=1,ymin=1,xmax=1e6,ymax=1e5,log="xy",add=T)
p=c(pz=5e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 200, m = 1e-04)
plane(xmin=1,ymin=1,xmax=1e6,ymax=1e5,log="xy",add=T)
dev.copy2pdf(file="~/Results/twoCities/phaseDeBoer3p8to4to5.pdf",width=4,height=3)

p=c(pz=4.0e3,Kz=1e3, py = 0.25, Ky = 1e+06, a = 2, rz = 200, m = 1e-04)
plane(xmin=1,ymin=1,xmax=1e6,ymax=1e5,log="xy",main="pz = 4000   rz = 200")
mid <- newton(c(y=615.8385,z=2498.4604),plot=T) 
mid2 <- newton(c(y=1615.8385,z=2498.4604),plot=T) 
high <- newton(c(y=1e6,z=2498.4604),plot=T) 
dev.copy2pdf(file="~/Results/twoCities/phaseDeBoer4p0stab.pdf",width=4,height=3) 



