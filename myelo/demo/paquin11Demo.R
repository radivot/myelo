# Paquin D, Kim PS, Lee PP, Levy D (2011) Strategic Treatment Interruptions During 
# Imatinib Treatment of Chronic Myelogenous Leukemia Bull Math Biol (2011) 73: 1082-1100
library(PBSddesolve) 	# PBSddesolve includes the delay-differential equation (DDE) solver
#Use 'install.packages("PBSddesolve")' to install the package into R

lam = 0.75  #fraction of leukemia death rates resulting from non-immune causes

#Vector of parameter defaults, values that do not change from patient-to-patient
pdef = c(d0 = 0.003*lam, 		#stem cell (SC) death rate
         d1 = 0.008*lam, 		#progenitor cell (PC) death rate
         d2 = 0.05*lam,			#differentiated cell (DC) death rate
         d3 = lam,				#terminal cell (TC) death rate
         ry = 0.008,				#growth rate for nonresistant stem cells
         ay = 1.6,				#PC growth rate without imatinib treatment
         by = 10,				#DC growth rate without imatinib treatment
         cy = 100,				#TC growth rate without imatinib treatment
         k = 1,					#kinetic coefficient
         p0 = 0.8,			  #probability T-cell engages cancer cell. ####Set to 0 to remove immune-response  
         qc = 0.75,				#probability cancer cell dies from encounter
         qt = 0.5,				#probability T-cell survives encounter
         tau = 1) 				#duration of one T-cell division

#parameters that change from patient-to-patient (average here)
pvar1=c(n = 1.2,				#average number of T-cell divisions
        dt = 0.001,			#anti-leukemia T-cell death rate
        st = 1.2e-6,			#anti-leukemia T-cell supply rate
        cn = 1,					#decay rate of immune responsivity
        y0 = 1.8e-5) 		#initial concentration of leukemia stem cells
pvarall = c(pdef,pvar1)

IC = function(parms){		#Function that produces the set of initial conditions, given the model parameters
  #Returns a vector of the given parameters concatenated with the initial conditions of 
  #nonresistant SC, PC, DC, TC cells (y0,y1,y2,y3),
  #and T-cells (Y).
  initialconditions = c(parms["y0"])
  initialconditions = c(initialconditions,parms["ay"]*initialconditions[1]/parms["d1"])
  initialconditions = c(initialconditions,parms["by"]*initialconditions[2]/parms["d2"])
  initialconditions = c(initialconditions,parms["cy"]*initialconditions[3]/parms["d3"])
  initialconditions = c(initialconditions,parms["st"]/parms["dt"])
  names(initialconditions) = c("y0","y1","y2","y3","Y")
  return(c(parms,initialconditions))
}

pvarall=IC(pvarall)

paquin11 = function(t,y,p){		#set of DDEs
  pf = function(C,T,parms){		#Function for the rate at which leukemic cells are destroyed by T-cells,
    #given the number of C (leukemic) and T-cells, and set of parameters.
    return(parms["p0"]*exp(-parms["cn"]*C)*parms["k"]*T)
  }
  C = sum(y[-5])
  const = p["qc"]*pf(C,y[5],p)
  dy0 = (p["ry"] - p["d0"])*y[1]-const*y[1]
  dy1 = (p["ay"]/100)*y[1]-p["d1"]*y[2]-const*y[2]
  dy2 = (p["by"]/750)*y[2]-p["d2"]*y[3]-const*y[3]
  dy3 = p["cy"]*y[3]-p["d3"]*y[4]-const*y[4]
  
  if(t < p["n"]*p["tau"]) lag = p[19:23] #the initial conditions
  else lag = pastvalue(t - p["n"]*p["tau"]) #change this to remove the time delay, done in yprimenolag
  Cnt = sum(lag[-5])
  dT = p["st"] - p["dt"]*y[5] - pf(C,y[5],p)*C + 2^(p["n"])*pf(Cnt,lag[5],p)*p["qt"]*Cnt
  return(c(dy0,dy1,dy2,dy3,dT))
}

timerange = 0:1200 #the range of time to simulate the model (in days)
results = dde(pvarall[19:23],timerange,paquin11,pvarall)
data = cbind(results[,1],apply(results[,2:5],1,sum),results[,6])

#Plot Fig. 2
par(mfrow = c(1,2))
plot(data[,1],250*data[,2],type = "l",col = "red",log="y",
     main = "Leukemic Cells",xlab = "Time (days)",ylab = "Concentration")
plot(data[,1],data[,3],type = "l",col = "blue",
     main = "T-Cells",xlab = "Time (days)",ylab = "Concentration")



paquin11Fig3 = function(t,y,p){		#set of DDEs
  pf = function(C,T,parms){		#Function for the rate at which leukemic cells are destroyed by T-cells,
    #given the number of C (leukemic) and T-cells, and set of parameters.
    return(parms["p0"]*exp(-parms["cn"]*C)*parms["k"]*T)
  }
  C = sum(y[-5])
  const = p["qc"]*pf(C,y[5],p)
  dy0 = (p["ry"] - p["d0"])*y[1]-const*y[1]
  if (t>300&t<315) {
    dy1 = (p["ay"])*y[1]-p["d1"]*y[2]-const*y[2]
    dy2 = (p["by"])*y[2]-p["d2"]*y[3]-const*y[3]
  } else {
    dy1 = (p["ay"]/100)*y[1]-p["d1"]*y[2]-const*y[2]
    dy2 = (p["by"]/750)*y[2]-p["d2"]*y[3]-const*y[3]
  }
  dy3 = p["cy"]*y[3]-p["d3"]*y[4]-const*y[4]
  
  if(t < p["n"]*p["tau"]) lag = p[19:23] #the initial conditions
  else lag = pastvalue(t - p["n"]*p["tau"]) #change this to remove the time delay, done in yprimenolag
  Cnt = sum(lag[-5])
  dT = p["st"] - p["dt"]*y[5] - pf(C,y[5],p)*C + 2^(p["n"])*pf(Cnt,lag[5],p)*p["qt"]*Cnt
  return(c(dy0,dy1,dy2,dy3,dT))
}

results = dde(pvarall[19:23],timerange,paquin11Fig3,pvarall)
data = cbind(results[,1],apply(results[,2:5],1,sum),results[,6])

#Plot Fig. 3
par(mfrow = c(1,2))
plot(data[,1],250*data[,2],type = "l",col = "red",log="y",
     main = "Leukemic Cells",xlab = "Time (days)",ylab = "Concentration")
plot(data[,1],data[,3],type = "l",col = "blue",
     main = "T-Cells",xlab = "Time (days)",ylab = "Concentration")


###### redo Fig.3 using deSolve
library(deSolve)

paq11 = function(t,y,p){		#set of DDEs
  pf = function(C,T,parms){		#Function for the rate at which leukemic cells are destroyed by T-cells,
    #given the number of C (leukemic) and T-cells, and set of parameters.
    return(parms["p0"]*exp(-parms["cn"]*C)*parms["k"]*T)
  }
  C = sum(y[-5])
  const = p["qc"]*pf(C,y[5],p)
  dy0 = (p["ry"] - p["d0"])*y[1]-const*y[1]
  if (t>300&t<315) {
    dy1 = (p["ay"])*y[1]-p["d1"]*y[2]-const*y[2]
    dy2 = (p["by"])*y[2]-p["d2"]*y[3]-const*y[3]
  } else {
    dy1 = (p["ay"]/100)*y[1]-p["d1"]*y[2]-const*y[2]
    dy2 = (p["by"]/750)*y[2]-p["d2"]*y[3]-const*y[3]
  }
  dy3 = p["cy"]*y[3]-p["d3"]*y[4]-const*y[4]
  
  if(t < p["n"]*p["tau"]) lag = p[19:23] #the initial conditions
  else lag = lagvalue(t - p["n"]*p["tau"]) 
  Cnt = sum(lag[-5])
  dT = p["st"] - p["dt"]*y[5] - pf(C,y[5],p)*C + 2^(p["n"])*pf(Cnt,lag[5],p)*p["qt"]*Cnt
  return(list(c(dy0,dy1,dy2,dy3,dT),C=sum(y[1:4])))
}

timerange = 0:1200 #the range of time to simulate the model (in days)
results = dede(pvarall[19:23],timerange,paq11,pvarall)
data = cbind(results[,1],results[,7],results[,6])

#Plot Fig. 3
par(mfrow = c(1,2))
plot(data[,1],250*data[,2],type = "l",col = "red",log="y",
     main = "Leukemic Cells",xlab = "Time (days)",ylab = "Concentration")
plot(data[,1],data[,3],type = "l",col = "blue",
     main = "T-Cells",xlab = "Time (days)",ylab = "Concentration")


