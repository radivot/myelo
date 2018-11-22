# Kim PS, Lee PP, Levy D (2008) Dynamics and Potential Impact of the Immune
# Response to Chronic Myelogenous Leukemia. PLoS Comput Biol 4(6): e1000095.
# doi:10.1371/journal.pcbi.1000095
library(PBSddesolve) 	# PBSddesolve includes the delay-differential equation (DDE) solver
						#Use 'install.packages("PBSddesolve")' to install the package into R

lam = 0.75  #fraction of leukemia death rates resulting from non-immune causes
            #changing this to 1 may be necessary when ignoring the immune-response

#Vector of parameter defaults, values that do not change from patient-to-patient
pdef = c(d0 = 0.003*lam, 		#stem cell (SC) death rate
	d1 = 0.008*lam, 		#progenitor cell (PC) death rate
	d2 = 0.05*lam,			#differentiated cell (DC) death rate
	d3 = lam,				#terminal cell (TC) death rate
	ry = 0.008,				#growth rate for nonresistant stem cells
	ay = 1.6,				#PC growth rate without imatinib treatment
	by = 10,				#DC growth rate without imatinib treatment
	cy = 100,				#TC growth rate without imatinib treatment
	ayp = 0.016,			#PC growth rate with imatinib treatment
	byp = 1/75,				#DC growth rate with imatinib treatment
	cyp = 100,				#TC growth rate with imatinib treatment
	rz = 0.023,				#growth rate for resistant cells
	az = 1.6,				#PC growth rate for resistant cells
	bz = 10,				#DC growth rate for resistant cells
	cz = 100,				#TC growth rate for resistant cells
	u = 0,					#mutation rate per division (nonresistant -> resistant)
	k = 1,					#kinetic coefficient
	p0 = 0.8,			  #probability T-cell engages cancer cell. ####Set to 0 to remove immune-response  
	qc = 0.75,				#probability cancer cell dies from encounter
	qt = 0.5,				#probability T-cell survives encounter
	tau = 1) 				#duration of one T-cell division
	
	
#parameters that change from patient-to-patient for patient P1
pvar1=c(n = 1.2,				#average number of T-cell divisions
	     dt = 0.001,			#anti-leukemia T-cell death rate
	     st = 1.2e-6,			#anti-leukemia T-cell supply rate
	     cn = 1,					#decay rate of immune responsivity
	     y0 = 7.6e-6) 		#initial concentration of leukemia stem cells
pvar4  = c(n = 2.2,   dt = 0.0022, st = 9e-7,    cn = 7,   y0 = 2.4e-6) #patient P4
pvar12 = c(n = 1.17,  dt = 0.007,  st = 3.08e-5, cn = 0.8, y0 = 1.2e-5) #patient P12

pvarall = list(P1 = c(pdef,pvar1),P4 = c(pdef,pvar4),P12 = c(pdef,pvar12))	
	#Puts the vectors of variable parameters for P1, P4, P12 into a list, and
	#concatenates in front of each vector the set of fixed parameters

IC = function(parms){		#Function that produces the set of initial conditions, given the model parameters
							#Returns a vector of the given parameters concatenated with the initial conditions of 
							#nonresistant SC, PC, DC, TC cells (y0,y1,y2,y3),
							#resistant SC, PC, DC, TC cells (z0,z1,z2,z3),
							#and T-cells (Y).
	initialconditions = c(parms["y0"])
	initialconditions = c(initialconditions,parms["ay"]*initialconditions[1]/parms["d1"])
	initialconditions = c(initialconditions,parms["by"]*initialconditions[2]/parms["d2"])
	initialconditions = c(initialconditions,parms["cy"]*initialconditions[3]/parms["d3"],0) #OR 1e-9
	initialconditions = c(initialconditions,parms["az"]*initialconditions[5]/parms["d1"])
	initialconditions = c(initialconditions,parms["bz"]*initialconditions[6]/parms["d2"])
	initialconditions = c(initialconditions,parms["cz"]*initialconditions[7]/parms["d3"],parms["st"]/parms["dt"])
	names(initialconditions) = c("y0","y1","y2","y3","z0","z1","z2","z3","Y")
	return(c(parms,initialconditions))
}

pvarall = lapply(pvarall,IC)	#Adds the initial conditions to each set of parameters


# kim08 = function(t,y,p){		#set of DDEs
#   pf = function(C,T,parms){		#Function for the rate at which leukemic cells are destroyed by T-cells,
#     #given the number of C (leukemic) and T-cells, and set of parameters.
#     return(parms["p0"]*exp(-parms["cn"]*C)*parms["k"]*T)
#   }
#   C = sum(y[-9])
#   const = p["qc"]*pf(C,y[9],p)
#   dy0 = (p["ry"]*(1-p["u"]) - p["d0"])*y[1]-const*y[1]
#   dy1 = p["ayp"]*y[1]-p["d1"]*y[2]-const*y[2]
#   dy2 = p["byp"]*y[2]-p["d2"]*y[3]-const*y[3]
#   dy3 = p["cyp"]*y[3]-p["d3"]*y[4]-const*y[4]
#   dz0 = (p["rz"] - p["d0"])*y[5]+p["ry"]*y[1]*p["u"] - const*y[5]
#   dz1 = p["az"]*y[5] - p["d1"]*y[6] - const*y[6]
#   dz2 = p["bz"]*y[6] - p["d2"]*y[7] - const*y[7]
#   dz3 = p["cz"]*y[7] - p["d3"]*y[8] - const*y[8]
#   
#   if(t < p["n"]*p["tau"]) lag = p[27:35] #the initial conditions
#   else lag = pastvalue(t - p["n"]*p["tau"]) #change this to remove the time delay, done in yprimenolag
#   Cnt = sum(lag[-9])
#   dT = p["st"] - p["dt"]*y[9] - pf(C,y[9],p)*C + 2^(p["n"])*pf(Cnt,lag[9],p)*p["qt"]*Cnt
#  return(c(dy0,dy1,dy2,dy3,dz0,dz1,dz2,dz3,dT))
# }

timerange = 0:1500 #the range of time to simulate the model (in days)
results = lapply(pvarall,function(x) dde(x[27:35],timerange,kim08,x))
      #uses the delay differential equation solver to simulate the model for each given
      #set of parameters in the list "pvarall"
      #returns a list, where each element in the list is a matrix that gives the 
      #resulting simulations with the corresponding parameters
resultsSum = lapply(results,function(x) cbind(x[,1],apply(x[,2:9],1,sum),x[,10]))
      #finds the sum of all leukemic cells (nonresistant and resistant, SC, PC, etc.)

#Plotting:
par(mfrow = c(ceiling(length(resultsSum)/2),2))
for(i in 1:length(resultsSum)){
	data = resultsSum[[i]]
	plot(data[,1],data[,2],ylim = c(0,0.06),type = "l",col = "red",
	     main = paste("Simulation for",names(resultsSum)[i]),xlab = "Time (days)",ylab = "Concentration")
	lines(data[,1],data[,3],col = "blue")
	legend("topright",legend = c("Leukemic Cells/2500","T-Cells"),col = c("red","blue"),lwd = c(1,1),cex=0.6,bty="n")
}
