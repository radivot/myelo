\name{rcj12}
\alias{rcj12}
\title{Granulopoiesis model of Robert C. Jackson (2012)}
\description{This function returns a list of length 2, as expected for uses with deSolve. The first list element
             is the right hand side of the ODEs.}
\usage{rcj12(Time, State, Pars)}

\arguments{
  \item{Time}{If the equations were time varying they would depend on this input. They are not. }
  \item{State}{Vector of named current states. In rcj12, with(as.list()) makes these available by name. }
  \item{Pars}{Vector of named parameter values. In rcj12, with(as.list()) makes these available by name. }
}

\details{This model is a work in progress. The current graph is depicted below.  
\if{html}{\figure{rcj.svg}}
\if{latex}{\figure{rcj.png}{options: width=5in}}
}

\value{A list of length 2 as expected by deSolve. The first list element is 
       the vector of derivatives comprised of the right hand side of the ODEs.
       The second element is a vector of auxiliary variables that we also wish to track over time. }

\note{ The intent is not to use this function directly, but rather to pass it as an argument to 
        \code{ode()} of the \code{deSolve} package.  }
        
\author{Robert C. Jackson (\email{bob@aol.com}) }
\seealso{\code{\link{myelo-package}, \link{fokas91}} }
\examples{
\dontrun{
library(myelo)
library(deSolve)
# state variables are capitalized
# P   = progenitors (old CFU-GM
# Q   = quiescent progenitors
# N   = Neutrophils
#TN   = Tissue Neutrophils
#AN   = Activated Neutrophils
# D   = dead cell
#GM   = GM-CSF   
#NAP2 = neutrophil activating peptide-2 
#RAC2 = G protein associated with increased ROS production
#Mcl1 = mantle cell lymphoma 1 (antiapoptosis factor)
#ROS  = reactive oxygen species

# Drugs used to manipulate the system
#mo = MOR103 = Ab blocker of GM   (2)
#sb = SB272844 = IL8 inhibitor    (4)
#cn = CNDAC = sapacitabine becomes this = SSB producer = DSBs in S via HR  (5)

#bu = Busulfan = DNA crosslinker  (1)  drop state variable P by 10^(2*bu/20)
#das= dasatinib                   (3)

pars=c(Vx.p=8.5e3, Vp.p=.0203, Kgm=2, Kcn=100, Vp.n=.026, Kg=1, Kmo=1,
	  Vn.tn=.781, Kil8=396, Ksb=250, Vtn.x=.038, Kmcl1=2.438, Vtn.an=.001, Knap2=.1, 
	  Van.tn=.09, Vq.p=.53778, Kp=7e4, Vx.gm=2.007, Kn=5e5, Vgm.x=.05,
      Vx.nap2=.002, Kabl=1, Kdas=5, Vnap2.x=.2, Vx.rac2=1.5, Kan=100, Vrac2.x=.1, 
	  Vn.x=2e-3, Vx.mcl1=1, Kabl2=1, Vmcl1.x=.5, Vx.ros=980, Krac2=250,
      Vros.x=1680, Kros=30, Kgsh=1e5, Vn.m=1.309e-3,
	  abl=0, g=1, il8=4,  gsh=1e6,  # boundary  conditions
	  # tnfa=0,   #bc not used 
	  Qf=.0373, volP=2900, volM=1400, volT=65700,
	  bu=0,mo=0,das=0,sb=0,cn=0) 

cellsIC=c(P=6.471e6, Q=5.663e5, N=4.561e6, TN=7.57e4, AN=76.9)
moleculesIC=c(GM=3.97,NAP2=.01004, RAC2=26.59,Mcl1=2.021, ROS=3.916)
X0=c(cellsIC,moleculesIC)
times <- seq(1, 3000, by = 1)
out   <- ode(X0, times, rcj12, pars)
head(out)
tail(out)  # Simulation 1 checks out OK
X0=out[3000,2:11]

# Simulation 2
bu=20 # set to MTD
(eventdat <- data.frame(var = "P", time = 0, value = 1/10^(2*bu/20), method = "mult"))
# Bob, this is a one time 100 fold drop in progenitor numbers due to busulfan
out   <- ode(X0, -20:600, rcj12, pars,events = list(data = eventdat))
tail(out) 

matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "Conc",
		main = "myelo_N", log="y", lwd = 2)
legend(300,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")
}
}
\keyword{IO}
