\name{rcj12}
\alias{rcj12}
\title{Granulopoiesis model of Robert C. Jackson (2012)}
\description{This function returns the right hand side of the ODEs of a model that Dr. Robert C. Jackson is 
             actively developing. The intent is not to use this function directly, but rather to pass it as an argument to 
        \code{ode()} of the \code{deSolve} package. }
\usage{rcj12(Time, State, Pars)}

\arguments{
  \item{Time}{The parameters of this model do not depend on time: this argument exists
              here as a dummy only because deSolve expects it to exist.  }
  \item{State}{Vector of current states. The elements of this vector must be named because within this 
              function's definition is \code{with(as.list(State,Pars), {XXX} )} 
              and code in XXX presupposes that it can call named elements as variables with those names. }
  \item{Pars}{Vector of parameter values with \emph{named} (see above) elements.  }
}

\details{This model is a work in progress. Its current form is depicted below as a graph.  
\if{html}{\figure{rcj.svg}}
\if{latex}{\figure{rcj.png}{options: width=5in}}
}

\value{A list of length 2 (as expected by deSolve) where the first list element is 
       a vector of the derivatives (i.e. the right hand side of the ODEs), and
       the second element of the list is a vector of auxiliary variables that one may wish to track over time. }

\author{Robert C. Jackson (\email{bob@aol.com}) }
\seealso{\code{\link{myelo-package}, \link{fokas91}} }
\examples{
\dontrun{
library(myelo)
pars=c(Vx.p=1.7e4, Kil3=1, Vp.p=.0203, Kgm=2, Kcn=100, Vp.n=.026, Kg=1, Kmo=1,
		Vn.tn=.781, Kil8=396, Ksb=250, Vtn.x=.038, Kmcl1=2.438, Vtn.an=.001, Knap2=.1, 
		Van.tn=.09, Vq.p=.53778, Kp=7e4, Vx.gm=2.007, Kn=5e5, Vgm.x=.05,
		Vx.nap2=.002, Kabl=1, Kdas=5, Vnap2.x=.2, Vx.rac2=1.5, Kan=100, Vrac2.x=.1, 
		Vn.x=2e-3, Vx.mcl1=1, Kabl2=1, Vmcl1.x=.5, Vx.ros=980, Krac2=250,
		Vgsh.gssg=1680, Vgssg.gsh=150,Kgssg=30,Kros=30, Kgsh=1e5, Vn.m=1.309e-3,Kros2=.026,
		Kseli=.3, Kifna=1, Vbact.x=5.756e-2, nHr=2.5,kkmax=320,Vx.abl=.22,Ktnfa=1,Kimab=1, Vabl.x=.0145,Vil8.x=5e-4,
		Vtnfa.x=5e-6,  
		g=1, il3=1, gsh=1e6,                          # these are boundary  conditions
		Qf=.0373, volP=2900, volM=1400, volT=65700,
		bu=0,mo=0,das=0,sb=0,cn=0,seli=0,imab=0,ifna=0,antiox=7)      # Drug concentrations

# Starting values for state variables
X0=c(P=6.466e6, Q=5.655e5, N=4.558e6, TN=7.564e4, AN=76.84,
		GM=3.97,NAP2=.01004, RAC2=26.59, Mcl1=2.021, ROS=1.967, bact=0, abl=0, il8=4, tnfa=0, GSSG=0)  
times <- seq(1,80, by = 1)
out   <- ode(X0, times, rcj12, pars)
matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
		main = "Simulation 1: Uninhibited control", log="y", lwd = 2)
legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")

# DNA cross-linking
X0=out[80,2:16]
bu=20 # set to MTD
(eventdat <- data.frame(var = "P", time = 0, value = 1/10^(2*bu/20), method = "mult"))
out   <- ode(X0, -20:600, rcj12, pars,events = list(data = eventdat))

matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
		main = "Simulation 2: busulfan treatment", log="y", lwd = 2)
legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")
} 
}
\keyword{IO}
