\name{raue10}
\alias{raue10}
\title{Epo model of Raue et al. 2010}
\description{This function returns the right hand side of an 
             ordinary differential equation 
             model of Epo-EpoR internalization and degradation as analyzed by 
             by Raue et al in Chaos in 2010.  
             The intended use of this function is as a quoted argument to \code{ode()} of the \code{deSolve} package. 
     }
\usage{raue10(Time, State, Pars)}

\arguments{
  \item{Time}{The parameters of this model do not depend on time: this argument exists
              here as a dummy only because deSolve expects it to exist.  }
  \item{State}{Vector of current states. The elements of this vector must be named because within this 
              function's definition is \code{with(as.list(State,Pars), {XXX} )} 
              and code in XXX presupposes that it can call named elements as variables with those names. }
  \item{Pars}{Vector of parameter values with \emph{named} (see above) elements.  }
}

\details{ A graph of the implemented model is shown below.  
This model exists as both R and C code. The C code is included in the doc directory as documentation.  
\if{html}{\figure{raue.svg}}
\if{latex}{\figure{raue.png}{options: width=5in}}
}

\value{The right hand side of the ODEs.  }

\note{  Coding of this model in C for R was supported by the National Cancer Institute and Tufts Integrative 
        Cancer Biology Program under U54CA149233-029689. }
\references{
Raue A, Becker V, Klingmüller U and Timmer J, 
Identifiability and Observability Analysis for Experimental Design in Non-Linear Dynamical Models.
 \emph{Chaos}, \bold{20}(4), 045105 (2010).
}
        
\author{Tom Radivoyevitch (\email{txr24@case.edu}) }
\seealso{\code{\link{myelo-package}, \link{rcj12}, \link{fokas91} }}
\examples{
\dontrun{
library(myelo)
dyn.load(file.path(system.file(paste("libs",Sys.getenv("R_ARCH"),sep=""), package = "myelo"),"myelo.dll"))
#test DLL
log10(c(2100,164,516))
tpts=unique(epo$time)
thetahat=c(Bmax = 2.809384, Epo0 = 3.338145,Kd= 2.583673, kde  = -1.871801,
	kdi=-2.738496, ke= -1.195154, kex=-2.557234,kon=-4.066954,kt=-1.701394,scale= 2.191042)    
p0=10^thetahat
np=as.numeric(p0) # strip names
yout <- ode(y = c(Epo=np[2],EpoR=np[1],EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), times = tpts, 
		func = "derivsc", parms = p0, 
		dllname = "myelo",initfunc = "parmsc",
		nout = 3, outnames = c("y1", "y2","y3")) 
yout  # checks out OK, so code in epo.c must be OK
}
}
\keyword{IO}
