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
Raue A, Becker V, Klingmuller U and Timmer J, 
Identifiability and Observability Analysis for Experimental Design in Non-Linear Dynamical Models.
 \emph{Chaos}, \bold{20}(4), 045105 (2010).
}
        
\author{Tom Radivoyevitch (\email{txr24@case.edu}) }
\seealso{\code{\link{myelo-package}, \link{rcj12}, \link{fokas91} }}
\examples{
\dontrun{
library(myelo)
thetahat=c(Bmax = 2.821, Epo0 = 3.322,Kd= 2.583, kde  = -1.884,
		kdi=-2.730, ke= -1.177, kex=-2.447,kon=-4.091,kt=-1.758,scale= 2.210)  # from Table 1 of Chaos 2010  
(p0=10^thetahat) # Log10s used above to keep p0 positive.
yout <- ode(y = c(Epo=2099,EpoR=662,EpoEpoR=0,EpoEpoRi=0,dEpoi=0,dEpoe=0), 
        times = c(0, 5, 20, 60, 120, 180, 240, 300), 
		func = raue10, parms = p0)
yout  
plot(yout)

par(mfrow=c(1,1))
with(epo,matplot(time,epo[,-1]))
matlines(yout[,1],yout[,-(1:7)])  # Table 1 fits look good

raue10  # see definition of system of ODEs
}
}
\keyword{IO}
