\name{fokas91}
\alias{fokas91}
\title{Granulopoiesis model of Fokas et al. 1991}
\description{This function returns a list of length 2 expected by deSolve 
             for the right hand side of the ODEs of Fokas et al 1991.}
\usage{fokas91(Time, State, Pars)}

\arguments{
  \item{Time}{If the equations were time varying they would depend on this input. They are not. }
  \item{State}{Vector of named current states. \code{with(as.list())} makes these available by name. }
  \item{Pars}{Vector of named parameter values. \code{with(as.list())} makes these available by name. }
}
\details{The paper offers normal and CML parameter values. Only the case of 7 actively dividing
 states and 7 parallel non-dividing maturation states was implemented. The extended model version that 
 includes stem cell dynamics was not implemented, nor was the model that considers 11 stages of each type. 
 The 7-7 model that was implemented is depicted below.  
\if{html}{\figure{fokas.svg}}
\if{latex}{\figure{fokas.png}{options: width=5in}}
}

\value{A list of length 2 as expected by deSolve. The first list element is 
       the vector of derivatives comprised of the right hand side of the ODEs of Fokas et al 1991.
       The second element is a vector of auxiliary variables that we also wish to track over time. }


\note{ The intent is not to use this function directly, but rather to pass it as an argument to 
        \code{ode()} of the \code{deSolve} package. 
        This work was supported by the National Cancer Institute and Tufts Integrative 
        Cancer Biology Program under U54CA149233-029689. }
\references{
A.S. Fokas, J.B. Keller and B.D. Clarkson, Mathematical Model of Granulocytopoiesis and
Chronic Myelogeneous Leukemia, Cancer Research, 51, 2084-2091 (1991)}
        
\author{Tom Radivoyevitch (\email{txr24@case.edu}) }
\seealso{\code{\link{myelo-package}, \link{rcj12}} }
\examples{
\dontrun{
library(myelo)
library(deSolve)
pars=c(
fb=0.8, # blast fraction remaining active per division (nb = 3 is assumed)
fp=0.6, # progenitor fraction remaining active per division (np=2 is assumed)
fm=0.5, # myelocyte fraction remaining active per division  (nm=2 is assumed)
T= 20, # hours per stage
q=72e6/16  # flux into the blast pool per hour per kg of marrow
)  # 

X0=rep(0,14)
names(X0)<-c(paste("A",1:7,sep=""),paste("M",1:7,sep=""))
times <- seq(1, 3000, by = 1)
out   <- ode(X0, times, fokas91, pars)
tail(out)  # Mismatch with Table 3: Nm/Np = 1.9 < 2.3, Nm = 2.1e9 < 2.6e9,=one problem, c(2.3/1.9,2.6/2.1)

X0=out[3000,2:15] # start at steady state, then hit with chemo dropping all cells to 1%
(eventdat <- data.frame(var = names(X0), time = 0, value = 0.01, method = "mult"))
out   <- ode(X0, -20:300, fokas91, pars,events = list(data = eventdat))
tail(out) 

matplot(out[ , 1], out[ , c("Nb","Np","Nm")], type = "l", xlab = "time", ylab = "Marrow Cells per kg of body mass",
		main = "Model of Fokas et al. 1991", log="y", lwd = 2)
legend(150,1e8, c("Blasts","Progenitors", "Myelocytes"), col = 1:3, lty = 1:3,bty="n",lwd=2)

}
}
\keyword{IO}
