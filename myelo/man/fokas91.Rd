\name{fokas91}
\alias{fokas91}
\title{Granulopoiesis model of Fokas et al. 1991}
\usage{
  fokas91(Time, State, Pars)
}
\arguments{
  \item{Time}{The parameters of this model do not depend on
  time: this argument exists here as a dummy only because
  deSolve expects it to exist.}

  \item{State}{Vector of current states. The elements of
  this vector must be named because within this function's
  definition is \code{with(as.list(State,Pars), {XXX} )}
  and code in XXX presupposes that it can call named
  elements as variables with those names.}

  \item{Pars}{Vector of parameter values with \emph{named}
  (see above) elements.}
}
\value{
  A list of length 2 (as expected by deSolve) where the
  first list element is the vector of derivatives, i.e. the
  right hand side of the ODEs of Fokas et al 1991, and the
  second element of the list is a vector of auxiliary
  variables that one may wish to track over time.
}
\description{
  This function returns the right hand side of an ordinary
  differential equation model of normal and CML
  granulopoiesis that was published by Fokas, Keller and
  Clarkson in Cancer Research in 1991.  The intended use of
  this function is as an argument to \code{ode()} of the
  \code{deSolve} package.
}
\details{
  The extended version of this model that includes stem
  cell dynamics, and the extension that considers 11 stages
  of granulopoiesis, were not implemented. Implemented here
  then is only the case of 7 actively dividing
  states/stages and 7 parallel non-dividing maturation
  states/stages.  This model is depicted below.
  \if{html}{\figure{fokas.svg}}
  \if{latex}{\figure{fokas.png}{options: width=5in}}
}
\note{
  This work was supported by the National Cancer Institute
  and Tufts Integrative Cancer Biology Program under
  U54CA149233-029689.
}
\examples{
\dontrun{
library(myelo)
library(deSolve)
pars=c(
fb=0.8,
fp=0.6,
fm=0.5,
T= 20,
q=72e6/16
)  #
# where
# fb=blast fraction remaining active per division (nb = 3 is assumed)
# fp=progenitor fraction remaining active per division (np=2 is assumed)
# fm=myelocyte fraction remaining active per division  (nm=2 is assumed)
# T=hours per stage
# q=flux into the blast pool per hour per kg of marrow

X0=rep(0,14)
names(X0)<-c(paste("A",1:7,sep=""),paste("M",1:7,sep=""))
times <- seq(1, 3000, by = 1)
out   <- ode(X0, times, fokas91, pars)
tail(out)  # Mismatch with Table 3: Nm/Np = 1.9 < 2.3, Nm = 2.1e9 < 2.6e9,=one problem, c(2.3/1.9,2.6/2.1)

X0=out[3000,2:15] # start at steady state, then hit with chemo dropping all cells to 1\%
(eventdat <- data.frame(var = names(X0), time = 0, value = 0.01, method = "mult"))
out   <- ode(X0, -20:300, fokas91, pars,events = list(data = eventdat))
tail(out)

matplot(out[ , 1], out[ , c("Nb","Np","Nm")], type = "l", xlab = "time", ylab = "Marrow Cells per kg of body mass",
		main = "Model of Fokas et al. 1991", log="y", lwd = 2)
legend(150,1e8, c("Blasts","Progenitors", "Myelocytes"), col = 1:3, lty = 1:3,bty="n",lwd=2)

}
}
\author{
  Tom Radivoyevitch (\email{txr24@case.edu})
}
\references{
  A.S. Fokas, J.B. Keller and B.D. Clarkson, Mathematical
  Model of Granulocytopoiesis and Chronic Myelogeneous
  Leukemia, \emph{Cancer Research}, \bold{51}, 2084-2091
  (1991)
}
\seealso{
  \code{\link{myelo-package}, \link{rcj12}}
}
\keyword{IO}

