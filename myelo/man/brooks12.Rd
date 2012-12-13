\name{brooks12}
\alias{brooks12}
\title{Neutrophil dynamics model of Brooks et al 2012}
\usage{
  brooks12(Time, State, Pars)
}
\arguments{
  \item{Time}{The current time.}

  \item{State}{Vector of current states. The elements of
  this vector must be named because within this function's
  definition is \code{with(as.list(State,Pars), {XXX} )}
  and code in XXX presupposes that it can call named
  elements as variables with those names.}

  \item{Pars}{Vector of parameter values with \emph{named}
  (see above) elements.}
}
\value{
  A list of length 1 where the first list element is the
  vector of derivatives, i.e. the right hand side of the
  ODEs
}
\description{
  This function returns the right hand side of an ordinary
  differential delay equation model of neutrophil dynamics
  published by Brooks et al. in 2012 The intended use of
  this function is as an argument to \code{dede()} of the
  \code{deSolve} package.
}
\details{
  In this model Q = quiescent stem cells and N =
  neutrophils
}
\note{
  This work was supported by the National Cancer Institute
  and Tufts Integrative Cancer Biology Program under
  U54CA149233-029689.
}
\examples{
\dontrun{
library(myelo)

times <- c(-20:-1,seq(-.1,1,by=0.02),2:14,seq(14.01,15,by=0.02),16:100)
yout <- dede(c(Q=brooksPars[["Qss"]],N=brooksPars[["Nss"]],EtaNP=0,Gam0=0,GamS=0,G=0,X=0,C=0),
		times = times, func = brooks12,	parms = brooksPars)
# this should make Figure 3A. The response appears to be off by a dose-amplitude parameter
par(mfrow=c(1,1),mar=c(4.7,6,3.3,1),oma=c(0,0,0,0),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
plot(times,yout[,3]/1e8,type="l",xlab="days",ylab="Neutrophils")


}
}
\author{
  Tom Radivoyevitch (\email{txr24@case.edu})
}
\references{
  Grace Brooks, Gabriel Provencher, Jinzhi Lei, Michael C.
  Mackey, Neutrophil dynamics after chemotherapy and G-CSF:
  The role of pharmacokinetics in shaping the response,
  \emph{Journal of Theoretical Biology} \bold{315} 97-109
  (2012).
}
\seealso{
  \code{\link{myelo-package}, \link{rcj12}}
}
\keyword{IO}

