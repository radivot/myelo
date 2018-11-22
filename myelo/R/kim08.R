#'CML-T-Cell Model of Kim et al 2008
#'
#'
#'This function returns the right hand side of a delay differential equation
#'model of T cells killing CML cells.  The intended use of this function is as
#'an argument to \code{dde()} of the \code{PBSdeSolve} package.
#'
#'
#'@param t The current time in the simulation.
#'@param y Vector of current states. 
#'@param p Vector of parameter values with \emph{named} (see above) elements.
#'@return A vector of derivatives (as expected by PBSddeSolve)
#'@author Tom Radivoyevitch and Rainer Sachs
#'@seealso \code{\link{myelo-package}, \link{clapp15}, \link{michor05}}
#'@references  Kim PS, Lee PP, Levy D, Dynamics and Potential Impact of the Immune
#' Response to Chronic Myelogenous Leukemia, \emph{PLoS Comput Biol}, \bold{4}, 6 (2008)
#'@keywords IO
#'@name kim08
#'@export
#'@import PBSddesolve

kim08 = function(t,y,p){		#set of DDEs
  pf = function(C,T,parms){		#Function for the rate at which leukemic cells are destroyed by T-cells,
    #given the number of C (leukemic) and T-cells, and set of parameters.
    return(parms["p0"]*exp(-parms["cn"]*C)*parms["k"]*T)
  }
  C = sum(y[-9])
  const = p["qc"]*pf(C,y[9],p)
  dy0 = (p["ry"]*(1-p["u"]) - p["d0"])*y[1]-const*y[1]
  dy1 = p["ayp"]*y[1]-p["d1"]*y[2]-const*y[2]
  dy2 = p["byp"]*y[2]-p["d2"]*y[3]-const*y[3]
  dy3 = p["cyp"]*y[3]-p["d3"]*y[4]-const*y[4]
  dz0 = (p["rz"] - p["d0"])*y[5]+p["ry"]*y[1]*p["u"] - const*y[5]
  dz1 = p["az"]*y[5] - p["d1"]*y[6] - const*y[6]
  dz2 = p["bz"]*y[6] - p["d2"]*y[7] - const*y[7]
  dz3 = p["cz"]*y[7] - p["d3"]*y[8] - const*y[8]
  
  if(t < p["n"]*p["tau"]) lag = p[27:35] #the initial conditions
  else lag = pastvalue(t - p["n"]*p["tau"]) 
  Cnt = sum(lag[-9])
  dT = p["st"] - p["dt"]*y[9] - pf(C,y[9],p)*C + 2^(p["n"])*pf(Cnt,lag[9],p)*p["qt"]*Cnt
  return(c(dy0,dy1,dy2,dy3,dz0,dz1,dz2,dz3,dT))
}
