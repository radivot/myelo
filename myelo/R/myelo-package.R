#'Erythropoietin (Epo) Dynamic Response Data of Becker et al. Science 2010
#'
#'This is the portion of the Becker et al. Science 2010 data that is used in
#'the profile likelihood analysis Raue et al. Chaos 2010.
#'
#'The values involve a scale factor.
#'
#'@name epo
#'@docType data
#'@format A data frame that includes time and 3 variables.  \describe{
#'\item{time}{The time point in minutes after Epo addition.}
#'\item{epo.e}{The counts per minute (cpm) of Epo external to the cell,
#'degraded or intact, but free, see \link{raue10}.} 
#'\item{epo.m}{The
#'cpm of Epo bound to the Epo Receptor (EpoR) on the cell's plasma membrane.}
#'\item{epo.i}{The cpm of Epo internal to the cell, either on an
#'internal EpoR or degraded and free.}  }
#'@note Coding supported by the National Cancer Institute and Tufts
#'Integrative Cancer Biology Program U54CA149233-029689.
#'@seealso \code{\link{myelo-package}, \link{rcj12}, \link{fokas91},
#'\link{raue10} }
#'@references Becker, V., Schilling, M., Bachmann, J., Baumann, U., Raue, A.,
#'Maiwald, T., Timmer, J., and Klingmuller, U., Covering a broad dynamic range:
#'Information processing at the erythropoietin receptor, \emph{Science}
#'\bold{328}, 1404-1408, 2010.
#'@source The data were taken from the excel file in the zip download of Epo
#'receptor model under \bold{Applications} at
#'\url{http://www.fdmold.uni-freiburg.de/~araue/pmwiki.php/Main/Software}
#'@keywords datasets
#'@examples
#'
#'library(myelo)
#'epo
#'with(epo,matplot(time,epo[,-1]))
#'
NULL


#'G-CSF Dynamic Response Data of van der Auwera et al. 2001
#'
#'In addition to the single injection data used in Scholz 2012 this data.frame
#'includes CFU-GM and CD34 cells counts in response to 7 daily filgrastim injections.
#'
#'@name auwera
#'@docType data
#'@format A data frame that includes time in hours and 8 variables.  \describe{
#'\item{hours}{Time points in hours after subcutanteous filgrastim injection.}
#'\item{gcsf5}{GCSF time course in ng/mL after a 5 ug/kg injection.} 
#'\item{gcsf10}{GCSF time course in ng/mL after a 10 ug/kg injection.} 
#'\item{cd34p5}{The number of CD34+ cells per cubic mm after 7 daily 5 ug/kg injection.}
#'\item{cd34p10}{The number of CD34+ cells per cubic mm after 7 daily 10 ug/kg injection.}
#'\item{cfugm5}{The number of CFUGM per 10^5 cells after 7 daily 5 ug/kg injections.}
#'\item{cfugm10}{The number of CFUGM per 10^5 cells after 7 daily 10 ug/kg injections.}
#'\item{anc5}{Absolute Neutrophil Count (ANC) per cubic mm time course after a 5 ug/kg injection.} 
#'\item{anc10}{ANC per cubic mm time course after a 10 ug/kg injection.} }
#'
#'@note Coding supported by the National Cancer Institute and Tufts
#'Integrative Cancer Biology Program U54CA149233-029689.
#'@seealso \code{\link{myelo-package}, \link{rcj12}, \link{fokas91},
#'\link{raue10} }
#'@references Phillippe van der Auwera, Erich Platzer, Z.-X. Xu, Rainer Schulz, Olivier Feugeas,
#' Renaud Capdeville, and David J. Edwards, Pharmacodynamics and Pharmacokinetics of Single
#'Doses of Subcutaneous Pegylated Human G-CSF Mutant
#'(Ro 25-8315) in Healthy Volunteers: Comparison With
#'Single and Multiple Daily Doses of Filgrastim, \emph{American Journal of Hematology}
#'\bold{66}, 245?251, 2001.  
#'@source The data were obtained from Figures 1,2,3, and 5 of vander Auwera et al. by first copying
#' into powerpoint, then saving as a gif, and then using PlotDigitizer. 
#'@keywords datasets
#'@examples
#'
#'library(myelo)
#'auwera
#'par(mfrow=c(2,2),oma=c(0,0,2,0))
#'with(auwera,{matplot(hours,auwera[,2:3],ylab="GCSF")
#'matplot(hours,auwera[,4:5],ylab="CD34")
#'matplot(hours,auwera[,6:7],ylab="CFU-GM")
#'matplot(hours,auwera[,8:9],ylab="ANC")}
#')
#'title("1 = 5 ug/kg           2 = 10 ug/kg",outer=TRUE)
#'
NULL





#'Model parameter values of Scholz et al. TBMM 2012
#'
#'This is a numeric vector of 128 parameter values of the model of Scholz et al. TBMM 2012. 
#'
#'@name scholzPars
#'@docType data
#'@format A vector of parameter values.  \describe{
#'See description in comments in sholzPars.R of the doc folder. }
#'@note Coding supported by the National Cancer Institute and Tufts
#'Integrative Cancer Biology Program U54CA149233-029689.
#'@seealso \code{\link{myelo-package}, \link{rcj12}, \link{fokas91},\link{scholz12},
#'\link{raue10} }
#'@references  M. Scholz, S. Schirm, M. Wetzler, C. Engel
#'and M. Loeffler, Pharmacokinetic and -dynamic modelling of
#'G-CSF derivatives in humans, \emph{Theoretical Biology and Medical Modelling} 
#' \bold{9} 32 (2012).
#'@source The parameter values were taken from the supplement PDF of the TBMM 2012 paper. 
#'@keywords datasets
#'@examples
#'
#'library(myelo)
#'scholzPars
#'
NULL


#'Model parameter values of Zhuge et al. JTB 2012
#'
#'This is a numeric vector that holds the parameter values of the model of Zhuge et al. JTB 2012. 
#'
#'@name zhugePars
#'@docType data
#'@format A vector of parameter values.  \describe{
#'See description in comments in zhugePars.R of the doc folder. }
#'@note Coding supported by the National Cancer Institute and Tufts
#'Integrative Cancer Biology Program U54CA149233-029689.
#'@seealso \code{\link{myelo-package}, \link{rcj12}, \link{fokas91},\link{scholz12},
#'\link{brooksPars},\link{raue10} }
#'@references  Changjing Zhuge, Jinzhi Lei, and Michael C. Mackey,
#'Neutrophil dynamics in response to chemotherapy and G-CSF, \emph{Journal of Theoretical Biology} 
#' \bold{293} 111?120 (2012).
#'@source The parameter values were taken from Table 1 of the PDF of the JTB 2012 paper. 
#'@keywords datasets
#'@examples
#'
#'library(myelo)
#'zhugePars
#'
NULL

#'Model parameter values of Brooks et al. JTB 2012
#'
#'This is a numeric vector that holds the parameter values of the model of Brooks et al. JTB 2012. 
#'
#'@name brooksPars
#'@docType data
#'@format A vector of parameter values.  \describe{
#'See description in comments in brooksPars.R of the doc folder. }
#'@note Coding supported by the National Cancer Institute and Tufts
#'Integrative Cancer Biology Program U54CA149233-029689.
#'@seealso \code{\link{myelo-package}, \link{rcj12}, \link{fokas91},\link{scholz12},
#'\link{zhugePars},\link{raue10} }
#'@references  Grace Brooks, Gabriel Provencher,  Jinzhi Lei, and Michael C. Mackey,
#'Neutrophil dynamics after chemotherapy and G-CSF: The role
#'of pharmacokinetics in shaping the response, \emph{Journal of Theoretical Biology} 
#'\bold{315} 97?109 (2012).
#'@source The parameter values were taken from Table 1 of the PDF of the JTB 2012 paper. 
#'@keywords datasets
#'@examples
#'
#'library(myelo)
#'brooksPars
#'
NULL



#'Models of myeloid hematopoiesis
#'
#'This package aims to make analyses of models of myeloid hematopoiesis more 
#'readily reproducible. In the case of ODE models it provides functions that
#'compute the right hand sides of systems of ODEs.
#'
#'
#'\tabular{ll}{ Package: \tab \pkg{myelo}\cr Type: \tab Package\cr Depends:
#'\tab deSolve\cr Suggests: \tab bbmle \cr License: \tab GPL-2\cr LazyLoad:
#'\tab yes\cr LazyData: \tab yes\cr URL: \tab
#'\url{http://epbi-radivot.cwru.edu/myelo/myelo.html}\cr }
#'
#'@name myelo-package
#'@docType package
#'@note This work was supported by the National Cancer Institute and Tufts
#'Integrative Cancer Biology Program under U54CA149233-029689.
#'@author Tom Radivoyevitch (\email{txr24@@case.edu})
#'@seealso \code{\link{fokas91},\link{rcj12}}
#'@keywords package
NULL

