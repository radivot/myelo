#'PKPD G-CSF model of Scholz et al. TBMM 2012, without pegfilgrastim or chemotoxicities
#'
#'This function returns the right hand side of an ordinary differential
#'equation model of G-CSF published by M. Scholz et al. in 2012 
#'(in Theoretical Biology and Medical Modelling). Subcutaneous injections are not modeled.
#'The intended use of this function is as an argument to \code{ode()} of the \code{deSolve} package.
#'
#'The model is depicted below  
#'\if{html}{\figure{scholz.png}} \if{latex}{\figure{scholz.png}{options: width=5in}}
#'
#' Here S = stem cells, CG = colony forming units of granulocytes
#' and macrophages, PGB = proliferating granulopoietic blasts, 
#' MGB = maturing granulopoietic blasts
#' The system is regulated by G-CSF.
#'
#'@param Time The parameters of this model do not depend on time: this argument
#'exists here as a dummy only because deSolve expects it to exist.
#'@param State Vector of current states. The elements of this vector must be
#'named because within this function's definition is
#'\code{with(as.list(State,Pars), {XXX} )} and code in XXX presupposes that it
#'can call named elements as variables with those names.
#'@param Pars Vector of parameter values with \emph{named} (see above)
#'elements.
#'@return A list of length 2 (as expected by deSolve) where the first list
#'element is the vector of derivatives, i.e. the right hand side of the ODEs 
#' and the second element of the list is a vector of auxiliary
#'variables that one may wish to track over time.
#'@note This work was supported by the National Cancer Institute and Tufts
#'Integrative Cancer Biology Program under U54CA149233-029689.
#'@author Tom Radivoyevitch (\email{txr24@@case.edu})
#'@seealso \code{\link{myelo-package}, \link{rcj12}}
#'@references  M. Scholz, S. Schirm, M. Wetzler, C. Engel
#'and M. Loeffler, Pharmacokinetic and -dynamic modelling of
#'G-CSF derivatives in humans, \emph{Theoretical Biology and Medical Modelling} 
#' \bold{9} 32 (2012).
#'@keywords IO
#'@export
#'@examples
#'
#'\dontrun{
#'library(myelo)
#' X0=rep(1,24)
#'names(X0)<-c("g","g1","g2","g3","g4","S","CG","PGB","G4a","G4b","G4c","G4d","G4e",
#'		"G5a","G5b","G5c","G5d","G5e","G6a","G6b","G6c","G6d","G6e","GRA")
#'# WARNING: the following is needed because the model is not completely specified! 
#'APGBout=APGBin=ACGout=ACGin=GSS=G6nor=GRAnor=1
#'out   <- ode(X0,1:300,scholz12, scholzPars)
#'rbind(head(out),tail(out))  
#'
#'}
#'

scholz12<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
 # to keep things simple we will leave out the pegfilgrastim and leave out chemo toxicity, i.e. we look at only G-CSF
				G6=G6a+G6b+G6c+G6d+G6e
				MGB=G4a+G4b+G4c+G4d+G4e+G5a+G5b+G5c+G5d+G5e+G6
				G=CG+PGB+MGB

				arg=wG6*G6+wGRA*GRA/(wG6*G6nor+wGRA*GRAnor)
	         Pendo=Pendomax - (Pendomax-Pendomin)*((Pendomax-Pendonor)/(Pendomax-Pendomin))^(arg^Pendob)
				Pref= VFD*gref*(kFu+vGRAFmax/(kGRAFm+VFD*gref))
				dg=Pref*Pendo-kFu*g-vGRAFmax*g*GRA/GRAnor/(kGRAFm+g)
				gCenRel=g/gref
				dg1=gCenRel-DFGCSF*g1
				dg2=DFGCSF*(g1-g2)
				dg3=DFGCSF*(g2-g3)
				dg4=DFGCSF*(g3-g4)
				gCenRelDelay=DFGCSF*g4
				
				Srel=S/Snor
            x=wG*log(G/GSS)+wS*ifelse(Srel>1,Srel-1,log(Srel))
				yS=(-0.5/log(2))*(log((aintS-amaxS)/(aminS-aintS))-log((anorS-amaxS)/(aminS-anorS)))*x+
						0.5*log((anorS-amaxS)/(aminS-anorS))
				yCG=(-0.5/log(2))*(log((aintCG-amaxCG)/(aminCG-aintCG))-log((anorCG-amaxCG)/(aminCG-anorCG)))*x+
						0.5*log((anorCG-amaxCG)/(aminCG-anorCG))
				aS=(amaxS*exp(-yS)+aminS*exp(yS))/(exp(-yS)+exp(yS))
				aCG=(amaxCG*exp(-yCG)+aminCG*exp(yCG))/(exp(-yCG)+exp(yCG))
				thetaS=2*ifelse(Srel>1,1,1/Srel^0.6)
				p=pdelta*tanh(-thetaS*(Srel-1)-thetaG*(Srel-1))+0.5
#				psiCX=ifelse(chemoOn,1,0)    # skip toxicity for now
#				psiCX=ifelse(chemoOn&first,ffc,0) 
				dS = (2*p-1)*S*aS/TS  # - kS*psiCX*S  # skip tox
				Sout=2*(1-p)*S*aS/TS	
	
				ACG=AmaxCGF - (AmaxCGF-AminCGF)*((AmaxCGF-AnorCGF)/(AmaxCGF-AminCGF))^(gCenRelDelay^AbCGF)
				TCG=TmaxCGF - (TmaxCGF-TminCGF)*((TmaxCGF-TnorCGF)/(TmaxCGF-TminCGF))^(gCenRelDelay^AbCGF)
				dCG = Sout*ACGin - CG*ACG/TCG # - kCG*psiCX*CG  # skip chemotox
				CGout=ACGout*CG*ACG/TCG
				
				APGB=AmaxPGBF - (AmaxPGBF-AminPGBF)*((AmaxPGBF-AnorPGBF)/(AmaxPGBF-AminPGBF))^(gCenRelDelay^AbPGBF)
				TPGB=TmaxPGBF - (TmaxPGBF-TminPGBF)*((TmaxPGBF-TnorPGBF)/(TmaxPGBF-TminPGBF))^(gCenRelDelay^AbPGBF)
				dPGB = CGout*APGBin - PGB*APGB/TPGB #- kPGB*psiCX*PGB
				PGBout=APGBout*PGB*APGB/TPGB

				AG4=AnorG4F
				TG4=TmaxG4F - (TmaxG4F-TminG4F)*((TmaxG4F-TnorG4F)/(TmaxG4F-TminG4F))^(gCenRelDelay^AbG6F)
				dG4a =PGBout - G4a*5/TG4 # - kMGB*psiCX*G4a
				G4aout=G4a*5*AG4/TG4
				dG4b =G4aout - G4b*5/TG4 #- kMGB*psiCX*G4b
				G4bout=G4b*5*AG4/TG4
				dG4c =G4bout - G4c*5/TG4 #- kMGB*psiCX*G4c
				G4cout=G4c*5*AG4/TG4
				dG4d =G4cout - G4d*5/TG4 #- kMGB*psiCX*G4d
				G4dout=G4d*5*AG4/TG4
				dG4e =G4dout - G4e*5/TG4 #- kMGB*psiCX*G4e
				G4out=G4e*5*AG4/TG4

				AG5=AnorG5F
				TG5=TmaxG5F - (TmaxG5F-TminG5F)*((TmaxG5F-TnorG5F)/(TmaxG5F-TminG5F))^(gCenRelDelay^AbG6F)
				dG5a =G4out - G5a*5/TG5 #- kMGB*psiCX*G5a
				G5aout=G5a*5*AG5/TG5
				dG5b =G5aout - G5b*5/TG5 #- kMGB*psiCX*G5b
				G5bout=G5b*5*AG5/TG5
				dG5c =G5bout - G5c*5/TG5 # - kMGB*psiCX*G5c
				G5cout=G5c*5*AG5/TG5
				dG5d =G5cout - G5d*5/TG5 #- kMGB*psiCX*G5d
				G5dout=G5d*5*AG5/TG5
				dG5e =G5dout - G5e*5/TG5 #- kMGB*psiCX*G5e
				G5out=G5e*5*AG5/TG5

				AG6=AmaxG6F - (AmaxG6F-AminG6F)*((AmaxG6F-AnorG6F)/(AmaxG6F-AminG6F))^(gCenRelDelay^AbG6F)
				TG6=TmaxG6F - (TmaxG6F-TminG6F)*((TmaxG6F-TnorG6F)/(TmaxG6F-TminG6F))^(gCenRelDelay^AbG6F)
				dG6a =G5out - G6a*5/TG6 #- kMGB*psiCX*G6a
				G6aout=G6a*5*AG6/TG6
				dG6b =G6aout - G6b*5/TG6 #- kMGB*psiCX*G6b
				G6bout=G6b*5*AG6/TG6
				dG6c =G6bout - G6c*5/TG6 #- kMGB*psiCX*G6c
				G6cout=G6c*5*AG6/TG6
				dG6d =G6cout - G6d*5/TG6 #- kMGB*psiCX*G6d
				G6dout=G6d*5*AG6/TG6
				dG6e =G6dout - G6e*5/TG6 #- kMGB*psiCX*G6e
				G6out=G6e*5*AG6/TG6

#				psiPred=ifelse(predOn,1,0) 
				TGRA=TnorGRA #*(1+TpredGRA*psiPred)
				dGRA=G6out-GRA/TGRA
				
				return(list(c(dg,dg1,dg2,dg3,dg4,dS,dCG,dPGB,dG4a,dG4b,dG4c,dG4d,dG4e,dG5a,dG5b,dG5c,dG5d,dG5e,
										dG6a,dG6b,dG6c,dG6d,dG6e,dGRA),
							c(G=G,CG=CG,PGB=PGB)))
			})
}

