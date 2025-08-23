## setup.R
gl=geom_line()
gx=xlab("Years Since CML Diagnosis")
gxTFR=xlab("Years Since TFR Attempt")
gx0=xlab("Years")

gy=ylab("BCR::ABL1/ABL1")
gyIS=ylab("BCR::ABL1/ABL1(%)")
gyP=ylab("Probability of Death by CML")
gyE=ylab("Excess Absolute Risk of Death")

tc=function(sz) theme_classic(base_size=sz)
top=theme(legend.position="top",legend.title = element_text(size=10))
sbb=theme(strip.background=element_blank()) 
geE=geom_errorbar(aes(ymin=LL,ymax=UL),width=0.2)#for absolute risks
bks=c(0.01,0.1,1,10,100)
sy=scale_y_log10(breaks=bks,labels=bks)
ssm=scale_shape_manual(values=c(16,1)) 
cc1=coord_cartesian(ylim=c(1e-2,100))
layout <- "
     A
     B
     B
     "

he=function(t,k1=0.35121676,k2=1.11192339,k3=1.61921238) {# from fit in Fig_2AB.R
  x1=exp(-k1*t)  
  x2=k1*(exp(-k1*t)-exp(-k2*t))/(k2-k1)
  x3B=k1*k2*(exp(-k1*t)/((k2-k1)*(k3-k1))+exp(-k2*t)/((k1-k2)*(k3-k2))+exp(-k3*t)/((k1-k3)*(k2-k3)))
  denom=(x1+x2+x3B)
  h2=k3*x3B/denom
  tibble(t,h2,x1,x2,x3B)
} #uses analytic excess hazard solutions valid for L = constant

#solve excess hazard markov model numerically when L is time varying in Figs. 3 & 4
pars=c(k1=0.35121676,k2=1.11192339,k3=1.61921238)/12 # updated to April 2025 values 
markov<-function(Time, State, Pars) { # time is in months here, in he time is in years 
  with(as.list(c(Time, State, Pars)),{
    dL   = 0
    dX1  = -k1*L*X1
    dX2  = k1*L*X1 - k2*X2
    dX3b = k2*X2-k3*X3b # H is the cumulative hazard of progression to death 
    dH   = k3*X3b/(X1+X2+X3b) # conditional on it not having happened yet  
    list(c(dL,dX1,dX2,dX3b,dH),c(P=1-exp(-H)))# P=1-RS where RS=exp(-H)
  }) #H is mean # of lethal hits (via CML), exp(-H) is Poisson prob of zero
}




