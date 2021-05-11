# NLMIXR project generated with IQRtools

# ==PROJECT HEADER START===================================================
# COMMENT             = 'BiExp Ke  NLMIXR version'
# TOOL                = 'NLMIXR'
# FILE                = 'project_NLMIXR.R'
# METHOD              = 'SAEM'
# DATA                = '../../data/biExport.csv'
# DOSINGTYPES         = 'BOLUS'
# TK0NAMES            = 'NA'
# COVNAMES            = ''
# CATNAMES            = ''
# REGRESSIONNAMES     = ''
# OUTPUTS             = 'Cc'
# ERRORMODELS         = 'rel'
# ERRORNAMES          = 'error_PROP1'
# PARAMNAMES          = 'Vc,Ke,Kpc,Kcp'
# PARAMGUESS          = '10,0.05,0.005,0.005'
# PARAMESTIMATE       = '1,1,1,1'
# PARAMTRANS          = 'exp(phi),exp(phi),exp(phi),exp(phi)'
# PARAMINVTRANS       = 'log(psi),log(psi),log(psi),log(psi)'
# COVARIATENAMES      = ''
# COVARIATESUSED      = ''
# BETACOVNAMES        = ''
# BETACOVTRANS        = ''
# BETACATNAMES        = ''
# BETACATREFERENCE    = ''
# BETACATCATEGORIES   = ''
# THETANAMES          = 'Vc,Ke,Kpc,Kcp,error_PROP1'
# THETAESTIMATE       = '1,1,1,1,1'
# ETANAMES            = 'omega(Vc),omega(Ke),omega(Kpc),omega(Kcp)'
# ETAESTIMATE         = '1,1,1,1'
# CORRELATIONNAMES    = ''
# CORRESTIMATE        = ''
# NROBSERVATIONS      = '435'
# NRPARAM_ESTIMATED   = '9'
# ==PROJECT HEADER END=====================================================

# ----------------------------------------
# Load dataset and convert to RxODE format
# ----------------------------------------
dataNLMIXR <- IQRtools::convertDataRxODE_IQRnlmixr("../../data/biExport.csv")
IQRtools::IQRsaveCSVdata(dataNLMIXR,"dataNLMIXR.csv")

# ----------------------------------------
# Define the NLMIXR model function
# ----------------------------------------
modelNLMIXR <- function() {
  ini({
# Initial guesses fixed effects
    transVc <- log(10)
    transKe <- log(0.05)
    transKpc <- log(0.005)
    transKcp <- log(0.005)

# Initial guesses random effects and covariance structure
    ETA_Vc ~ 0.25 # (variance)
    ETA_Ke ~ 0.25 # (variance)
    ETA_Kpc ~ 0.25 # (variance)
    ETA_Kcp ~ 0.25 # (variance)

# Initial guesses error model parameters
    error_PROP1 <- 0.3
  })
  model({
# Parameter transformations
    Vc <- exp(transVc + ETA_Vc)
    Ke <- exp(transKe + ETA_Ke)
    Kpc <- exp(transKpc + ETA_Kpc)
    Kcp <- exp(transKcp + ETA_Kcp)


# Model variables
    Cc = Ac/Vc;
    OUTPUT1 = Cc;

# Model differential equations
    d/dt(Ac) = -Kcp*Ac+Kpc*Ap1-Ke*Ac;
    d/dt(Ap1) = Kcp*Ac-Kpc*Ap1;

# Error model
    OUTPUT1 ~ prop(error_PROP1)
  })
}

# ----------------------------------------
# Define NLMIXR method/control settings
# ----------------------------------------
nlmixr.method <- "saem"
nlmixr.control <- nlmixr::saemControl()
nlmixr.control$mcmc$niter[1] <- 50
nlmixr.control$mcmc$niter[2] <- 20
nlmixr.control$seed <- 123456
nlmixr.control$mcmc$nmc <- 1

# ----------------------------------------
# Run the NLMIXR model
# ----------------------------------------
fit <- nlmixr::nlmixr(modelNLMIXR,dataNLMIXR,est=nlmixr.method,control=nlmixr.control,table=nlmixr::tableControl(npde=TRUE,cwres=TRUE))

# ----------------------------------------
# Store fit as Rdata object
# ----------------------------------------
save(fit,file="RESULTSORIG/project.fit")

