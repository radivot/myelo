# NLMIXR project generated with IQRtools

# ==PROJECT HEADER START===================================================
# COMMENT             = 'BiExp CL  NLMIXR version'
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
# PARAMNAMES          = 'CL,Vc,Q1,Vp1'
# PARAMGUESS          = '10,10,0.5,10'
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
# THETANAMES          = 'CL,Vc,Q1,Vp1,error_PROP1'
# THETAESTIMATE       = '1,1,1,1,1'
# ETANAMES            = 'omega(CL),omega(Vc),omega(Q1),omega(Vp1)'
# ETAESTIMATE         = '1,1,1,1'
# CORRELATIONNAMES    = ''
# CORRESTIMATE        = ''
# NROBSERVATIONS      = '340'
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
    transCL <- log(10)
    transVc <- log(10)
    transQ1 <- log(0.5)
    transVp1 <- log(10)

# Initial guesses random effects and covariance structure
    ETA_CL ~ 0.25 # (variance)
    ETA_Vc ~ 0.25 # (variance)
    ETA_Q1 ~ 0.25 # (variance)
    ETA_Vp1 ~ 0.25 # (variance)

# Initial guesses error model parameters
    error_PROP1 <- 0.3
  })
  model({
# Parameter transformations
    CL <- exp(transCL + ETA_CL)
    Vc <- exp(transVc + ETA_Vc)
    Q1 <- exp(transQ1 + ETA_Q1)
    Vp1 <- exp(transVp1 + ETA_Vp1)


# Model variables
    Cc = Ac/Vc;
    OUTPUT1 = Cc;

# Model differential equations
    d/dt(Ac) = -Q1/Vc*Ac+Q1/Vp1*Ap1-CL/Vc*Ac;
    d/dt(Ap1) = Q1/Vc*Ac-Q1/Vp1*Ap1;

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
nlmixr.control$mcmc$nmc <- 3

# ----------------------------------------
# Run the NLMIXR model
# ----------------------------------------
fit <- nlmixr::nlmixr(modelNLMIXR,dataNLMIXR,est=nlmixr.method,control=nlmixr.control,table=nlmixr::tableControl(npde=TRUE,cwres=TRUE))

# ----------------------------------------
# Store fit as Rdata object
# ----------------------------------------
save(fit,file="RESULTSORIG/project.fit")

