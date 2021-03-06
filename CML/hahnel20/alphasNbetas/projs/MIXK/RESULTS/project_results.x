list(type = "NLMIXR", method = "SAEM", path = structure("projs/MIXK", class = "IQRnlmeProject", absProjectPath = "/IQDESKTOP/SHARE/projs/MIXK", absModelPath = "/IQDESKTOP/SHARE/projs/MIXK/model.txt", absDataPath = "/IQDESKTOP/SHARE/data/biExport.csv"), 
    objectivefunction = list(OBJ = 91.8117190399498, AIC = 109.81171903995, 
        BIC = 144.272229598442), residualerrormodels = "rel", 
    trans_randeffects = c("exp(phi)", "exp(phi)", "exp(phi)", 
    "exp(phi)"), inv_trans_randeffects = c("log(psi)", "log(psi)", 
    "log(psi)", "log(psi)"), PROJECTINFO = list(COMMENT = "BiExp Ke  NLMIXR version", 
        TOOL = "NLMIXR", FILE = "project_NLMIXR.R", METHOD = "SAEM", 
        DATA = "../../data/biExport.csv", DOSINGTYPES = "BOLUS", 
        TK0NAMES = "NA", COVNAMES = "", CATNAMES = "", REGRESSIONNAMES = "", 
        OUTPUTS = "Cc", ERRORMODELS = "rel", ERRORNAMES = "error_PROP1", 
        PARAMNAMES = c("Vc", "Ke", "Kpc", "Kcp"), PARAMGUESS = c("10", 
        "10", "0.05", "0.05"), PARAMESTIMATE = c("1", "1", "1", 
        "1"), PARAMTRANS = c("exp(phi)", "exp(phi)", "exp(phi)", 
        "exp(phi)"), PARAMINVTRANS = c("log(psi)", "log(psi)", 
        "log(psi)", "log(psi)"), COVARIATENAMES = "", COVARIATESUSED = "", 
        BETACOVNAMES = "", BETACOVTRANS = "", BETACATNAMES = "", 
        BETACATREFERENCE = "", BETACATCATEGORIES = "", THETANAMES = c("Vc", 
        "Ke", "Kpc", "Kcp", "error_PROP1"), THETAESTIMATE = c("1", 
        "1", "1", "1", "1"), ETANAMES = c("omega(Vc)", "omega(Ke)", 
        "omega(Kpc)", "omega(Kcp)"), ETAESTIMATE = c("1", "1", 
        "1", "1"), CORRELATIONNAMES = "", CORRESTIMATE = "", 
        NROBSERVATIONS = "340", NRPARAM_ESTIMATED = "9"), parameters = list(
        names = c(Vc = "Vc", Ke = "Ke", Kpc = "Kpc", Kcp = "Kcp", 
        error_PROP1 = "error_PROP1"), FLAGestimated = c(1, 1, 
        1, 1, 1), transformation = c("exp(phi)", "exp(phi)", 
        "exp(phi)", "exp(phi)", ""), values = c(Vc = 2.30699612146821, 
        Ke = -0.102902491476166, Kpc = -3.7421719152121, Kcp = -4.38204301441174, 
        error_PROP1 = 0.909161666811688), stderrors = c(0.359831588903864, 
        0.117799923253774, 0.360645078258779, 0.486296729350955, 
        0), covariancematrix = structure(c(0.129478772373079, 
        -0.0131319025363414, -0.00101116224543276, 0.0352315415399722, 
        0, -0.0131319025363414, 0.0138768219185951, 0.000947997593902832, 
        -0.00492305052779681, 0, -0.00101116224543276, 0.000947997593902832, 
        0.130064872472281, -0.0246205471998219, 0, 0.0352315415399722, 
        -0.00492305052779681, -0.0246205471998219, 0.236484508977436, 
        0, 0, 0, 0, 0, 0), .Dim = c(5L, 5L)), correlationmatrix = structure(c(1, 
        -0.309801386244893, -0.00779186733104176, 0.201340438783613, 
        0, -0.309801386244893, 1, 0.0223142460638847, -0.0859385332255897, 
        0, -0.00779186733104176, 0.0223142460638847, 1, -0.140383589449516, 
        0, 0.201340438783613, -0.0859385332255897, -0.140383589449516, 
        1, 0, 0, 0, 0, 0, 1), .Dim = c(5L, 5L))), rawParameterInfo = list(
        fixedEffects = list(names = c(Vc = "Vc", Ke = "Ke", Kpc = "Kpc", 
        Kcp = "Kcp"), trans = c("", "", "", ""), invtrans = c("", 
        "", "", ""), estimated = c(1, 1, 1, 1), values = c(10.0442077138041, 
        0.902214942839823, 0.0237025672194387, 0.0124997952549536
        ), stderr = c(Vc = 3.61422322093857, Ke = 0.106280851024939, 
        Kpc = 0.00854821420978844, Kcp = 0.00607860955004053), 
            rse = c(Vc = 35.9831588903864, Ke = 11.7799923253774, 
            Kpc = 36.0645078258779, Kcp = 48.6296729350955), 
            distribution_info = c("log(psi)", "log(psi)", "log(psi)", 
            "log(psi)")), fixedEffects_transformed = list(names = c(Vc = "Vc", 
        Ke = "Ke", Kpc = "Kpc", Kcp = "Kcp"), trans = c("exp(phi)", 
        "exp(phi)", "exp(phi)", "exp(phi)"), invtrans = c("log(psi)", 
        "log(psi)", "log(psi)", "log(psi)"), estimated = c(1, 
        1, 1, 1), values = c(Vc = 2.30699612146821, Ke = -0.102902491476166, 
        Kpc = -3.7421719152121, Kcp = -4.38204301441174), stderr = c(0.359831588903864, 
        0.117799923253774, 0.360645078258779, 0.486296729350955
        ), rse = c(Vc = 15.5974076226388, Ke = 114.477231371077, 
        Kpc = 9.63731988882553, Kcp = 11.0974887227627), distribution_info = c("log(psi)", 
        "log(psi)", "log(psi)", "log(psi)")), randomEffects = list(
            names = c("omega(Vc)", "omega(Ke)", "omega(Kpc)", 
            "omega(Kcp)"), estimated = c(1, 1, 1, 1), values = c(1.18374558992207, 
            0.456784238402092, 1.36468851864467, 1.95615707263334
            ), covariance = structure(c(1.40125362165996, 0, 
            0, 0, 0, 0.208651840452579, 0, 0, 0, 0, 1.86237475292059, 
            0, 0, 0, 0, 3.82655049281342), .Dim = c(4L, 4L)), 
            stderr = c(NA, NA, NA, NA), rse = c(NA, NA, NA, NA
            )), errorParameter = list(names = c(error_PROP1 = "error_PROP1"), 
            estimated = 1, values = c(error_PROP1 = 0.909161666811688), 
            stderr = NA, rse = NA)))
