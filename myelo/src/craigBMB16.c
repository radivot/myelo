#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
static double parms[43];
#define Qss parms[0]
#define Nrss parms[1]
#define Nss parms[2]
#define G1ss parms[3]
#define G2ss parms[4]
#define Aqss parms[5]
#define Anss parms[6]
#define gamNMss parms[7] //death rates
#define gamNrss parms[8]
#define gamNss parms[9]
#define tauS parms[10]
#define tauNP parms[11]
#define fQ   parms[12]
#define the2   parms[13]
#define s2    parms[14]
#define kapss   parms[15]
#define kapDel parms[16]
#define etaNPss   parms[17]  //growth rates
#define etaNPmin    parms[18]
#define bNP parms[19]
#define ka parms[20]
#define kren parms[21]
#define kint parms[22]
#define Gprod parms[23]
#define k12g parms[24]
#define k21g parms[25]
#define Pow parms[26]
#define Vmax parms[27]
#define bV parms[28]
#define phiNrss parms[29]
#define phiNrmax parms[30]
#define bG parms[31]
#define V parms[32]
#define k12 parms[33]
#define k21 parms[34]
#define k24 parms[35]
#define k42 parms[36]
#define k13 parms[37]
#define k31 parms[38]
#define kelC parms[39]
#define hQ parms[40]
#define EM50 parms[41]
#define Sc parms[42]

/* Interface to dede utility functions in package deSolve */
void lagvalue(double T, int *nr, int N, double *ytau) {
    static void(*fun)(double, int*, int, double*) = NULL;
    if (fun == NULL)
        fun = (void(*)(double, int*, int, double*))R_GetCCallable("deSolve", "lagvalue");
    return fun(T, nr, N, ytau);
}

/*
 void lagderiv(double T, int *nr, int N, double *ytau) {
 static void(*fun)(double, int*, int, double*) = NULL;
 if (fun == NULL)
 fun = (void(*)(double, int*, int, double*))R_GetCCallable("deSolve", "lagderiv"); return fun(T, nr, N, ytau);
 } */


void parmsCraig16(void (* odeparms)(int *, double *))
{   int N=43;
    odeparms(&N, parms);
}


void derivsCraig16(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    double Q,Nr,N,G1,G2,Tn,An,Aq,Cp,Cf,Cs1,Cs2,Gs,Ic,Ig;
    double dQ,dNr,dN,dG1,dG2,dTn,dAn,dAq,dCp,dCf,dCs1,dCs2,dGs,dIc,dIg;
    double Vn,VnTnm,Vrat,GBF,GBFss,phiNr;
    double beta,betaTs,eta,etaTn,etaTnm, tauNM; //,etaNP,etaNPchemo,Vn,GBF,GBFss,phiNr;
    if (ip[0] <1) error("nout should be at least 1");
    Q  =y[0]; Nr =y[1];N  =y[2]; G1 =y[3]; G2 =y[4];     Tn =y[5]; An =y[6];Aq =y[7];
    Cp =y[8]; Cf =y[9];Cs1=y[10];Cs2=y[11];              Gs =y[12];Ic=y[13];Ig =y[14];
    
    int NoutTs  = 2;  //number of returned lags for tauS
    int nrTs[2] = {0, 8}; // positions of Q and Cp  
//    double ytauTs[2] = {Qss,0.0};
    double ytauTs[2] = {0.0,0.0};
    
    int NoutTn  = 3;  //number of returned lags for tauN (i.e. state variable Tn)
    int nrTn[3] = {0, 3, 8}; // positions of Q and Cp  
 //   double ytauTn[3] = {Qss,G1ss,0.0};
    double ytauTn[3] = {0.0,0.0,0.0};
    
    int NoutTnm  = 2;  //number of returned lags for tauS
    int nrTnm[2] = {3, 8}; // positions of Q and Cp  
 //   double ytauTnm[2] = {G1ss,0.0};
    double ytauTnm[2] = {0.0,0.0};
    double ts,tn,tnm;
    double Qts,Qtn,G1tn,G1tnm,Cpts,Cptn,Cptnm;
    
    beta=fQ/(1+pow(Q/the2,s2));
 //   etaNP=etaNPss  + (etaNPss -etaNPmin)*(bNP/G1ss)*(G1-G1ss)/(G1+bNP);
 //   etaNPchemo=etaNP/(1+pow(Cp/EM50,Sc));
    eta=etaNPss  + (etaNPss - etaNPmin)*(bNP/G1ss)*(G1-G1ss)/(G1+bNP);
    eta=eta/(1+pow(Cp/EM50,Sc));
    GBF=G2/(V*(Nr+N));
    GBFss=G2ss/(V*(Nrss+Nss));
    phiNr=phiNrss  + (phiNrmax-phiNrss)*(GBF-GBFss)/(GBF-GBFss+bG);
    dIc=0; 
    dIg=0; // same for GCSF.  Changes only by events
    if (*t<0) {
        dQ=-(beta+kapss+kapDel)*Q + Aq*beta*Q;
        dNr=-gamNrss*Nr-phiNr*Nr + (An/1e3)*kapss*Q;
        dN= phiNr*Nr-gamNss*N;
        dG1=Gprod - kren*G1 - k12g*((Nr+N)*V-G2)*pow(G1,Pow)+k21g*G2; //unbound circulating GCSF
        dG2=k12g*((Nr+N)*V-G2)*pow(G1,Pow) -kint*G2 - k21g*G2;   // bound GCSF
        dTn = 0;    
        dAn = 0; 
        etaTnm=eta;
        etaTn=eta;
        dAq = 0;    
        dCp = 0;    
        dCf = 0;    
        dCs1= 0;    
        dCs2= 0;    
        dGs = 0; //skin pool that feeds into G1  
        
        ytauTs[0]=Qss; 
        ytauTs[1]=0;
        ytauTn[0]=Qss;
        ytauTn[1]=G1ss;
        ytauTn[2]=0;
        ytauTnm[0]=G1ss;
        ytauTnm[1]=0;
        
        Qts=ytauTs[0]; 
        Cpts=ytauTs[1];
        Qtn=ytauTn[0];
        G1tn=ytauTn[1];
        Cptn=ytauTn[2];
        G1tnm=ytauTnm[0];
        Cptnm=ytauTnm[1];
        
        
        
    }	else {
        
//        ts= *t - tauS;
        ts= *t - round(1000*tauS)/1000;
//        ts= round(1000*(*t - tauS))/1000;
//        tn= *t - Tn;
        tauNM = Tn - tauNP; //  Tn is time to reserves
//        tnm= *t - tauNM;
    //    tn= round(100*(*t - Tn))/100;
        tn= *t - round(1000*Tn)/1000;
        tnm= *t - round(1000*tauNM)/1000;
//        tnm= round(100*(*t - tauNM))/100;
        
        lagvalue(ts, nrTs, NoutTs, ytauTs);
        lagvalue(tn, nrTn, NoutTn, ytauTn);
        lagvalue(tnm, nrTnm, NoutTnm, ytauTnm);
        
        Qts=ytauTs[0]; 
        Cpts=ytauTs[1];
        Qtn=ytauTn[0];
        G1tn=ytauTn[1];
        Cptn=ytauTn[2];
        G1tnm=ytauTnm[0];
        Cptnm=ytauTnm[1];
        betaTs=fQ/(1+pow(Qts/the2,s2));
        dQ=-(beta+kapss+kapDel)*Q + Aq*betaTs*Qts;
        Vn= 1 + (Vmax-1)*(G1-G1ss)/(G1-G1ss+bV);
        VnTnm= 1 + (Vmax-1)*(G1tnm-G1ss)/(G1tnm-G1ss+bV);
        Vrat=Vn/VnTnm;
        dTn=1-Vrat;
        dNr=(An/1e3)*kapss*Qtn*Vrat-gamNrss*Nr-phiNr*Nr; // 1e3 maps units of Q to N
        dN = phiNr*Nr-gamNss*N;
        
        dG1=Gprod + ka*Gs + Ig - kren*G1 - k12g*((Nr+N)*V-G2)*pow(G1,Pow) + k21g*G2;
        dG2=k12g*((Nr+N)*V-G2)*pow(G1,Pow) - kint*G2 - k21g*G2;
        
        etaTn=etaNPss  + (etaNPss -etaNPmin)*(bNP/G1ss)*(G1tn-G1ss)/(G1tn+bNP);
        etaTn=etaTn/(1+pow(Cptn/EM50,Sc));
        etaTnm=etaNPss  + (etaNPss -etaNPmin)*(bNP/G1ss)*(G1tnm-G1ss)/(G1tnm+bNP);
        etaTnm=etaTnm/(1+pow(Cptnm/EM50,Sc));
        
        dAn=An*((1-dTn)*(etaTnm-etaTn)-gamNMss*dTn);
        dAq=Aq*hQ*(Cpts-Cp);
        dCp=k21*Cf+k31*Cs1-(k12+k13+kelC)*Cp  + Ic;
        dCf=k12*Cp+k42*Cs2-(k21+k24)*Cf;
        dCs1=k13*Cp-k31*Cs1;
        dCs2=k24*Cf-k42*Cs2;
        dGs=-ka*Gs + Ig; // G in skin, add F*Dg/Vd to Gs at each injection
    }
    ydot[0] = dQ  ; ydot[1] = dNr ; ydot[2] = dN  ; ydot[3] = dG1 ; ydot[4] = dG2 ;
    ydot[5] = dTn ; ydot[6] = dAn ; ydot[7] = dAq ;
    ydot[8] = dCp ; ydot[9] = dCf ; ydot[10]= dCs1; ydot[11]= dCs2;
    ydot[12]= dGs ; ydot[13]= dIc ; ydot[14]= dIg ;
    yout[0] = N*8.19;  // N in 1e9 cell/kg, ANC in 1e3 cells/uL
    yout[1] = Qts;  
    yout[2] = Cpts;
    yout[3] = Qtn;
    yout[4] = G1tn;
    yout[5] = Cptn;
    yout[6] = G1tnm;
    yout[7] = Cptnm;
    yout[8] = dAn;
    yout[9] = dTn;
    yout[10] = eta;
    yout[11] = etaTnm;
    yout[12] = etaTn;
    yout[13] = dTn*gamNMss;
}

