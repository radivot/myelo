#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
static double parms[16];
#define Qss parms[0]
#define Aqss parms[1]
#define tauS parms[2]
#define fQ   parms[3]
#define the2   parms[4]
#define s2    parms[5]
#define kapss   parms[6]
#define kapDel parms[7]
#define k12 parms[8]
#define k21 parms[9]
#define k24 parms[10]
#define k42 parms[11]
#define k13 parms[12]
#define k31 parms[13]
#define kelC parms[14]
#define hQ parms[15]


void lagvalue(double T, int *nr, int N, double *ytau); // defined in craigBMB16.c

/* Interface to dede utility functions in package deSolve  
void lagvalue2(double T, int *nr, int N, double *ytau) {
    static void(*fun)(double, int*, int, double*) = NULL;
    if (fun == NULL)
        fun = (void(*)(double, int*, int, double*))R_GetCCallable("deSolve", "lagvalue");
    return fun(T, nr, N, ytau);
}
*/


void parmsCraig16Q(void (* odeparms)(int *, double *))
{   int N=16;
    odeparms(&N, parms);
}


void derivsCraig16Q(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    double  Q, Aq, Cp, Cf, Cs1, Cs2, Ic;
    double dQ,dAq,dCp,dCf,dCs1,dCs2,dIc;
    double beta,betaTs; //,etaNP,etaNPchemo,Vn,GBF,GBFss,phiNr;
    if (ip[0] <1) error("nout should be at least 1");
    Q=y[0]; Aq=y[1];    Cp=y[2]; Cf=y[3]; Cs1=y[4]; Cs2=y[5]; Ic=y[6];
    
    int NoutTs  = 2;  //number of returned lags for tauS
    int nrTs[2] = {0, 2}; // positions of Q and Cp  
//    double ytauTs[2] = {Qss,0.0};
    double ytauTs[2] = {0.0,0.0};
    
    double ts;
    double Qts,Cpts;
    
    beta=fQ/(1+pow(Q/the2,s2));

    dIc=0; // Changes only by events
    if (*t<0) {
        dQ=-(beta+kapss+kapDel)*Q + Aq*beta*Q;
        dAq = 0;  
        
        dCp = 0;    
        dCf = 0;    
        dCs1= 0;    
        dCs2= 0;    
        
        ytauTs[0]=Qss; 
        ytauTs[1]=0;
        
        Qts=ytauTs[0]; 
        Cpts=ytauTs[1];
        
    }	else {
//        ts= *t - tauS;
        ts= *t - round(1e4*tauS)/1e4;
//        lagvalue2(ts, nrTs, NoutTs, ytauTs);
        lagvalue(ts, nrTs, NoutTs, ytauTs);

        Qts=ytauTs[0]; 
        Cpts=ytauTs[1];

        betaTs=fQ/(1+pow(Qts/the2,s2));
        dQ=-(beta+kapss+kapDel)*Q + Aq*betaTs*Qts;
        dAq=Aq*hQ*(Cpts-Cp);
        dCp=k21*Cf+k31*Cs1-(k12+k13+kelC)*Cp  + Ic;
        dCf=k12*Cp+k42*Cs2-(k21+k24)*Cf;
        dCs1=k13*Cp-k31*Cs1;
        dCs2=k24*Cf-k42*Cs2;
    }
    ydot[0] = dQ  ; ydot[1] = dAq ;
    ydot[2] = dCp ; ydot[3] = dCf ; ydot[4] = dCs1; ydot[5] = dCs2;
    ydot[6] = dIc ; 

    yout[0] = Qts;  
    yout[1] = Cpts;
}

