#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
static double parms[7];
#define Qss parms[0]
#define Aqss parms[1]
#define tauS parms[2]
#define fQ   parms[3]
#define the2   parms[4]
#define s2    parms[5]
#define kapDel parms[6] 



void lagvalue(double T, int *nr, int N, double *ytau); // defined in craigBMB16.c



void parmsCraig16Qnc(void (* odeparms)(int *, double *))
{   int N=7;
    odeparms(&N, parms);
}


void derivsCraig16Qnc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    double  Q, Aq;
    double dQ,dAq;
    double beta,betaTs; //,etaNP,etaNPchemo,Vn,GBF,GBFss,phiNr;
    if (ip[0] <1) error("nout should be at least 1");
    Q=y[0]; Aq=y[1];   
    
//    int NoutTs  = 2;  //number of returned lags for tauS
    int NoutTs  = 1;  //reducing to one state grab back in time made no diff
//    int nrTs[2] = {0, 1}; // positions of Q and Aq  
    int nrTs[1] = {0}; // positions of Q and Aq  
//    double ytauTs[2] = {Qss,Aqss};
//    double ytauTs[2] = {0.0,0.0};
    double ytauTs[1] = {0.0};
    
    double ts;
    double Qts;
//    double Aqts;
    
    beta=fQ/(1+pow(Q/the2,s2));

    dAq = 0;   // Changes only by events
    if (*t<0) {
        dQ=-(beta+kapDel)*Q + Aq*beta*Q;
        ytauTs[0]=Qss; 
//        ytauTs[1]=Aqss;
        
        Qts=ytauTs[0]; 
//        Aqts=ytauTs[1];
        
    }	else {
//        ts= *t - tauS;
        ts= *t - round(1e4*tauS)/1e4;
//        lagvalue2(ts, nrTs, NoutTs, ytauTs);
        lagvalue(ts, nrTs, NoutTs, ytauTs);

        Qts=ytauTs[0]; 
//        Aqts=ytauTs[1];

        betaTs=fQ/(1+pow(Qts/the2,s2));
//        dQ=-(beta+kapDel)*Q + ((Aq+Aqts)/2)*betaTs*Qts;
//        dQ=-(beta+kapDel)*Q + Aqts*betaTs*Qts;
        dQ=-(beta+kapDel)*Q + Aq*betaTs*Qts;
    }
    ydot[0] = dQ  ; 
    ydot[1] = dAq ;

    yout[0] = Qts;  
//    yout[1] = Aqts;
}

