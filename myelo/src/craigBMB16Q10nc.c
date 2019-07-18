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



void parmsCraig16Q10nc(void (* odeparms)(int *, double *))
{   int N=7;
    odeparms(&N, parms);
}


void derivsCraig16Q10nc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    double  Q, S1, S2, S3, S4, S5, S6, S7, S8, S9, Aq;
    double dQ,dS1,dS2,dS3,dS4,dS5,dS6,dS7,dS8,dS9,dAq;
    double beta, k;
    if (ip[0] <1) error("nout should be at least 1");
    Q=y[0]; S1=y[1];S2=y[2];S3=y[3];S4=y[4];S5=y[5];S6=y[6];S7=y[7];S8=y[8];S9=y[9]; Aq=y[10];   
    k=9.0/tauS;
    
    beta=fQ/(1+pow(Q/the2,s2));
    dAq = 0;   // Changes only by events
    dQ=-(beta+kapDel)*Q + Aq*k*S9;
    dS1=beta*Q - k*S1;
    dS2=k*S1 - k*S2;
    dS3=k*S2 - k*S3;
    dS4=k*S3 - k*S4;
    dS5=k*S4 - k*S5;
    dS6=k*S5 - k*S6;
    dS7=k*S6 - k*S7;
    dS8=k*S7 - k*S8;
    dS9=k*S8 - k*S9;
    
    ydot[0] = dQ ; 
    ydot[1] = dS1; 
    ydot[2] = dS2; 
    ydot[3] = dS3; 
    ydot[4] = dS4; 
    ydot[5] = dS5; 
    ydot[6] = dS6; 
    ydot[7] = dS7; 
    ydot[8] = dS8; 
    ydot[9] = dS9; 
    ydot[10] = dAq;
    
    yout[0] = beta;  
    
}

