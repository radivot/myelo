#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
static double  parms[5];
#define tauS   parms[0]
#define fQ     parms[1]
#define thresh parms[2]
#define betaSS parms[3]
#define kapDel parms[4]

void parmsQdb(void (* odeparms)(int *, double *))
{   int N=5;
    odeparms(&N, parms);
}


void derivsQdb(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    double  Q, S1, S2, S3, S4, Aq;
    double dQ,dS1,dS2,dS3,dS4,dAq;
    double beta, k;
    if (ip[0] <1) error("nout should be at least 1");
    Q=y[0]; S1=y[1];S2=y[2];S3=y[3];S4=y[4]; Aq=y[5];   
    k=4.0/tauS;
    
    if (Q>thresh) {
        beta=betaSS;
    } else { 
        beta = fQ-((fQ-betaSS)/thresh)*Q;
    }
    dAq = 0;   // Changes only by events
    dQ=-(beta+kapDel)*Q + Aq*k*S4;
    dS1=beta*Q - k*S1;
    dS2=k*S1   - k*S2;
    dS3=k*S2   - k*S3;
    dS4=k*S3   - k*S4;
    
    ydot[0] = dQ ; 
    ydot[1] = dS1; 
    ydot[2] = dS2; 
    ydot[3] = dS3; 
    ydot[4] = dS4; 
    ydot[5] = dAq;
    
    yout[0] = beta;  
    
}

