#include <R.h>   // qdb1.c with parameters reordered and renamed (for beta0)
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
static double  parms[4]; 
#define beta0  parms[0]
#define betaSS parms[1]
#define kapDel parms[2]
#define thresh parms[3]

void parmsQdb1r(void (* odeparms)(int *, double *))
{   int N=4;
    odeparms(&N, parms);
}


void derivsQdb1r(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    double  Q,  Aq;
    double dQ, dAq;
    double beta;
    if (ip[0] <1) error("nout should be at least 1");
    Q=y[0];  Aq=y[1];   
    
    if (Q>thresh) {
        beta=betaSS;
    } else { 
        beta = beta0-((beta0-betaSS)/thresh)*Q;
    }
    dAq = 0;   // Changes only by events
    dQ=-kapDel*Q + (Aq-1)*beta*Q;

    ydot[0] = dQ ; 
    ydot[1] = dAq;
    
    yout[0] = beta;  
    
}

