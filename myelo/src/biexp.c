#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
static double  parms[2];
#define Ke    parms[0]
#define Kpc   parms[1]


void parmsBiExp(void (* odeparms)(int *, double *))
{   int N=2;
    odeparms(&N, parms);
}


void derivsBiExp(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    double Ac, Ap;
    double dAc,dAp;

    if (ip[0] <1) error("nout should be at least 1");
    
    Ac=y[0];  Ap=y[1];    
    
    dAc = Kpc*Ap - Ke*Ac;
    dAp =  - Kpc*Ap;
    ydot[0] = dAc; 
    ydot[1] = dAp;

    yout[0]=(Ac+Ap)/(Ac+Ap+2);
}
