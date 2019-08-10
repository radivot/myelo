#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
static double  parms[11];
#define Circ0  parms[0]
#define ktr    parms[1]
#define gam    parms[2]
#define slope  parms[3]
#define k12    parms[4]
#define k21    parms[5]
#define k13    parms[6]
#define k31    parms[7]
#define k10    parms[8]
#define V1     parms[9]
#define mw     parms[10]


void parmsFri02(void (* odeparms)(int *, double *))
{   int N=11;
    odeparms(&N, parms);
}


void derivsFri02(int *neq, double *t, double *y, double *ydot)
{
    double C1, C2,C3,Prol,Trans1,Trans2,Trans3,Circ, Edrug, Cp;
    double dC1, dC2,dC3,dProl,dTrans1,dTrans2,dTrans3,dCirc;
    C1=y[0];  C2=y[1];   C3=y[2];  Prol=y[3];   
    Trans1=y[4];  Trans2=y[5];   Trans3=y[6];  Circ=y[7];   
    
    dC1=-(k12+k13+k10)*C1 + k21*C2 + k31*C3;
    dC2=k12*C1 - k21*C2;
    dC3=k13*C1 - k31*C3;
    Cp=C1/V1/mw;
    Edrug=slope*Cp;
    dProl = ktr*Prol*(1-Edrug)*pow(Circ0/Circ,gam) - ktr*Prol;
    dTrans1=ktr*Prol-ktr*Trans1;
    dTrans2=ktr*Trans1-ktr*Trans2;
    dTrans3=ktr*Trans2-ktr*Trans3;
    dCirc=ktr*Trans3-ktr*Circ;
    
    ydot[0] = dC1 ; 
    ydot[1] = dC2;
    ydot[2] = dC3;
    ydot[3] = dProl;
    ydot[4] = dTrans1;
    ydot[5] = dTrans2;
    ydot[6] = dTrans3;
    ydot[7] = dCirc;
}
