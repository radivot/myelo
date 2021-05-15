#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
static double  parms[15];
#define d0    parms[0]
#define d1    parms[1]
#define d2    parms[2]
#define d3    parms[3]
#define ax    parms[4]
#define bx    parms[5]
#define cx    parms[6]
#define D0a   parms[7]
#define D0b   parms[8]
#define az    parms[9]
#define bz    parms[10]
#define rz    parms[11]
#define ry    parms[12]
#define u     parms[13]
#define sp    parms[14]


void parmsMichor(void (* odeparms)(int *, double *))
{   int N=15;
    odeparms(&N, parms);
}


void derivsMichor(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    double X0, X1, X2, X3, Y0, Y1, Y2, Y3, Z0, Z1, Z2, Z3, D;
    double dX0,dX1,dX2,dX3,dY0,dY1,dY2,dY3,dZ0,dZ1,dZ2,dZ3, dD;
    double lambda;
    
    if (ip[0] <1) error("nout should be at least 1");
    
    X0=y[0];  X1=y[1];   X2=y[2];  X3=y[3];   
    Y0=y[4];  Y1=y[5];   Y2=y[6];  Y3=y[7];   
    Z0=y[8];  Z1=y[9];   Z2=y[10]; Z3=y[11];   
    
    lambda=-0.5*(X0-sp);
    dX0 = (lambda-d0)*X0;
    dX1 = ax*X0-d1*X1;
    dX2 = bx*X1-d2*X2;
    dX3 = cx*X2-d3*X3;
    dY0 = (ry*(1-u)-d0)*Y0;
    dY1 = Y0*az/(1+D/D0a) - d1*Y1;
    dY2 = Y1*bz/(1+D/D0b) - d2*Y2;
    dY3 = cx*Y2-d3*Y3;
    dZ0 = (rz-d0)*Z0+ry*u*Y0;
    dZ1 = az*Z0-d1*Z1;
    dZ2 = bz*Z1-d2*Z2;
    dZ3 = cx*Z2-d3*Z3;
    dD = 0;
    ydot[0] = dX0; 
    ydot[1] = dX1;
    ydot[2] = dX2;
    ydot[3] = dX3;
    ydot[4] = dY0;
    ydot[5] = dY1;
    ydot[6] = dY2;
    ydot[7] = dY3;
    ydot[8] = dZ0;
    ydot[9] = dZ1;
    ydot[10]= dZ2;
    ydot[11]= dZ3;
    ydot[12]= dD;
}
