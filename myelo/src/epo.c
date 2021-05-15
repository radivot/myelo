#include <R.h>
static double parms[10];
#define Bmax  parms[0]
#define Epo0  parms[1]
#define Kd    parms[2]
#define kde   parms[3]
#define kdi   parms[4]
#define ke    parms[5]
#define kex   parms[6]
#define kon   parms[7]
#define kt    parms[8]
#define scale parms[9]

void parmsEpo(void (* odeparms)(int *, double *))
{   int N=10;
    odeparms(&N, parms); 
}


void derivsEpo(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
	double v1,v2,v3,v4,v5,v6,v7,v8;
	double Epo,EpoR,EpoEpoR,EpoEpoRi,dEpoi,dEpoe;
    Epo=y[0];
    EpoR=y[1];
    EpoEpoR=y[2];
    EpoEpoRi=y[3];
    dEpoi=y[4];
    dEpoe=y[5];
    v1=kon*Epo*EpoR;
    v2=kon*Kd*EpoEpoR;
    v3=kt*Bmax;
    v4=kt*EpoR;
    v5=ke*EpoEpoR;
    v6=kex*EpoEpoRi;
    v7=kdi*EpoEpoRi;
    v8=kde*EpoEpoRi;
    ydot[0] = -v1+v2+v6;
    ydot[1] = -v1+v2+v3-v4+v6;
    ydot[2] = v1-v2-v5;
    ydot[3] = v5-v6-v7-v8;
    ydot[4] = v7;
    ydot[5] = v8;
    yout[0] = scale*(Epo+dEpoe);
    yout[1] = scale*EpoEpoR;
    yout[2] = scale*(EpoEpoRi+dEpoi);
}

