#include <R.h>
static double parms[7];
#define kde   parms[0]
#define kdi   parms[1]
#define ke    parms[2]
#define kex   parms[3]
#define kon   parms[4]
#define kt    parms[5]
#define scale parms[6]

void parmscHC(void (* odeparms)(int *, double *))
{   int N=7;
    odeparms(&N, parms);
}


void derivscHC(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
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
    v2=kon*164*EpoEpoR;
    v3=kt*516;
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

