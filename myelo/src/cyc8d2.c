/* file cyc8d2.c (with maturation out of G1)*/
#include <R.h>
static double parms[9];
#define kS parms[0]
#define k2 parms[1]
#define kM parms[2]
#define fSurv parms[3]
#define k1 parms[4]
#define kMat parms[5]
#define kLeaveMarrow parms[6]
#define kdG1mat parms[7]
#define kdSmat parms[8]

/* initializer  */
void parmsCyc8d2(void (* odeparms)(int *, double *))
{
    int N=9;
    odeparms(&N, parms);
}

/* Derivatives */
void derivsCyc8d2 (int *neq, double *t, double *y, double *ydot)
{   double M,G1,Mat,Sa,Sb,Sc,Sd,G2,drug;
    M=y[0];
    G1=y[1];Mat=y[2];
    Sa=y[3];Sb=y[4];Sc=y[5];Sd=y[6];
    G2=y[7];drug=y[8];
    ydot[0] = k2*G2-kM*M;
    if (drug == 0) {
        ydot[1] = 2*kM*fSurv*M-(k1+kMat)*G1;
        ydot[2] = kMat*G1-kLeaveMarrow*Mat;
        ydot[3] = k1*G1-kS*Sa;
        ydot[4] = kS*Sa-kS*Sb;
        ydot[5] = kS*Sb-kS*Sc;
        ydot[6] = kS*Sc-kS*Sd;
        ydot[7] = kS*Sd-k2*G2;
    }
    if (drug == 1) {
        ydot[1] = 2*kM*fSurv*M-(k1+kMat)*G1;
        ydot[2] = kMat*G1+2*fSurv*kdSmat*(Sa+Sb+Sc+Sd)-kLeaveMarrow*Mat;
        ydot[3] = k1*G1-(kS+kdSmat)*Sa;
        ydot[4] = kS*Sa-(kS+kdSmat)*Sb;
        ydot[5] = kS*Sb-(kS+kdSmat)*Sc;
        ydot[6] = kS*Sc-(kS+kdSmat)*Sd;
        ydot[7] = kS*Sd-k2*G2;
    }
    if (drug == 2) {
        ydot[1] = 2*kM*fSurv*M-(k1+kMat+kdG1mat)*G1;
        ydot[2] = kMat*G1+2*fSurv*kdG1mat*G1-kLeaveMarrow*Mat;
        ydot[3] = k1*G1-kS*Sa;
        ydot[4] = kS*Sa-kS*Sb;
        ydot[5] = kS*Sb-kS*Sc;
        ydot[6] = kS*Sc-kS*Sd;
        ydot[7] = kS*Sd-k2*G2;
    }
    ydot[8]= 0;
}

