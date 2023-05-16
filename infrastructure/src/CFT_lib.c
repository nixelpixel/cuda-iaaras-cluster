#include "CFT_lib.h"
#include <math.h>

int CFT_FindPower( int N ) {
    int i=0, j=1;
    while (j<N) { j*=2; i++; }
    return j==N ?i :-1;
}

typedef struct {
	double re;
	double im;
} complex;


void CFT(int sign, double t, double *xre, double *xim, int npow)
{
    complex cxcs,hold,xa;
    int nmax,i,nn,mm,nw,layer,ii,loc,ll,nw1,j,ij;
    int msk[32];
    double  zz,delta,w;

    for(i = 1, nmax = 1; i <= npow; i++)
        nmax*=2;

    zz = 2.*M_PI*(double)sign/(double)nmax;
    if(sign > 0)
        delta = 1./(t*(double)nmax);
    else
        delta = t;

    msk[1] = nmax/2;
    for(i = 2; i <= npow; i++)
        msk[i]=msk[i-1]/2;

    nn = nmax;
    mm = 2;
    for(layer = 1; layer <= npow; layer++)
    {
        nn = nn/2;
        nw = 0;
        for(i = 1; i<= mm; i+=2)
        {
            ii = nn*i;
            w = zz*(double)nw;
            cxcs.re = cos(w);
            cxcs.im = sin(w);
            for(j = 1; j <= nn; j++)
            {
                ii++;
                ij = ii - nn;
                xa.re = cxcs.re*xre[ii-1] - cxcs.im*xim[ii-1];
                xa.im = cxcs.re*xim[ii-1] + cxcs.im*xre[ii-1];
                xre[ii-1] = xre[ij-1] - xa.re;
                xim[ii-1] = xim[ij-1] - xa.im;
                xre[ij-1] = xre[ij-1] + xa.re;
                xim[ij-1] = xim[ij-1] + xa.im;
            }
            for(loc = 2; loc <= npow; loc++)
            {
                ll = nw - msk[loc];
                if(ll < 0)
                {
                    nw += msk[loc];
                    goto L_40;
                }
                if(ll == 0)
                {
                    nw = msk[loc+1];
                    goto L_40;
                }
                nw = ll;
            }
            L_40:;
        }
        mm *= 2;
    }
    nw = 0;
    for(i = 1; i <= nmax; i++)
    {
        nw1 = nw + 1;
        hold.re = xre[nw1-1];
        hold.im = xim[nw1-1];
        if(nw1-i < 0) goto a;
        if(nw1-i > 0)
        {
            xre[nw1-1] = xre[i-1]*delta;
            xim[nw1-1] = xim[i-1]*delta;
        }
        xre[i-1] = hold.re*delta;
        xim[i-1] = hold.im*delta;
        a:		for(loc = 1; loc <= npow; loc++)
    {
        ll = nw - msk[loc];
        if(ll < 0)
        {
            nw += msk[loc];
            goto L_80;
        }
        if(ll == 0)
        {
            nw = msk[loc+1];
            goto L_80;
        }
        nw = ll;
    }
        L_80:;
    }
}
