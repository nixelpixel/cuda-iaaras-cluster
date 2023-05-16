//CFT.h  Complex Fourier Transform

#ifndef _IAARAS_CFT_LIBRARY_H
#define _IAARAS_CFT_LIBRARY_H

#ifndef M_PI
#	define M_PI		3.14159265358979323846	/* <math.h>'s pi constant. */
#endif

//=========================================================================================================================================

int CFT_FindPower( int N );

void CFT( int sign, double t, double *xre, double *xim, int npow );

//------------------------------------------------------------------------------------------------------------------------------------------

#endif // of #ifndef _IAARAS_CFT_LIBRARY_H
