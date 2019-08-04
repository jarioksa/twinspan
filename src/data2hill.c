/*********************************************************
 Changes a data matrix into a form used internally in Mark Hill's
 TWINSPAN. This is based on similar function in the vegan package used
 there to translate rectangular R data into internal form used by
 decorana. The internal form is different in TWINSPAN: data are in an
 integer vector in pairs where the first item gives the index of
 species, and the second its abundance, and the end of the SU is
 indicated with "species number" -1. To handle decimal data,
 abundances are multiplied with 1000 and rounded to integers. For NZ
 non-zero items and NR rows, the length of this vector is 2*NZ + NR.
**********************************************************/

#include <R.h>

void data2hill(double *x,
	       int *mi, int *n, int *nid, int *ibegin, int *idat)
{
    int nr, nc, i, j, ij, now;
    
    nr = *mi;
    nc = *n;
    if (nr <= 0 || nc <= 0)
	error("zero extent dimensions");
     
    for (i=0, now=0; i<nr; i++) {
	ibegin[i] = now + 1;
	for (j=0; j<nc; j++) {
	    ij = i + nr*j;
	    if (x[ij] > 0.0) {
		idat[now++] = j+1;
		idat[now++] = (int) (1000.0 * x[ij] + 0.5);
	    }
	}
	idat[now++] = -1;
    }
}
