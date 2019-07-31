/*********************************************************
 Changes a data matrix into a form used internally in Mark Hill's
 programs DECORANA and TWINSPAN.  This is a part of the function in
 the vegan package used there to translate rectangular R data into
 internal form used by decorana.
**********************************************************/

#include <R.h>

void data2hill(double *x,
	       int *mi, int *n, int *nid, int *ibegin, int *iend,
	       int *idat, double *qidat)
{
    int nr, nc, i, j, ij, now;
    
    nr = *mi;
    nc = *n;
    if (nr <= 0 || nc <= 0)
	error("zero extent dimensions");
     
    now=0;
    for (i=0; i<nr; i++) {
	for (j=0; j<nc; j++) {
	    ij = i+nr*j;
	    if (x[ij] > 0.0) {
		idat[now] = j+1;
		qidat[now] = x[ij];
		now++;
	    }
	}
	iend[i] = now;
    }
    ibegin[0] = 1;
    for (i=1; i<nr; i++)
	ibegin[i] = iend[i-1] + 1;
    *mi = nr;
    *n = nc;
    *nid = now;
}
