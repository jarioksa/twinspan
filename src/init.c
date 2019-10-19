#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> /* for NULL */
#include <R_ext/Rdynload.h>

/* extern definitions */

extern void data2hill(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pseudo)(void *, void *,void *, void *, void *, void *,
    void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(class)(void *, void *, void *, void *, void *, void *,
			    void *, void *, void *, void *, void *, void *,
			    void *, void *, void *, void *, void *, void *,
			    void *, void *, void *, void *, void *, void *,
			    void *, void *, void *, void *, void *, void *,
			    void *, void *, void *, void *, void *, void *,
			    void *, void *, void *, void *);
extern void F77_NAME(revspec)(void *, void *, void *, void *, void *, void *,
			      void *, void *, void *);

/* C */

static const R_CMethodDef CCalls[] = {
    {"data2hill", (DL_FUNC) &data2hill, 6},
    {NULL, NULL, 0}
};

/* FORTRAN */

static const R_FortranMethodDef FortranCalls[] = {
    {"pseudo", (DL_FUNC) &F77_NAME(pseudo), 13},
    {"class", (DL_FUNC) &F77_NAME(class), 40},
    {"revspec", (DL_FUNC) &F77_NAME(revspec), 9},
    {NULL, NULL, 0}
};

/* registration */

void R_init_twinspan(DllInfo *dll)
{
    R_registerRoutines(dll, CCalls, NULL, FortranCalls, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
