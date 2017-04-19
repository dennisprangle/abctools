#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* Symbol registration table */
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

void nnone(double *x, int *lx, int *k, double *D);
void nnk(double *x, int *lx, int *cols, int *k, double *D);

static const R_CMethodDef R_CDef[] = {
   CALLDEF(nnone, 4),
   CALLDEF(nnk, 5),
   {NULL, NULL, 0}
};

void R_init_abctools(DllInfo *dll)
{
    R_registerRoutines(dll, R_CDef, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
