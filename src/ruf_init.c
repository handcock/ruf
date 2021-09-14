#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 *    Check these declarations against the C/Fortran source code.
 *    */

/* .Fortran calls */
extern void F77_NAME(matcov)(void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"matcov", (DL_FUNC) &F77_NAME(matcov), 4},
    {NULL, NULL, 0}
};

void R_init_ruf(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
