#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void fast_inner_prod(void *, void *, void *, void *, void *);
extern void multi_across_fip(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"fast_inner_prod",  (DL_FUNC) &fast_inner_prod,   5},
    {"multi_across_fip", (DL_FUNC) &multi_across_fip, 15},
    {NULL, NULL, 0}
};

void R_init_wbsts(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
