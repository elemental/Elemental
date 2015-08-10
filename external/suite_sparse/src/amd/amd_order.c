/* ========================================================================= */
/* === AMD_order =========================================================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: DrTimothyAldenDavis@gmail.com                                      */
/* ------------------------------------------------------------------------- */

/* User-callable AMD minimum degree ordering routine.  See amd.h for
 * documentation.
 */

#include "ElSuiteSparse/amd_internal.h"

/* ========================================================================= */
/* === AMD_order =========================================================== */
/* ========================================================================= */

GLOBAL Int ElAMD_order
(
  Int n,
  const Int Ap [ ],
  const Int Ai [ ],
  Int P [ ],
  double Control [ ],
  double Info [ ]
)
{
    Int *Len, *S, nz, i, *Pinv, info, status, *Rp, *Ri, *Cp, *Ci, ok ;
    size_t nzaat, slen ;
    double mem = 0 ;

#ifndef NDEBUG
    ElAMD_debug_init ("amd") ;
#endif

    /* clear the Info array, if it exists */
    info = Info != (double *) NULL ;
    if (info)
    {
	for (i = 0 ; i < EL_AMD_INFO ; i++)
	{
	    Info [i] = EMPTY ;
	}
	Info [EL_AMD_N] = n ;
	Info [EL_AMD_STATUS] = EL_AMD_OK ;
    }

    if (n == 0)
    {
	return EL_AMD_OK;	    /* n is 0 so there's nothing to do */
    }

    /* make sure inputs exist and n is >= 0 */
    if (Ai == (Int *) NULL || Ap == (Int *) NULL || P == (Int *) NULL || n < 0)
    {
	if (info) Info [EL_AMD_STATUS] = EL_AMD_INVALID ;
	return (EL_AMD_INVALID) ;	    /* arguments are invalid */
    }

    nz = Ap [n] ;
    if (info)
    {
	Info [EL_AMD_NZ] = nz ;
    }
    if (nz < 0)
    {
	if (info) Info [EL_AMD_STATUS] = EL_AMD_INVALID ;
	return (EL_AMD_INVALID) ;
    }

    /* check if n or nz will cause size_t overflow */
    if (((size_t) n) >= SIZE_T_MAX / sizeof (Int)
     || ((size_t) nz) >= SIZE_T_MAX / sizeof (Int))
    {
	if (info) Info [EL_AMD_STATUS] = EL_AMD_OUT_OF_MEMORY ;
	return (EL_AMD_OUT_OF_MEMORY) ;	    /* problem too large */
    }

    /* check the input matrix:	AMD_OK, AMD_INVALID, or AMD_OK_BUT_JUMBLED */
    status = ElAMD_valid (n, n, Ap, Ai) ;

    if (status == EL_AMD_INVALID)
    {
	if (info) Info [EL_AMD_STATUS] = EL_AMD_INVALID ;
	return (EL_AMD_INVALID) ;	    /* matrix is invalid */
    }

    /* allocate two size-n integer workspaces */
    Len  = ElSuiteSparse_malloc (n, sizeof (Int)) ;
    Pinv = ElSuiteSparse_malloc (n, sizeof (Int)) ;
    mem += n ;
    mem += n ;
    if (!Len || !Pinv)
    {
	/* :: out of memory :: */
	ElSuiteSparse_free (Len) ;
	ElSuiteSparse_free (Pinv) ;
	if (info) Info [EL_AMD_STATUS] = EL_AMD_OUT_OF_MEMORY ;
	return (EL_AMD_OUT_OF_MEMORY) ;
    }

    if (status == EL_AMD_OK_BUT_JUMBLED)
    {
	/* sort the input matrix and remove duplicate entries */
	EL_AMD_DEBUG1 (("Matrix is jumbled\n")) ;
	Rp = ElSuiteSparse_malloc (n+1, sizeof (Int)) ;
	Ri = ElSuiteSparse_malloc (nz,  sizeof (Int)) ;
	mem += (n+1) ;
	mem += MAX (nz,1) ;
	if (!Rp || !Ri)
	{
	    /* :: out of memory :: */
	    ElSuiteSparse_free (Rp) ;
	    ElSuiteSparse_free (Ri) ;
	    ElSuiteSparse_free (Len) ;
	    ElSuiteSparse_free (Pinv) ;
	    if (info) Info [EL_AMD_STATUS] = EL_AMD_OUT_OF_MEMORY ;
	    return (EL_AMD_OUT_OF_MEMORY) ;
	}
	/* use Len and Pinv as workspace to create R = A' */
	ElAMD_preprocess (n, Ap, Ai, Rp, Ri, Len, Pinv) ;
	Cp = Rp ;
	Ci = Ri ;
    }
    else
    {
	/* order the input matrix as-is.  No need to compute R = A' first */
	Rp = NULL ;
	Ri = NULL ;
	Cp = (Int *) Ap ;
	Ci = (Int *) Ai ;
    }

    /* --------------------------------------------------------------------- */
    /* determine the symmetry and count off-diagonal nonzeros in A+A' */
    /* --------------------------------------------------------------------- */

    nzaat = ElAMD_aat (n, Cp, Ci, Len, P, Info) ;
    EL_AMD_DEBUG1 (("nzaat: %g\n", (double) nzaat)) ;
    ASSERT ((MAX (nz-n, 0) <= nzaat) && (nzaat <= 2 * (size_t) nz)) ;

    /* --------------------------------------------------------------------- */
    /* allocate workspace for matrix, elbow room, and 6 size-n vectors */
    /* --------------------------------------------------------------------- */

    S = NULL ;
    slen = nzaat ;			/* space for matrix */
    ok = ((slen + nzaat/5) >= slen) ;	/* check for size_t overflow */
    slen += nzaat/5 ;			/* add elbow room */
    for (i = 0 ; ok && i < 7 ; i++)
    {
	ok = ((slen + n) > slen) ;	/* check for size_t overflow */
	slen += n ;			/* size-n elbow room, 6 size-n work */
    }
    mem += slen ;
    ok = ok && (slen < SIZE_T_MAX / sizeof (Int)) ; /* check for overflow */
    ok = ok && (slen < Int_MAX) ;	/* S[i] for Int i must be OK */
    if (ok)
    {
	S = ElSuiteSparse_malloc (slen, sizeof (Int)) ;
    }
    EL_AMD_DEBUG1 (("slen %g\n", (double) slen)) ;
    if (!S)
    {
	/* :: out of memory :: (or problem too large) */
	ElSuiteSparse_free (Rp) ;
	ElSuiteSparse_free (Ri) ;
	ElSuiteSparse_free (Len) ;
	ElSuiteSparse_free (Pinv) ;
	if (info) Info [EL_AMD_STATUS] = EL_AMD_OUT_OF_MEMORY ;
	return (EL_AMD_OUT_OF_MEMORY) ;
    }
    if (info)
    {
	/* memory usage, in bytes. */
	Info [EL_AMD_MEMORY] = mem * sizeof (Int) ;
    }

    /* --------------------------------------------------------------------- */
    /* order the matrix */
    /* --------------------------------------------------------------------- */

    ElAMD_1 (n, Cp, Ci, P, Pinv, Len, slen, S, Control, Info) ;

    /* --------------------------------------------------------------------- */
    /* free the workspace */
    /* --------------------------------------------------------------------- */

    ElSuiteSparse_free (Rp) ;
    ElSuiteSparse_free (Ri) ;
    ElSuiteSparse_free (Len) ;
    ElSuiteSparse_free (Pinv) ;
    ElSuiteSparse_free (S) ;
    if (info) Info [EL_AMD_STATUS] = status ;
    return (status) ;	    /* successful ordering */
}
