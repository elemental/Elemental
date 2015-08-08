/* ========================================================================= */
/* === amd_internal.h ====================================================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,                                      */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: DrTimothyAldenDavis@gmail.com                                      */
/* ------------------------------------------------------------------------- */

/* This file is for internal use in AMD itself, and does not normally need to
 * be included in user code (it is included in UMFPACK, however).   All others
 * should use amd.h instead.
 */

/* ========================================================================= */
/* === NDEBUG ============================================================== */
/* ========================================================================= */

/*
 * Turning on debugging takes some work (see below).   If you do not edit this
 * file, then debugging is always turned off, regardless of whether or not
 * -DNDEBUG is specified in your compiler options.
 *
 * If AMD is being compiled as a mexFunction, then MATLAB_MEX_FILE is defined,
 * and mxAssert is used instead of assert.  If debugging is not enabled, no
 * MATLAB include files or functions are used.  Thus, the AMD library libamd.a
 * can be safely used in either a stand-alone C program or in another
 * mexFunction, without any change.
 */

/*
    AMD will be exceedingly slow when running in debug mode.  The next three
    lines ensure that debugging is turned off.
*/
#ifndef NDEBUG
#define NDEBUG
#endif

/*
    To enable debugging, uncomment the following line:
#undef NDEBUG
*/

/* ------------------------------------------------------------------------- */
/* ANSI include files */
/* ------------------------------------------------------------------------- */

/* from stdlib.h:  size_t, malloc, free, realloc, and calloc */
#include <stdlib.h>

#if !defined(NPRINT) || !defined(NDEBUG)
/* from stdio.h:  printf.  Not included if NPRINT is defined at compile time.
 * fopen and fscanf are used when debugging. */
#include <stdio.h>
#endif

/* from limits.h:  INT_MAX and LONG_MAX */
#include <limits.h>

/* from math.h: sqrt */
#include <math.h>

/* ------------------------------------------------------------------------- */
/* MATLAB include files (only if being used in or via MATLAB) */
/* ------------------------------------------------------------------------- */

#ifdef MATLAB_MEX_FILE
#include "matrix.h"
#include "mex.h"
#endif

/* ------------------------------------------------------------------------- */
/* basic definitions */
/* ------------------------------------------------------------------------- */

#ifdef FLIP
#undef FLIP
#endif

#ifdef MAX
#undef MAX
#endif

#ifdef MIN
#undef MIN
#endif

#ifdef EMPTY
#undef EMPTY
#endif

#ifdef GLOBAL
#undef GLOBAL
#endif

#ifdef PRIVATE
#undef PRIVATE
#endif

/* FLIP is a "negation about -1", and is used to mark an integer i that is
 * normally non-negative.  FLIP (EMPTY) is EMPTY.  FLIP of a number > EMPTY
 * is negative, and FLIP of a number < EMTPY is positive.  FLIP (FLIP (i)) = i
 * for all integers i.  UNFLIP (i) is >= EMPTY. */
#define EMPTY (-1)
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) ((i < EMPTY) ? FLIP (i) : (i))

/* for integer MAX/MIN, or for doubles when we don't care how NaN's behave: */
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* logical expression of p implies q: */
#define IMPLIES(p,q) (!(p) || (q))

/* Note that the IBM RS 6000 xlc predefines TRUE and FALSE in <types.h>. */
/* The Compaq Alpha also predefines TRUE and FALSE. */
#ifdef TRUE
#undef TRUE
#endif
#ifdef FALSE
#undef FALSE
#endif

#define TRUE (1)
#define FALSE (0)
#define PRIVATE static
#define GLOBAL
#define EMPTY (-1)

/* Note that Linux's gcc 2.96 defines NULL as ((void *) 0), but other */
/* compilers (even gcc 2.95.2 on Solaris) define NULL as 0 or (0).  We */
/* need to use the ANSI standard value of 0. */
#ifdef NULL
#undef NULL
#endif

#define NULL 0

/* largest value of size_t */
#ifndef SIZE_T_MAX
#ifdef SIZE_MAX
/* C99 only */
#define SIZE_T_MAX SIZE_MAX
#else
#define SIZE_T_MAX ((size_t) (-1))
#endif
#endif

/* ------------------------------------------------------------------------- */
/* integer type for AMD: int or SuiteSparse_long */
/* ------------------------------------------------------------------------- */

#include "ElSuiteSparse/amd.h"

#if defined (DLONG) || defined (ZLONG)

#define Int ElSuiteSparse_long
#define ID  ElSuiteSparse_long_id
#define Int_MAX ElSuiteSparse_long_max

#define ElAMD_order El_amd_l_order
#define ElAMD_defaults El_amd_l_defaults
#define ElAMD_control El_amd_l_control
#define ElAMD_info El_amd_l_info
#define ElAMD_1 El_amd_l1
#define ElAMD_2 El_amd_l2
#define ElAMD_valid El_amd_l_valid
#define ElAMD_aat El_amd_l_aat
#define ElAMD_postorder El_amd_l_postorder
#define ElAMD_post_tree El_amd_l_post_tree
#define ElAMD_dump El_amd_l_dump
#define ElAMD_debug El_amd_l_debug
#define ElAMD_debug_init El_amd_l_debug_init
#define ElAMD_preprocess El_amd_l_preprocess

#else

// NOTE: This could conflict with Elemental's Int typedef, but hopefully not
//       since this is an internal header
#define Int int
#define ID "%d"
#define Int_MAX INT_MAX

#define ElAMD_order El_amd_order
#define ElAMD_defaults El_amd_defaults
#define ElAMD_control El_amd_control
#define ElAMD_info El_amd_info
#define ElAMD_1 El_amd_1
#define ElAMD_2 El_amd_2
#define ElAMD_valid El_amd_valid
#define ElAMD_aat El_amd_aat
#define ElAMD_postorder El_amd_postorder
#define ElAMD_post_tree El_amd_post_tree
#define ElAMD_dump El_amd_dump
#define ElAMD_debug El_amd_debug
#define ElAMD_debug_init El_amd_debug_init
#define ElAMD_preprocess El_amd_preprocess

#endif

/* ------------------------------------------------------------------------- */
/* AMD routine definitions (not user-callable) */
/* ------------------------------------------------------------------------- */

GLOBAL size_t ElAMD_aat
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int Len [ ],
    Int Tp [ ],
    double Info [ ]
) ;

GLOBAL void ElAMD_1
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int P [ ],
    Int Pinv [ ],
    Int Len [ ],
    Int slen,
    Int S [ ],
    double Control [ ],
    double Info [ ]
) ;

GLOBAL void ElAMD_postorder
(
    Int nn,
    Int Parent [ ],
    Int Npiv [ ],
    Int Fsize [ ],
    Int Order [ ],
    Int Child [ ],
    Int Sibling [ ],
    Int Stack [ ]
) ;

GLOBAL Int ElAMD_post_tree
(
    Int root,
    Int k,
    Int Child [ ],
    const Int Sibling [ ],
    Int Order [ ],
    Int Stack [ ]
#ifndef NDEBUG
    , Int nn
#endif
) ;

GLOBAL void ElAMD_preprocess
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int Rp [ ],
    Int Ri [ ],
    Int W [ ],
    Int Flag [ ]
) ;

/* ------------------------------------------------------------------------- */
/* debugging definitions */
/* ------------------------------------------------------------------------- */

#ifndef NDEBUG

/* from assert.h:  assert macro */
#include <assert.h>

#ifndef EXTERN
#define EXTERN extern
#endif

EXTERN Int ElAMD_debug ;

GLOBAL void ElAMD_debug_init ( char *s ) ;

GLOBAL void ElAMD_dump
(
    Int n,
    Int Pe [ ],
    Int Iw [ ],
    Int Len [ ],
    Int iwlen,
    Int pfree,
    Int Nv [ ],
    Int Next [ ],
    Int Last [ ],
    Int Head [ ],
    Int Elen [ ],
    Int Degree [ ],
    Int W [ ],
    Int nel
) ;

#ifdef ASSERT
#undef ASSERT
#endif

/* Use mxAssert if AMD is compiled into a mexFunction */
#ifdef MATLAB_MEX_FILE
#define ASSERT(expression) (mxAssert ((expression), ""))
#else
#define ASSERT(expression) (assert (expression))
#endif

#define EL_AMD_DEBUG0(params) { SUITESPARSE_PRINTF (params) ; }
#define EL_AMD_DEBUG1(params) \
 { if (ElAMD_debug >= 1) EL_SUITESPARSE_PRINTF (params) ; }
#define EL_AMD_DEBUG2(params) \
 { if (ElAMD_debug >= 2) EL_SUITESPARSE_PRINTF (params) ; }
#define EL_AMD_DEBUG3(params) \
 { if (ElAMD_debug >= 3) EL_SUITESPARSE_PRINTF (params) ; }
#define EL_AMD_DEBUG4(params) \
 { if (ElAMD_debug >= 4) EL_SUITESPARSE_PRINTF (params) ; }

#else

/* no debugging */
#define ASSERT(expression)
#define EL_AMD_DEBUG0(params)
#define EL_AMD_DEBUG1(params)
#define EL_AMD_DEBUG2(params)
#define EL_AMD_DEBUG3(params)
#define EL_AMD_DEBUG4(params)

#endif
