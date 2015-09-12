/* ========================================================================= */
/* === AMD_defaults ======================================================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: DrTimothyAldenDavis@gmail.com                                      */
/* ------------------------------------------------------------------------- */

/* User-callable.  Sets default control parameters for AMD.  See amd.h
 * for details.
 */

#include "ElSuiteSparse/amd_internal.h"

/* ========================================================================= */
/* === AMD defaults ======================================================== */
/* ========================================================================= */

GLOBAL void ElAMD_defaults
(
    double Control [ ]
)
{
    Int i ;

    if (Control != (double *) NULL)
    {
	for (i = 0 ; i < EL_AMD_CONTROL ; i++)
	{
	    Control [i] = 0 ;
	}
	Control [EL_AMD_DENSE] = EL_AMD_DEFAULT_DENSE ;
	Control [EL_AMD_AGGRESSIVE] = EL_AMD_DEFAULT_AGGRESSIVE ;
    }
}
