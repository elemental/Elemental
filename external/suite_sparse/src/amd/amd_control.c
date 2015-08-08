/* ========================================================================= */
/* === AMD_control ========================================================= */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: DrTimothyAldenDavis@gmail.com                                      */
/* ------------------------------------------------------------------------- */

/* User-callable.  Prints the control parameters for AMD.  See amd.h
 * for details.  If the Control array is not present, the defaults are
 * printed instead.
 */

#include "ElSuiteSparse/amd_internal.h"

GLOBAL void ElAMD_control
(
    double Control [ ]
)
{
    double alpha ;
    Int aggressive ;

    if (Control != (double *) NULL)
    {
	alpha = Control [EL_AMD_DENSE] ;
	aggressive = Control [EL_AMD_AGGRESSIVE] != 0 ;
    }
    else
    {
	alpha = EL_AMD_DEFAULT_DENSE ;
	aggressive = EL_AMD_DEFAULT_AGGRESSIVE ;
    }

    EL_SUITESPARSE_PRINTF ((
        "\nAMD version %d.%d.%d, %s: approximate minimum degree ordering\n"
	"    dense row parameter: %g\n", EL_AMD_MAIN_VERSION, 
        EL_AMD_SUB_VERSION, EL_AMD_SUBSUB_VERSION, EL_AMD_DATE, alpha)) ;

    if (alpha < 0)
    {
	EL_SUITESPARSE_PRINTF (("    no rows treated as dense\n")) ;
    }
    else
    {
	EL_SUITESPARSE_PRINTF ((
	"    (rows with more than max (%g * sqrt (n), 16) entries are\n"
	"    considered \"dense\", and placed last in output permutation)\n",
	alpha)) ;
    }

    if (aggressive)
    {
	EL_SUITESPARSE_PRINTF (("    aggressive absorption:  yes\n")) ;
    }
    else
    {
	EL_SUITESPARSE_PRINTF (("    aggressive absorption:  no\n")) ;
    }

    EL_SUITESPARSE_PRINTF (("    size of AMD integer: %d\n\n", sizeof (Int))) ;
}
