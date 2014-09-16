/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.
   Copyright (c) 2014, Jed Brown
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_ELEMENT_C_H
#define EL_ELEMENT_C_H

#ifdef __cplusplus
extern "C" {
typedef struct { float real, imag; } complex_float;
typedef struct { double real, imag; } complex_double;
#else
#include <complex.h>
typedef complex float complex_float;
typedef complex double complex_double;
#endif

/* TODO: Better interface between C and C++ complex */

/* Basic element manipulation
   ========================== */

/* Return the real/imaginary part of a complex number
   -------------------------------------------------- */
ElError ElRealPart_c( complex_float alpha, float* result );
ElError ElRealPart_z( complex_double alpha, double* result );

ElError ElImagPart_c( complex_float alpha, float* result );
ElError ElImagPart_z( complex_double alpha, double* result );

/* Set the real/imaginary part of a complex number
   ----------------------------------------------- */
ElError ElSetRealPart_c( complex_float* alpha, float alphaReal );
ElError ElSetRealPart_z( complex_double* alpha, double alphaReal );

ElError ElSetImagPart_c( complex_float* alpha, float alphaImag );
ElError ElSetImagPart_z( complex_double* alpha, double alphaImag );

/* Update the real/imaginary part of a complex number
   -------------------------------------------------- */
ElError ElUpdateRealPart_c( complex_float* alpha, float betaReal );
ElError ElUpdateImagPart_z( complex_double* alpha, double betaReal );

ElError ElUpdateImagPart_c( complex_float* alpha, float betaImag );
ElError ElUpdateImagPart_z( complex_double* alpha, double betaImag );

/* Conjugate a complex number
   -------------------------- */
ElError ElConj_c( complex_float alpha, complex_float* result );
ElError ElConj_z( complex_double alpha, complex_double* result );

/* Return the complex argument of a scalar
   --------------------------------------- */
ElError ElArg_s( float alpha, float* result );
ElError ElArg_d( double alpha, double* result );
ElError ElArg_c( complex_float alpha, float* result );
ElError ElArg_z( complex_double alpha, double* result );

/* Construct a complex number from its polar coordinates
   ----------------------------------------------------- */
ElError ElComplexFromPolar_c( float r, float theta, complex_float* result );
ElError ElComplexFromPolar_z( double r, double theta, complex_double* result );

/* Magnitude and sign
   ================== */
ElError ElAbs_i( ElInt alpha, ElInt* result );
ElError ElAbs_s( float alpha, float* result );
ElError ElAbs_d( double alpha, double* result );
ElError ElAbs_c( complex_float alpha, float* result );
ElError ElAbs_z( complex_double alpha, double* result );

ElError ElSafeAbs_c( complex_float alpha, float* result );
ElError ElSafeAbs_z( complex_double alpha, double* result );

ElError ElFastAbs_c( complex_float alpha, float* result );
ElError ElFastAbs_z( complex_double alpha, double* result );

ElError ElSgn_i( ElInt alpha, bool symmetric, ElInt* result );
ElError ElSgn_s( float alpha, bool symmetric, float* result );
ElError ElSgn_d( double alpha, bool symmetric, double* result );

/* Exponentiation
   ============== */
ElError ElExp_s( float alpha, float* result );
ElError ElExp_d( double alpha, double* result );
ElError ElExp_c( complex_float alpha, complex_float* result );
ElError ElExp_z( complex_double alpha, complex_double* result );

ElError ElPow_s
( float alpha, float beta, float* result );
ElError ElPow_d
( double alpha, double beta, double* result );
ElError ElPow_c
( complex_float apha, complex_float beta, complex_float* result );
ElError ElPow_z
( complex_double alpha, complex_double beta, complex_double* result );

/* Inverse exponentiation
   ---------------------- */
ElError ElLog_s( float alpha, float* result );
ElError ElLog_d( double alpha, double* result );
ElError ElLog_c( complex_float alpha, complex_float* result );
ElError ElLog_z( complex_double alpha, complex_double* result );

ElError ElSqrt_s( float alpha, float* result );
ElError ElSqrt_d( double alpha, double* result );
ElError ElSqrt_c( complex_float alpha, complex_float* result );
ElError ElSqrt_z( complex_double alpha, complex_double* result );

/* Trigonometric
   ============= */
ElError ElCos_s( float alpha, float* result );
ElError ElCos_d( double alpha, double* result );
ElError ElCos_c( complex_float alpha, complex_float* result );
ElError ElCos_z( complex_double alpha, complex_double* result );

ElError ElSin_s( float alpha, float* result );
ElError ElSin_d( double alpha, double* result );
ElError ElSin_c( complex_float alpha, complex_float* result );
ElError ElSin_z( complex_double alpha, complex_double* result );

ElError ElTan_s( float alpha, float* result );
ElError ElTan_d( double alpha, double* result );
ElError ElTan_c( complex_float alpha, complex_float* result );
ElError ElTan_z( complex_double alpha, complex_double* result );

/* Inverse trigonometric
   --------------------- */
ElError ElAcos_s( float alpha, float* result );
ElError ElAcos_d( double alpha, double* result );
ElError ElAcos_c( complex_float alpha, complex_float* result );
ElError ElAcos_z( complex_double alpha, complex_double* result );

ElError ElAsin_s( float alpha, float* result );
ElError ElAsin_d( double alpha, double* result );
ElError ElAsin_c( complex_float alpha, complex_float* result );
ElError ElAsin_z( complex_double alpha, complex_double* result );

ElError ElAtan_s( float alpha, float* result );
ElError ElAtan_d( double alpha, double* result );
ElError ElAtan_c( complex_float alpha, complex_float* result );
ElError ElAtan_z( complex_double alpha, complex_double* result );

ElError ElAtan2_s( float y, float x, float* result );
ElError ElAtan2_d( double y, double x, double* result );

/* Hyperbolic
   ========== */
ElError ElCosh_s( float alpha, float* result );
ElError ElCosh_d( double alpha, double* result );
ElError ElCosh_c( complex_float alpha, complex_float* result );
ElError ElCosh_z( complex_double alpha, complex_double* result );

ElError ElSinh_s( float alpha, float* result );
ElError ElSinh_d( double alpha, double* result );
ElError ElSinh_c( complex_float alpha, complex_float* result );
ElError ElSinh_z( complex_double alpha, complex_double* result );

ElError ElTanh_s( float alpha, float* result );
ElError ElTanh_d( double alpha, double* result );
ElError ElTanh_c( complex_float alpha, complex_float* result );
ElError ElTanh_z( complex_double alpha, complex_double* result );

/* Inverse hyperbolic
   ------------------ */
ElError ElAcosh_s( float alpha, float* result );
ElError ElAcosh_d( double alpha, double* result );
ElError ElAcosh_c( complex_float alpha, complex_float* result );
ElError ElAcosh_z( complex_double alpha, complex_double* result );

ElError ElAsinh_s( float alpha, float* result );
ElError ElAsinh_d( double alpha, double* result );
ElError ElAsinh_c( complex_float alpha, complex_float* result );
ElError ElAsinh_z( complex_double alpha, complex_double* result );

ElError ElAtanh_s( float alpha, float* result );
ElError ElAtanh_d( double alpha, double* result );
ElError ElAtanh_c( complex_float alpha, complex_float* result );
ElError ElAtanh_z( complex_double alpha, complex_double* result );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_ELEMENT_C_H */
