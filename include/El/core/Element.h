/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.
   Copyright (c) 2014, Jed Brown
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
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

/* Return the complex argument of a scalar
   --------------------------------------- */
EL_EXPORT ElError ElArg_s( float alpha, float* result );
EL_EXPORT ElError ElArg_d( double alpha, double* result );
EL_EXPORT ElError ElArg_c( complex_float alpha, float* result );
EL_EXPORT ElError ElArg_z( complex_double alpha, double* result );

/* Construct a complex number from its polar coordinates
   ----------------------------------------------------- */
EL_EXPORT ElError ElComplexFromPolar_c
( float r, float theta, complex_float* result );
EL_EXPORT ElError ElComplexFromPolar_z
( double r, double theta, complex_double* result );

/* Magnitude and sign
   ================== */
EL_EXPORT ElError ElAbs_i( ElInt alpha, ElInt* result );
EL_EXPORT ElError ElAbs_s( float alpha, float* result );
EL_EXPORT ElError ElAbs_d( double alpha, double* result );
EL_EXPORT ElError ElAbs_c( complex_float alpha, float* result );
EL_EXPORT ElError ElAbs_z( complex_double alpha, double* result );

EL_EXPORT ElError ElSafeAbs_c( complex_float alpha, float* result );
EL_EXPORT ElError ElSafeAbs_z( complex_double alpha, double* result );

EL_EXPORT ElError ElFastAbs_c( complex_float alpha, float* result );
EL_EXPORT ElError ElFastAbs_z( complex_double alpha, double* result );

EL_EXPORT ElError ElSgn_i( ElInt alpha, bool symmetric, ElInt* result );
EL_EXPORT ElError ElSgn_s( float alpha, bool symmetric, float* result );
EL_EXPORT ElError ElSgn_d( double alpha, bool symmetric, double* result );

/* Exponentiation
   ============== */
EL_EXPORT ElError ElExp_s( float alpha, float* result );
EL_EXPORT ElError ElExp_d( double alpha, double* result );
EL_EXPORT ElError ElExp_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElExp_z( complex_double alpha, complex_double* result );

EL_EXPORT ElError ElPow_s
( float alpha, float beta, float* result );
EL_EXPORT ElError ElPow_d
( double alpha, double beta, double* result );
EL_EXPORT ElError ElPow_c
( complex_float apha, complex_float beta, complex_float* result );
EL_EXPORT ElError ElPow_z
( complex_double alpha, complex_double beta, complex_double* result );

/* Inverse exponentiation
   ---------------------- */
EL_EXPORT ElError ElLog_s( float alpha, float* result );
EL_EXPORT ElError ElLog_d( double alpha, double* result );
EL_EXPORT ElError ElLog_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElLog_z( complex_double alpha, complex_double* result );

EL_EXPORT ElError ElSqrt_s( float alpha, float* result );
EL_EXPORT ElError ElSqrt_d( double alpha, double* result );
EL_EXPORT ElError ElSqrt_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElSqrt_z( complex_double alpha, complex_double* result );

/* Trigonometric
   ============= */
EL_EXPORT ElError ElCos_s( float alpha, float* result );
EL_EXPORT ElError ElCos_d( double alpha, double* result );
EL_EXPORT ElError ElCos_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElCos_z( complex_double alpha, complex_double* result );

EL_EXPORT ElError ElSin_s( float alpha, float* result );
EL_EXPORT ElError ElSin_d( double alpha, double* result );
EL_EXPORT ElError ElSin_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElSin_z( complex_double alpha, complex_double* result );

EL_EXPORT ElError ElTan_s( float alpha, float* result );
EL_EXPORT ElError ElTan_d( double alpha, double* result );
EL_EXPORT ElError ElTan_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElTan_z( complex_double alpha, complex_double* result );

/* Inverse trigonometric
   --------------------- */
EL_EXPORT ElError ElAcos_s( float alpha, float* result );
EL_EXPORT ElError ElAcos_d( double alpha, double* result );
EL_EXPORT ElError ElAcos_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElAcos_z( complex_double alpha, complex_double* result );

EL_EXPORT ElError ElAsin_s( float alpha, float* result );
EL_EXPORT ElError ElAsin_d( double alpha, double* result );
EL_EXPORT ElError ElAsin_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElAsin_z( complex_double alpha, complex_double* result );

EL_EXPORT ElError ElAtan_s( float alpha, float* result );
EL_EXPORT ElError ElAtan_d( double alpha, double* result );
EL_EXPORT ElError ElAtan_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElAtan_z( complex_double alpha, complex_double* result );

EL_EXPORT ElError ElAtan2_s( float y, float x, float* result );
EL_EXPORT ElError ElAtan2_d( double y, double x, double* result );

/* Hyperbolic
   ========== */
EL_EXPORT ElError ElCosh_s( float alpha, float* result );
EL_EXPORT ElError ElCosh_d( double alpha, double* result );
EL_EXPORT ElError ElCosh_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElCosh_z( complex_double alpha, complex_double* result );

EL_EXPORT ElError ElSinh_s( float alpha, float* result );
EL_EXPORT ElError ElSinh_d( double alpha, double* result );
EL_EXPORT ElError ElSinh_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElSinh_z( complex_double alpha, complex_double* result );

EL_EXPORT ElError ElTanh_s( float alpha, float* result );
EL_EXPORT ElError ElTanh_d( double alpha, double* result );
EL_EXPORT ElError ElTanh_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElTanh_z( complex_double alpha, complex_double* result );

/* Inverse hyperbolic
   ------------------ */
EL_EXPORT ElError ElAcosh_s( float alpha, float* result );
EL_EXPORT ElError ElAcosh_d( double alpha, double* result );
EL_EXPORT ElError ElAcosh_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElAcosh_z( complex_double alpha, complex_double* result );

EL_EXPORT ElError ElAsinh_s( float alpha, float* result );
EL_EXPORT ElError ElAsinh_d( double alpha, double* result );
EL_EXPORT ElError ElAsinh_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElAsinh_z( complex_double alpha, complex_double* result );

EL_EXPORT ElError ElAtanh_s( float alpha, float* result );
EL_EXPORT ElError ElAtanh_d( double alpha, double* result );
EL_EXPORT ElError ElAtanh_c( complex_float alpha, complex_float* result );
EL_EXPORT ElError ElAtanh_z( complex_double alpha, complex_double* result );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_ELEMENT_C_H */
