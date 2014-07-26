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
#ifndef EL_SCALAR_C_H
#define EL_SCALAR_C_H

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

/* Basic complex entry manipulation
   ================================ */
ElError ElRealPart_c( complex_float alpha, float* alphaReal );
ElError ElRealPart_z( complex_double alpha, double* alphaReal );

ElError ElImagPart_c( complex_float alpha, float* alphaImag );
ElError ElImagPart_z( complex_double alpha, double* alphaImag );

ElError ElSetRealPart_c( complex_float* alpha, float alphaReal );
ElError ElSetRealPart_z( complex_double* alpha, double alphaReal );

ElError ElSetImagPart_c( complex_float* alpha, float alphaImag );
ElError ElSetImagPart_z( complex_double* alpha, double alphaImag );

ElError ElUpdateRealPart_c( complex_float* alpha, float betaReal );
ElError ElUpdateImagPart_z( complex_double* alpha, double betaReal );

ElError ElUpdateImagPart_c( complex_float* alpha, float betaImag );
ElError ElUpdateImagPart_z( complex_double* alpha, double betaImag );

ElError ElConj_c( complex_float alpha, complex_float* alphaConj );
ElError ElConj_z( complex_double alpha, complex_double* alphaConj );

ElError ElArg_s( float alpha, float* alphaArg );
ElError ElArg_d( double alpha, double* alphaArg );
ElError ElArg_c( complex_float alpha, float* alphaArg );
ElError ElArg_z( complex_double alpha, double* alphaArg );

/* NOTE: This is not simply named 'Polar' due to a conflict with the
         polar decomposition function */
ElError ElComplexFromPolar_c( float r, float theta, complex_float* alpha );
ElError ElComplexFromPolar_z( double r, double theta, complex_double* alpha );

/* Size measurements
   ================= */
ElError ElAbs_i( ElInt alpha, ElInt* alphaAbs );
ElError ElAbs_s( float alpha, float* alphaAbs );
ElError ElAbs_d( double alpha, double* alphaAbs );
ElError ElAbs_c( complex_float alpha, float* alphaAbs );
ElError ElAbs_z( complex_double alpha, double* alphaAbs );

ElError ElSafeAbs_c( complex_float alpha, float* alphaAbs );
ElError ElSafeAbs_z( complex_double alpha, double* alphaAbs );

ElError ElFastAbs_c( complex_float alpha, float* alphaAbs );
ElError ElFastAbs_z( complex_double alpha, double* alphaAbs );

/* Exponentiation
   ============== */
ElError ElExp_s( float alpha, float* alphaExp );
ElError ElExp_d( double alpha, double* alphaExp );
ElError ElExp_c( complex_float alpha, complex_float* alphaExp );
ElError ElExp_z( complex_double alpha, complex_double* alphaExp );

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
ElError ElLog_s( float alpha, float* alphaLog );
ElError ElLog_d( double alpha, double* alphaLog );
ElError ElLog_c( complex_float alpha, complex_float* alphaLog );
ElError ElLog_z( complex_double alpha, complex_double* alphaLog );

ElError ElSqrt_s( float alpha, float* alphaSqrt );
ElError ElSqrt_d( double alpha, double* alphaSqrt );
ElError ElSqrt_c( complex_float alpha, complex_float* alphaSqrt );
ElError ElSqrt_z( complex_double alpha, complex_double* alphaSqrt );

/* Trigonometric
   ============= */
ElError ElCos_s( float alpha, float* alphaCos );
ElError ElCos_d( double alpha, double* alphaCos );
ElError ElCos_c( complex_float alpha, complex_float* alphaCos );
ElError ElCos_z( complex_double alpha, complex_double* alphaCos );

ElError ElSin_s( float alpha, float* alphaSin );
ElError ElSin_d( double alpha, double* alphaSin );
ElError ElSin_c( complex_float alpha, complex_float* alphaSin );
ElError ElSin_z( complex_double alpha, complex_double* alphaSin );

ElError ElTan_s( float alpha, float* alphaTan );
ElError ElTan_d( double alpha, double* alphaTan );
ElError ElTan_c( complex_float alpha, complex_float* alphaTan );
ElError ElTan_z( complex_double alpha, complex_double* alphaTan );

/* Inverse trigonometric
   --------------------- */
ElError ElAcos_s( float alpha, float* alphaAcos );
ElError ElAcos_d( double alpha, double* alphaAcos );
ElError ElAcos_c( complex_float alpha, complex_float* alphaAcos );
ElError ElAcos_z( complex_double alpha, complex_double* alphaAcos );

ElError ElAsin_s( float alpha, float* alphaAsin );
ElError ElAsin_d( double alpha, double* alphaAsin );
ElError ElAsin_c( complex_float alpha, complex_float* alphaAsin );
ElError ElAsin_z( complex_double alpha, complex_double* alphaAsin );

ElError ElAtan_s( float alpha, float* alphaAtan );
ElError ElAtan_d( double alpha, double* alphaAtan );
ElError ElAtan_c( complex_float alpha, complex_float* alphaAtan );
ElError ElAtan_z( complex_double alpha, complex_double* alphaAtan );

ElError ElAtan2_s( float y, float x, float* result );
ElError ElAtan2_d( double y, double x, double* result );

/* Hyperbolic
   ========== */
ElError ElCosh_s( float alpha, float* alphaCosh );
ElError ElCosh_d( double alpha, double* alphaCosh );
ElError ElCosh_c( complex_float alpha, complex_float* alphaCosh );
ElError ElCosh_z( complex_double alpha, complex_double* alphaCosh );

ElError ElSinh_s( float alpha, float* alphaSinh );
ElError ElSinh_d( double alpha, double* alphaSinh );
ElError ElSinh_c( complex_float alpha, complex_float* alphaSinh );
ElError ElSinh_z( complex_double alpha, complex_double* alphaSinh );

ElError ElTanh_s( float alpha, float* alphaTanh );
ElError ElTanh_d( double alpha, double* alphaTanh );
ElError ElTanh_c( complex_float alpha, complex_float* alphaTanh );
ElError ElTanh_z( complex_double alpha, complex_double* alphaTanh );

/* Inverse hyperbolic
   ------------------ */
ElError ElAcosh_s( float alpha, float* alphaAcosh );
ElError ElAcosh_d( double alpha, double* alphaAcosh );
ElError ElAcosh_c( complex_float alpha, complex_float* alphaAcosh );
ElError ElAcosh_z( complex_double alpha, complex_double* alphaAcosh );

ElError ElAsinh_s( float alpha, float* alphaAsinh );
ElError ElAsinh_d( double alpha, double* alphaAsinh );
ElError ElAsinh_c( complex_float alpha, complex_float* alphaAsinh );
ElError ElAsinh_z( complex_double alpha, complex_double* alphaAsinh );

ElError ElAtanh_s( float alpha, float* alphaAtanh );
ElError ElAtanh_d( double alpha, double* alphaAtanh );
ElError ElAtanh_c( complex_float alpha, complex_float* alphaAtanh );
ElError ElAtanh_z( complex_double alpha, complex_double* alphaAtanh );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_SCALAR_C_H */
