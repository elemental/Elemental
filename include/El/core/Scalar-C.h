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

float  ElRealPart_c( const complex_float alpha );
double ElRealPart_z( const complex_double alpha );

float  ElImagPart_c( const complex_float alpha );
double ElImagPart_z( const complex_double alpha );

void ElSetRealPart_c( complex_float* alpha, float alphaReal );
void ElSetRealPart_z( complex_double* alpha, double alphaReal );

void ElSetImagPart_c( complex_float* alpha, float alphaImag );
void ElSetImagPart_z( complex_double* alpha, double alphaImag );

void ElUpdateRealPart_c( complex_float* alpha, float betaReal );
void ElUpdateImagPart_z( complex_double* alpha, double betaReal );

void ElUpdateImagPart_c( complex_float* alpha, float betaImag );
void ElUpdateImagPart_z( complex_double* alpha, double betaImag );

float  ElAbs_s( float alpha );
double ElAbs_d( double alpha );
float  ElAbs_c( complex_float alpha );
double ElAbs_z( complex_double alpha );

float  ElSafeAbs_c( complex_float alpha );
double ElSafeAbs_z( complex_double alpha );

float  ElFastAbs_c( complex_float alpha );
double ElFastAbs_z( complex_double alpha );

float  ElConj_c( complex_float alpha );
double ElConj_z( complex_double alpha );

float          ElSqrt_s( float alpha );
double         ElSqrt_d( double alpha );
complex_float  ElSqrt_c( complex_float alpha );
complex_double ElSqrt_z( complex_double alpha );

float          ElCos_s( float alpha );
double         ElCos_d( double alpha );
complex_float  ElCos_c( complex_float alpha );
complex_double ElCos_z( complex_double alpha );

float          ElSin_s( float alpha );
double         ElSin_d( double alpha );
complex_float  ElSin_c( complex_float alpha );
complex_double ElSin_z( complex_double alpha );

float          ElTan_s( float alpha );
double         ElTan_d( double alpha );
complex_float  ElTan_c( complex_float alpha );
complex_double ElTan_z( complex_double alpha );

float          ElCosh_s( float alpha );
double         ElCosh_d( double alpha );
complex_float  ElCosh_c( complex_float alpha );
complex_double ElCosh_z( complex_float alpha );

float          ElSinh_s( float alpha );
double         ElSinh_d( double alpha );
complex_float  ElSinh_c( complex_float alpha );
complex_double ElSinh_z( complex_double alpha );

float          ElAcos_s( float alpha );
double         ElAcos_d( double alpha );
complex_float  ElAcos_c( complex_float alpha );
complex_double ElAcos_z( complex_double alpha );

float          ElAsin_s( float alpha );
double         ElAsin_d( double alpha );
complex_float  ElAsin_c( complex_float alpha );
complex_double ElAsin_z( complex_double alpha );

float          ElAtan_s( float alpha );
double         ElAtan_d( double alpha );
complex_float  ElAtan_c( complex_float alpha );
complex_double ElAtan_z( complex_double alpha );

float  ElAtan2_s( float y, float x );
double ElAtan2_d( double y, double x );

float          ElAcosh_s( float alpha );
double         ElAcosh_d( double alpha );
complex_float  ElAcosh_c( complex_float alpha );
complex_double ElAcosh_z( complex_double alpha );

float          ElAsinh_s( float alpha );
double         ElAsinh_d( double alpha );
complex_float  ElAsinh_c( complex_float alpha );
complex_double ElAsinh_z( complex_double alpha );

float          ElAtanh_s( float alpha );
double         ElAtanh_d( double alpha );
complex_float  ElAtanh_c( complex_float alpha );
complex_double ElAtanh_z( complex_double alpha );

float  ElArg_s( float alpha );
double ElArg_d( double alpha );
float  ElArg_c( complex_float alpha );
double ElArg_z( complex_double alpha );

complex_float  ElPolar_c( float r, float theta );
complex_double ElPolar_z( double r, double theta );

float          ElExp_s( float alpha );
double         ElExp_d( double alpha );
complex_float  ElExp_c( complex_float alpha );
complex_double ElExp_z( complex_double alpha );

float          ElPow_s( float alpha, float beta );
double         ElPow_d( double alpha, double beta );
complex_float  ElPow_c( complex_float apha, complex_float beta );
complex_double ElPow_z( complex_double alpha, complex_double beta );

float          ElLog_s( float alpha );
double         ElLog_d( double alpha );
complex_double ElLog_c( complex_float alpha );
complex_double ElLog_z( complex_double alpha );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_SCALAR_C_H */
