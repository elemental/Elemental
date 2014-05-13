/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SCALAR_CINT_H
#define EL_SCALAR_CINT_H

#ifdef __cplusplus
extern "C" {
#endif

float  ElRealPart_c( const void* alpha );
double ElRealPart_z( const void* alpha );

float  ElImagPart_c( const void* alpha );
double ElImagPart_z( const void* alpha );

void ElSetRealPart_c( void* alpha, float alphaReal );
void ElSetRealPart_z( void* alpha, double alphaReal );

void ElSetImagPart_c( void* alpha, float alphaImag );
void ElSetImagPart_z( void* alpha, double alphaImag );

void ElUpdateRealPart_c( void* alpha, float betaReal );
void ElUpdateImagPart_z( void* alpha, double betaReal );

void ElUpdateImagPart_c( void* alpha, float betaImag );
void ElUpdateImagPart_z( void* alpha, double betaImag );

float  ElAbs_s( float alpha );
double ElAbs_d( double alpha );
float  ElAbs_c( const void* alpha );
double ElAbs_z( const void* alpha );

float  ElSafeAbs_c( const void* alpha );
double ElSafeAbs_z( const void* alpha );

float  ElFastAbs_c( const void* alpha );
double ElFastAbs_z( const void* alpha );

float  ElConj_c( void* alpha );
double ElConj_z( void* alpha );

float  ElSqrt_s( float alpha );
double ElSqrt_d( double alpha );
void   ElSqrt_c( const void* alpha, void* alphaSqrt );
void   ElSqrt_z( const void* alpha, void* alphaSqrt );

float  ElCos_s( float alpha );
double ElCos_d( double alpha );
void   ElCos_c( const void* alpha, void* alphaCos );
void   ElCos_z( const void* alpha, void* alphaCos );

float  ElSin_s( float alpha );
double ElSin_d( double alpha );
void   ElSin_c( const void* alpha, void* alphaSin );
void   ElSin_z( const void* alpha, void* alphaSin );

float  ElTan_s( float alpha );
double ElTan_d( double alpha );
void   ElTan_c( const void* alpha, void* alphaTan );
void   ElTan_z( const void* alpha, void* alphaTan );

float  ElCosh_s( float alpha );
double ElCosh_d( double alpha );
void   ElCosh_c( const void* alpha, void* alphaCosh );
void   ElCosh_z( const void* alpha, void* alphaCosh );

float  ElSinh_s( float alpha );
double ElSinh_d( double alpha );
void   ElSinh_c( const void* alpha, void* alphaSinh );
void   ElSinh_z( const void* alpha, void* alphaSinh );

float  ElAcos_s( float alpha );
double ElAcos_d( double alpha );
void   ElAcos_c( const void* alpha, void* alphaAcos );
void   ElAcos_z( const void* alpha, void* alphaAcos );

float  ElAsin_s( float alpha );
double ElAsin_d( double alpha );
void   ElAsin_c( const void* alpha, void* alphaAsin );
void   ElAsin_z( const void* alpha, void* alphaAsin );

float  ElAtan_s( float alpha );
double ElAtan_d( double alpha );
void   ElAtan_c( const void* alpha, void* alphaAtan );
void   ElAtan_z( const void* alpha, void* alphaAtan );

float  ElAtan2_s( float y, float x );
double ElAtan2_d( double y, double x );

float  ElAcosh_s( float alpha );
double ElAcosh_d( double alpha );
void   ElAcosh_c( const void* alpha, void* alphaAcosh );
void   ElAcosh_z( const void* alpha, void* alphaAcosh );

float  ElAsinh_s( float alpha );
double ElAsinh_d( double alpha );
void   ElAsinh_c( const void* alpha, void* alphaAsinh );
void   ElAsinh_z( const void* alpha, void* alphaAsinh );

float  ElAtanh_s( float alpha );
double ElAtanh_d( double alpha );
void   ElAtanh_c( const void* alpha, void* alphaAtanh );
void   ElAtanh_z( const void* alpha, void* alphaAtanh );

float  ElArg_s( float alpha );
double ElArg_d( double alpha );
float  ElArg_c( const void* alpha );
double ElArg_z( const void* alpha );

void ElPolar_c( float r, float theta, void* alpha );
void ElPolar_z( double r, double theta, void* alpha );

float  ElExp_s( float alpha );
double ElExp_d( double alpha );
void   ElExp_c( const void* alpha, void* alphaExp );
void   ElExp_z( const void* alpha, void* alphaExp );

float  ElPow_s( float alpha, float beta );
double ElPow_d( double alpha, double beta );
void   ElPow_c( const void* alpha, const void* beta, void* gamma );
void   ElPow_z( const void* alpha, const void* beta, void* gamma );

float  ElLog_s( float alpha );
double ElLog_d( double alpha );
void   ElLog_c( const void* alpha, void* alphaLog );
void   ElLog_z( const void* alpha, void* alphaLog );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_SCALAR_CINT_H */
