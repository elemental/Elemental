/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISTMATRIX_C_H
#define EL_DISTMATRIX_C_H

#ifdef __cplusplus
extern "C" {
#endif

struct ElDistMatrix_s; typedef struct ElDistMatrix_s ElDistMatrix_s;
struct ElDistMatrix_d; typedef struct ElDistMatrix_d ElDistMatrix_d;
struct ElDistMatrix_c; typedef struct ElDistMatrix_c ElDistMatrix_c;
struct ElDistMatrix_z; typedef struct ElDistMatrix_z ElDistMatrix_z;

ElDistMatrix_s* ElDistMatrixCreate_s( const ElGrid* g );
ElDistMatrix_d* ElDistMatrixCreate_d( const ElGrid* g );
ElDistMatrix_c* ElDistMatrixCreate_c( const ElGrid* g );
ElDistMatrix_z* ElDistMatrixCreate_z( const ElGrid* g );

void ElDistMatrixDestroy_s( const ElDistMatrix_s* A );
void ElDistMatrixDestroy_d( const ElDistMatrix_d* A );
void ElDistMatrixDestroy_c( const ElDistMatrix_c* A );
void ElDistMatrixDestroy_z( const ElDistMatrix_z* A );

void ElDistMatrixEmpty_s( ElDistMatrix_s* A );
void ElDistMatrixEmpty_d( ElDistMatrix_d* A );
void ElDistMatrixEmpty_c( ElDistMatrix_c* A );
void ElDistMatrixEmpty_z( ElDistMatrix_z* A );

void ElDistMatrixResize_s( ElDistMatrix_s* A, ElInt height, ElInt width );
void ElDistMatrixResize_d( ElDistMatrix_d* A, ElInt height, ElInt width );
void ElDistMatrixResize_c( ElDistMatrix_c* A, ElInt height, ElInt width );
void ElDistMatrixResize_z( ElDistMatrix_z* A, ElInt height, ElInt width );

void ElDistMatrixResizeWithLDim_s
( ElDistMatrix_s* A, ElInt height, ElInt width, ElInt ldim );
void ElDistMatrixResizeWithLDim_d
( ElDistMatrix_d* A, ElInt height, ElInt width, ElInt ldim );
void ElDistMatrixResizeWithLDim_c
( ElDistMatrix_c* A, ElInt height, ElInt width, ElInt ldim );
void ElDistMatrixResizeWithLDim_z
( ElDistMatrix_z* A, ElInt height, ElInt width, ElInt ldim );

float  ElDistMatrixGet_s( const ElDistMatrix_s* A, ElInt i, ElInt j );
double ElDistMatrixGet_d( const ElDistMatrix_d* A, ElInt i, ElInt j );
void   ElDistMatrixGet_c
( const ElDistMatrix_c* A, ElInt i, ElInt j, void* alpha );
void   ElDistMatrixGet_z
( const ElDistMatrix_z* A, ElInt i, ElInt j, void* alpha );

void ElDistMatrixSet_s( ElDistMatrix_s* A, ElInt i, ElInt j, float alpha );
void ElDistMatrixSet_d( ElDistMatrix_d* A, ElInt i, ElInt j, double alpha );
void ElDistMatrixSet_c( ElDistMatrix_c* A, ElInt i, ElInt j, void* alpha );
void ElDistMatrixSet_z( ElDistMatrix_z* A, ElInt i, ElInt j, void* alpha );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_DISTMATRIX_C_H */
