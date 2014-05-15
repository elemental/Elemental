/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_MATRIX_C_H
#define EL_MATRIX_C_H

#ifdef __cplusplus
extern "C" {
#endif

struct ElMatrix_s; typedef struct ElMatrix_s ElMatrix_s;
struct ElMatrix_d; typedef struct ElMatrix_d ElMatrix_d;
struct ElMatrix_c; typedef struct ElMatrix_c ElMatrix_c;
struct ElMatrix_z; typedef struct ElMatrix_z ElMatrix_z;

ElMatrix_s* ElMatrixCreate_s();
ElMatrix_d* ElMatrixCreate_d();
ElMatrix_c* ElMatrixCreate_c();
ElMatrix_z* ElMatrixCreate_z();

void ElMatrixDestroy_s( const ElMatrix_s* A );
void ElMatrixDestroy_d( const ElMatrix_d* A );
void ElMatrixDestroy_c( const ElMatrix_c* A );
void ElMatrixDestroy_z( const ElMatrix_z* A );

void ElMatrixEmpty_s( ElMatrix_s* A );
void ElMatrixEmpty_d( ElMatrix_d* A );
void ElMatrixEmpty_c( ElMatrix_c* A );
void ElMatrixEmpty_z( ElMatrix_z* A );

void ElMatrixResize_s( ElMatrix_s* A, ElInt height, ElInt width );
void ElMatrixResize_d( ElMatrix_d* A, ElInt height, ElInt width );
void ElMatrixResize_c( ElMatrix_c* A, ElInt height, ElInt width );
void ElMatrixResize_z( ElMatrix_z* A, ElInt height, ElInt width );

void ElMatrixResizeWithLDim_s
( ElMatrix_s* A, ElInt height, ElInt width, ElInt ldim );
void ElMatrixResizeWithLDim_d
( ElMatrix_d* A, ElInt height, ElInt width, ElInt ldim );
void ElMatrixResizeWithLDim_c
( ElMatrix_c* A, ElInt height, ElInt width, ElInt ldim );
void ElMatrixResizeWithLDim_z
( ElMatrix_z* A, ElInt height, ElInt width, ElInt ldim );

void ElMatrixAttach_s
( ElMatrix_s* A, ElInt height, ElInt width, float* buffer, ElInt ldim );
void ElMatrixAttach_d
( ElMatrix_d* A, ElInt height, ElInt width, double* buffer, ElInt ldim );
void ElMatrixAttach_c
( ElMatrix_c* A, ElInt height, ElInt width, void* buffer, ElInt ldim );
void ElMatrixAttach_z
( ElMatrix_z* A, ElInt height, ElInt width, void* buffer, ElInt ldim );

void ElMatrixLockedAttach_s
( ElMatrix_s* A, ElInt height, ElInt width, const float* buffer, ElInt ldim );
void ElMatrixLockedAttach_d
( ElMatrix_d* A, ElInt height, ElInt width, const double* buffer, ElInt ldim );
void ElMatrixLockedAttach_c
( ElMatrix_c* A, ElInt height, ElInt width, const void* buffer, ElInt ldim );
void ElMatrixLockedAttach_z
( ElMatrix_z* A, ElInt height, ElInt width, const void* buffer, ElInt ldim );

void ElMatrixControl_s
( ElMatrix_s* A, ElInt height, ElInt width, float* buffer, ElInt ldim );
void ElMatrixControl_d
( ElMatrix_d* A, ElInt height, ElInt width, double* buffer, ElInt ldim );
void ElMatrixControl_c
( ElMatrix_c* A, ElInt height, ElInt width, void* buffer, ElInt ldim );
void ElMatrixControl_z
( ElMatrix_z* A, ElInt height, ElInt width, void* buffer, ElInt ldim );

void ElMatrixCopy_s( const ElMatrix_s* A, ElMatrix_s* B );
void ElMatrixCopy_d( const ElMatrix_d* A, ElMatrix_d* B );
void ElMatrixCopy_c( const ElMatrix_c* A, ElMatrix_c* B );
void ElMatrixCopy_z( const ElMatrix_z* A, ElMatrix_z* B );

ElInt ElMatrixHeight_s( const ElMatrix_s* A );
ElInt ElMatrixHeight_d( const ElMatrix_d* A );
ElInt ElMatrixHeight_c( const ElMatrix_c* A );
ElInt ElMatrixHeight_z( const ElMatrix_z* A );

ElInt ElMatrixWidth_s( const ElMatrix_s* A );
ElInt ElMatrixWidth_d( const ElMatrix_d* A );
ElInt ElMatrixWidth_c( const ElMatrix_c* A );
ElInt ElMatrixWidth_z( const ElMatrix_z* A );

ElInt ElMatrixLDim_s( const ElMatrix_s* A );
ElInt ElMatrixLDim_d( const ElMatrix_d* A );
ElInt ElMatrixLDim_c( const ElMatrix_c* A );
ElInt ElMatrixLDim_z( const ElMatrix_z* A );

ElInt ElMatrixMemorySize_s( const ElMatrix_s* A );
ElInt ElMatrixMemorySize_d( const ElMatrix_d* A );
ElInt ElMatrixMemorySize_c( const ElMatrix_c* A );
ElInt ElMatrixMemorySize_z( const ElMatrix_z* A );

ElInt ElMatrixDiagonalLength_s( const ElMatrix_s* A, ElInt offset );
ElInt ElMatrixDiagonalLength_d( const ElMatrix_d* A, ElInt offset );
ElInt ElMatrixDiagonalLength_c( const ElMatrix_c* A, ElInt offset );
ElInt ElMatrixDiagonalLength_z( const ElMatrix_z* A, ElInt offset );

float*  ElMatrixBuffer_s( ElMatrix_s* A );
double* ElMatrixBuffer_d( ElMatrix_d* A );
void*   ElMatrixBuffer_c( ElMatrix_c* A );
void*   ElMatrixBuffer_z( ElMatrix_z* A );

const float*  ElMatrixLockedBuffer_s( const ElMatrix_s* A );
const double* ElMatrixLockedBuffer_d( const ElMatrix_d* A );
const void*   ElMatrixLockedBuffer_c( const ElMatrix_c* A );
const void*   ElMatrixLockedBuffer_z( const ElMatrix_z* A );

bool ElMatrixViewing_s( const ElMatrix_s* A );
bool ElMatrixViewing_d( const ElMatrix_d* A );
bool ElMatrixViewing_c( const ElMatrix_c* A );
bool ElMatrixViewing_z( const ElMatrix_z* A );

bool ElMatrixFixedSize_s( const ElMatrix_s* A );
bool ElMatrixFixedSize_d( const ElMatrix_d* A );
bool ElMatrixFixedSize_c( const ElMatrix_c* A );
bool ElMatrixFixedSize_z( const ElMatrix_z* A );

bool ElMatrixLocked_s( const ElMatrix_s* A );
bool ElMatrixLocked_d( const ElMatrix_d* A );
bool ElMatrixLocked_c( const ElMatrix_c* A );
bool ElMatrixLocked_z( const ElMatrix_z* A );

float  ElMatrixGet_s( const ElMatrix_s* A, ElInt i, ElInt j );
double ElMatrixGet_d( const ElMatrix_d* A, ElInt i, ElInt j );
void   ElMatrixGet_c( const ElMatrix_c* A, ElInt i, ElInt j, void* alpha );
void   ElMatrixGet_z( const ElMatrix_z* A, ElInt i, ElInt j, void* alpha );

float  ElMatrixGetRealPart_c( const ElMatrix_c* A, ElInt i, ElInt j );
double ElMatrixGetRealPart_z( const ElMatrix_z* A, ElInt i, ElInt j );

float  ElMatrixGetImagPart_c( const ElMatrix_c* A, ElInt i, ElInt j );
double ElMatrixGetImagPart_z( const ElMatrix_z* A, ElInt i, ElInt j );

void ElMatrixSet_s( ElMatrix_s* A, ElInt i, ElInt j, float alpha );
void ElMatrixSet_d( ElMatrix_d* A, ElInt i, ElInt j, double alpha );
void ElMatrixSet_c( ElMatrix_c* A, ElInt i, ElInt j, void* alpha );
void ElMatrixSet_z( ElMatrix_z* A, ElInt i, ElInt j, void* alpha );

void ElMatrixSetRealPart_c( ElMatrix_c* A, ElInt i, ElInt j, float alpha );
void ElMatrixSetRealPart_z( ElMatrix_z* A, ElInt i, ElInt j, double alpha );

void ElMatrixSetImagPart_c( ElMatrix_c* A, ElInt i, ElInt j, float alpha );
void ElMatrixSetImagPart_z( ElMatrix_z* A, ElInt i, ElInt j, double alpha );

void ElMatrixUpdate_s( ElMatrix_s* A, ElInt i, ElInt j, float alpha );
void ElMatrixUpdate_d( ElMatrix_d* A, ElInt i, ElInt j, double alpha );
void ElMatrixUpdate_c( ElMatrix_c* A, ElInt i, ElInt j, void* alpha );
void ElMatrixUpdate_z( ElMatrix_z* A, ElInt i, ElInt j, void* alpha );

void ElMatrixUpdateRealPart_c( ElMatrix_c* A, ElInt i, ElInt j, float alpha );
void ElMatrixUpdateRealPart_z( ElMatrix_z* A, ElInt i, ElInt j, double alpha );

void ElMatrixUpdateImagPart_c( ElMatrix_c* A, ElInt i, ElInt j, float alpha );
void ElMatrixUpdateImagPart_z( ElMatrix_z* A, ElInt i, ElInt j, double alpha );

void ElMatrixMakeReal_c( ElMatrix_c* A, ElInt i, ElInt j );
void ElMatrixMakeReal_z( ElMatrix_z* A, ElInt i, ElInt j );

void ElMatrixConjugate_c( ElMatrix_c* A, ElInt i, ElInt j );
void ElMatrixConjugate_z( ElMatrix_z* A, ElInt i, ElInt j );

ElMatrix_s* ElMatrixGetDiagonal_s( const ElMatrix_s* A, ElInt offset );
ElMatrix_d* ElMatrixGetDiagonal_d( const ElMatrix_d* A, ElInt offset );
ElMatrix_c* ElMatrixGetDiagonal_c( const ElMatrix_c* A, ElInt offset );
ElMatrix_z* ElMatrixGetDiagonal_z( const ElMatrix_z* A, ElInt offset );

ElMatrix_s* ElMatrixGetRealPartOfDiagonal_c
( const ElMatrix_c* A, ElInt offset );
ElMatrix_d* ElMatrixGetRealPartOfDiagonal_z
( const ElMatrix_z* A, ElInt offset );

ElMatrix_s* ElMatrixGetImagPartOfDiagonal_c
( const ElMatrix_c* A, ElInt offset );
ElMatrix_d* ElMatrixGetImagPartOfDiagonal_z
( const ElMatrix_z* A, ElInt offset );

void ElMatrixSetDiagonal_s
( ElMatrix_s* A, const ElMatrix_s* d, ElInt offset );
void ElMatrixSetDiagonal_d
( ElMatrix_d* A, const ElMatrix_d* d, ElInt offset );
void ElMatrixSetDiagonal_c
( ElMatrix_c* A, const ElMatrix_c* d, ElInt offset );
void ElMatrixSetDiagonal_z
( ElMatrix_z* A, const ElMatrix_z* d, ElInt offset );

void ElMatrixSetRealPartOfDiagonal_c
( ElMatrix_c* A, const ElMatrix_s* d, ElInt offset );
void ElMatrixSetRealPartOfDiagonal_z
( ElMatrix_z* A, const ElMatrix_d* d, ElInt offset );

void ElMatrixSetImagPartOfDiagonal_c
( ElMatrix_c* A, const ElMatrix_s* d, ElInt offset );
void ElMatrixSetImagPartOfDiagonal_z
( ElMatrix_z* A, const ElMatrix_d* d, ElInt offset );

void ElMatrixUpdateDiagonal_s
( ElMatrix_s* A, const ElMatrix_s* d, ElInt offset );
void ElMatrixUpdateDiagonal_d
( ElMatrix_d* A, const ElMatrix_d* d, ElInt offset );
void ElMatrixUpdateDiagonal_c
( ElMatrix_c* A, const ElMatrix_c* d, ElInt offset );
void ElMatrixUpdateDiagonal_z
( ElMatrix_z* A, const ElMatrix_z* d, ElInt offset );

void ElMatrixUpdateRealPartOfDiagonal_c
( ElMatrix_c* A, const ElMatrix_s* d, ElInt offset );
void ElMatrixUpdateRealPartOfDiagonal_z
( ElMatrix_z* A, const ElMatrix_d* d, ElInt offset );

void ElMatrixUpdateImagPartOfDiagonal_c
( ElMatrix_c* A, const ElMatrix_s* d, ElInt offset );
void ElMatrixUpdateImagPartOfDiagonal_z
( ElMatrix_z* A, const ElMatrix_d* d, ElInt offset );

void ElMatrixMakeDiagonalReal_c( ElMatrix_c* A, ElInt offset );
void ElMatrixMakeDiagonalReal_z( ElMatrix_z* A, ElInt offset );

void ElMatrixConjugateDiagonal_c( ElMatrix_c* A, ElInt offset );
void ElMatrixConjugateDiagonal_z( ElMatrix_z* A, ElInt offset );

ElMatrix_s* ElMatrixGetSubmatrix_s
( const ElMatrix_s* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );
ElMatrix_d* ElMatrixGetSubmatrix_d
( const ElMatrix_d* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );
ElMatrix_c* ElMatrixGetSubmatrix_c
( const ElMatrix_c* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );
ElMatrix_z* ElMatrixGetSubmatrix_z
( const ElMatrix_z* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );

ElMatrix_s* ElMatrixGetRealPartOfSubmatrix_c
( const ElMatrix_c* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );
ElMatrix_d* ElMatrixGetRealPartOfSubmatrix_z
( const ElMatrix_z* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );

ElMatrix_s* ElMatrixGetImagPartOfSubmatrix_c
( const ElMatrix_c* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );
ElMatrix_d* ElMatrixGetImagPartOfSubmatrix_z
( const ElMatrix_z* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );

void ElMatrixSetSubmatrix_s
( ElMatrix_s* A, ElInt* rowInds, ElInt* colInds, const ElMatrix_s* ASub );
void ElMatrixSetSubmatrix_d
( ElMatrix_d* A, ElInt* rowInds, ElInt* colInds, const ElMatrix_d* ASub );
void ElMatrixSetSubmatrix_c
( ElMatrix_c* A, ElInt* rowInds, ElInt* colInds, const ElMatrix_c* ASub );
void ElMatrixSetSubmatrix_z
( ElMatrix_z* A, ElInt* rowInds, ElInt* colInds, const ElMatrix_z* ASub );

void ElMatrixSetRealPartOfSubmatrix_c
( ElMatrix_c* A, ElInt* rowInds, ElInt* colInds, const ElMatrix_s* ASub );
void ElMatrixSetRealPartOfSubmatrix_z
( ElMatrix_z* A, ElInt* rowInds, ElInt* colInds, const ElMatrix_d* ASub );

void ElMatrixSetImagPartOfSubmatrix_c
( ElMatrix_c* A, ElInt* rowInds, ElInt* colInds, const ElMatrix_s* ASub );
void ElMatrixSetImagPartOfSubmatrix_z
( ElMatrix_z* A, ElInt* rowInds, ElInt* colInds, const ElMatrix_d* ASub );

void ElMatrixUpdateSubmatrix_s
( ElMatrix_s* A, ElInt* rowInds, ElInt* colInds, 
  float alpha, const ElMatrix_s* ASub );
void ElMatrixUpdateSubmatrix_d
( ElMatrix_d* A, ElInt* rowInds, ElInt* colInds, 
  double alpha, const ElMatrix_d* ASub );
void ElMatrixUpdateSubmatrix_c
( ElMatrix_c* A, ElInt* rowInds, ElInt* colInds, 
  void* alpha, const ElMatrix_c* ASub );
void ElMatrixUpdateSubmatrix_z
( ElMatrix_z* A, ElInt* rowInds, ElInt* colInds, 
  void* alpha, const ElMatrix_z* ASub );

void ElMatrixUpdateRealPartOfSubmatrix_c
( ElMatrix_c* A, ElInt* rowInds, ElInt* colInds, 
  void* alpha, const ElMatrix_s* ASub );
void ElMatrixUpdateRealPartOfSubmatrix_z
( ElMatrix_z* A, ElInt* rowInds, ElInt* colInds, 
  void* alpha, const ElMatrix_d* ASub );

void ElMatrixUpdateImagPartOfSubmatrix_c
( ElMatrix_c* A, ElInt* rowInds, ElInt* colInds, 
  void* alpha, const ElMatrix_s* ASub );
void ElMatrixUpdateImagPartOfSubmatrix_z
( ElMatrix_z* A, ElInt* rowInds, ElInt* colInds, 
  void* alpha, const ElMatrix_d* ASub );

void ElMatrixMakeSubmatrixReal_c
( ElMatrix_c* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );
void ElMatrixMakeSubmatrixReal_z
( ElMatrix_z* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );

void ElMatrixConjugateSubmatrix_c
( ElMatrix_c* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );
void ElMatrixConjugateSubmatrix_z
( ElMatrix_z* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_MATRIX_C_H */
