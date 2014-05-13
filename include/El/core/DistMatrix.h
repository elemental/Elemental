/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISTMATRIX_CINT_H
#define EL_DISTMATRIX_CINT_H

#ifdef __cplusplus
extern "C" {
#endif

struct ElDistMatrix_s; typedef struct ElDistMatrix_s ElDistMatrix_s;
struct ElDistMatrix_d; typedef struct ElDistMatrix_d ElDistMatrix_d;
struct ElDistMatrix_c; typedef struct ElDistMatrix_c ElDistMatrix_c;
struct ElDistMatrix_z; typedef struct ElDistMatrix_z ElDistMatrix_z;

ElDistMatrix_s* ElDistMatrixCreate_s();
ElDistMatrix_d* ElDistMatrixCreate_d();
ElDistMatrix_c* ElDistMatrixCreate_c();
ElDistMatrix_z* ElDistMatrixCreate_z();

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

void ElDistMatrixAttach_s
( ElDistMatrix_s* A, ElInt height, ElInt width, float* buffer, ElInt ldim );
void ElDistMatrixAttach_d
( ElDistMatrix_d* A, ElInt height, ElInt width, double* buffer, ElInt ldim );
void ElDistMatrixAttach_c
( ElDistMatrix_c* A, ElInt height, ElInt width, void* buffer, ElInt ldim );
void ElDistMatrixAttach_z
( ElDistMatrix_z* A, ElInt height, ElInt width, void* buffer, ElInt ldim );

void ElDistMatrixLockedAttach_s
( ElDistMatrix_s* A, ElInt height, ElInt width, const float* buffer, ElInt ldim );
void ElDistMatrixLockedAttach_d
( ElDistMatrix_d* A, ElInt height, ElInt width, const double* buffer, ElInt ldim );
void ElDistMatrixLockedAttach_c
( ElDistMatrix_c* A, ElInt height, ElInt width, const void* buffer, ElInt ldim );
void ElDistMatrixLockedAttach_z
( ElDistMatrix_z* A, ElInt height, ElInt width, const void* buffer, ElInt ldim );

void ElDistMatrixControl_s
( ElDistMatrix_s* A, ElInt height, ElInt width, float* buffer, ElInt ldim );
void ElDistMatrixControl_d
( ElDistMatrix_d* A, ElInt height, ElInt width, double* buffer, ElInt ldim );
void ElDistMatrixControl_c
( ElDistMatrix_c* A, ElInt height, ElInt width, void* buffer, ElInt ldim );
void ElDistMatrixControl_z
( ElDistMatrix_z* A, ElInt height, ElInt width, void* buffer, ElInt ldim );

void ElDistMatrixCopy_s( const ElDistMatrix_s* A, ElDistMatrix_s* B );
void ElDistMatrixCopy_d( const ElDistMatrix_d* A, ElDistMatrix_d* B );
void ElDistMatrixCopy_c( const ElDistMatrix_c* A, ElDistMatrix_c* B );
void ElDistMatrixCopy_z( const ElDistMatrix_z* A, ElDistMatrix_z* B );

Int ElDistMatrixHeight_s( const ElDistMatrix_s* A );
Int ElDistMatrixHeight_d( const ElDistMatrix_d* A );
Int ElDistMatrixHeight_c( const ElDistMatrix_c* A );
Int ElDistMatrixHeight_z( const ElDistMatrix_z* A );

Int ElDistMatrixWidth_s( const ElDistMatrix_s* A );
Int ElDistMatrixWidth_d( const ElDistMatrix_d* A );
Int ElDistMatrixWidth_c( const ElDistMatrix_c* A );
Int ElDistMatrixWidth_z( const ElDistMatrix_z* A );

Int ElDistMatrixLDim_s( const ElDistMatrix_s* A );
Int ElDistMatrixLDim_d( const ElDistMatrix_d* A );
Int ElDistMatrixLDim_c( const ElDistMatrix_c* A );
Int ElDistMatrixLDim_z( const ElDistMatrix_z* A );

Int ElDistMatrixMemorySize_s( const ElDistMatrix_s* A );
Int ElDistMatrixMemorySize_d( const ElDistMatrix_d* A );
Int ElDistMatrixMemorySize_c( const ElDistMatrix_c* A );
Int ElDistMatrixMemorySize_z( const ElDistMatrix_z* A );

Int ElDistMatrixDiagonalLength_s( const ElDistMatrix_s* A, ElInt offset );
Int ElDistMatrixDiagonalLength_d( const ElDistMatrix_d* A, ElInt offset );
Int ElDistMatrixDiagonalLength_c( const ElDistMatrix_c* A, ElInt offset );
Int ElDistMatrixDiagonalLength_z( const ElDistMatrix_z* A, ElInt offset );

float*  ElDistMatrixBuffer_s( ElDistMatrix_s* A );
double* ElDistMatrixBuffer_d( ElDistMatrix_d* A );
void*   ElDistMatrixBuffer_c( ElDistMatrix_c* A );
void*   ElDistMatrixBuffer_z( ElDistMatrix_z* A );

const float*  ElDistMatrixLockedBuffer_s( const ElDistMatrix_s* A );
const double* ElDistMatrixLockedBuffer_d( const ElDistMatrix_d* A );
const void*   ElDistMatrixLockedBuffer_c( const ElDistMatrix_c* A );
const void*   ElDistMatrixLockedBuffer_z( const ElDistMatrix_z* A );

bool ElDistMatrixViewing_s( const ElDistMatrix_s* A );
bool ElDistMatrixViewing_d( const ElDistMatrix_d* A );
bool ElDistMatrixViewing_c( const ElDistMatrix_c* A );
bool ElDistMatrixViewing_z( const ElDistMatrix_z* A );

bool ElDistMatrixFixedSize_s( const ElDistMatrix_s* A );
bool ElDistMatrixFixedSize_d( const ElDistMatrix_d* A );
bool ElDistMatrixFixedSize_c( const ElDistMatrix_c* A );
bool ElDistMatrixFixedSize_z( const ElDistMatrix_z* A );

bool ElDistMatrixLocked_s( const ElDistMatrix_s* A );
bool ElDistMatrixLocked_d( const ElDistMatrix_d* A );
bool ElDistMatrixLocked_c( const ElDistMatrix_c* A );
bool ElDistMatrixLocked_z( const ElDistMatrix_z* A );

float ElDistMatrixGet_s( const ElDistMatrix_s* A, ElInt i, ElInt j );
double ElDistMatrixGet_d( const ElDistMatrix_d* A, ElInt i, ElInt j );
void ElDistMatrixGet_c( const ElDistMatrix_c* A, ElInt i, ElInt j, void* alpha );
void ElDistMatrixGet_z( const ElDistMatrix_z* A, ElInt i, ElInt j, void* alpha );

float  ElDistMatrixGetRealPart_c( const ElDistMatrix_c* A, ElInt i, ElInt j );
double ElDistMatrixGetRealPart_z( const ElDistMatrix_z* A, ElInt i, ElInt j );

float  ElDistMatrixGetImagPart_c( const ElDistMatrix_c* A, ElInt i, ElInt j );
double ElDistMatrixGetImagPart_z( const ElDistMatrix_z* A, ElInt i, ElInt j );

void ElDistMatrixSet_s( ElDistMatrix_s* A, ElInt i, ElInt j, float alpha );
void ElDistMatrixSet_d( ElDistMatrix_d* A, ElInt i, ElInt j, double alpha );
void ElDistMatrixSet_c( ElDistMatrix_c* A, ElInt i, ElInt j, void* alpha );
void ElDistMatrixSet_z( ElDistMatrix_z* A, ElInt i, ElInt j, void* alpha );

void ElDistMatrixSetRealPart_c( ElDistMatrix_c* A, ElInt i, ElInt j, float alpha );
void ElDistMatrixSetRealPart_z( ElDistMatrix_z* A, ElInt i, ElInt j, double alpha );

void ElDistMatrixSetImagPart_c( ElDistMatrix_c* A, ElInt i, ElInt j, float alpha );
void ElDistMatrixSetImagPart_z( ElDistMatrix_z* A, ElInt i, ElInt j, double alpha );

void ElDistMatrixUpdate_s( ElDistMatrix_s* A, ElInt i, ElInt j, float alpha );
void ElDistMatrixUpdate_d( ElDistMatrix_d* A, ElInt i, ElInt j, double alpha );
void ElDistMatrixUpdate_c( ElDistMatrix_c* A, ElInt i, ElInt j, void* alpha );
void ElDistMatrixUpdate_z( ElDistMatrix_z* A, ElInt i, ElInt j, void* alpha );

void ElDistMatrixUpdateRealPart_c( ElDistMatrix_c* A, ElInt i, ElInt j, float alpha );
void ElDistMatrixUpdateRealPart_z( ElDistMatrix_z* A, ElInt i, ElInt j, double alpha );

void ElDistMatrixUpdateImagPart_c( ElDistMatrix_c* A, ElInt i, ElInt j, float alpha );
void ElDistMatrixUpdateImagPart_z( ElDistMatrix_z* A, ElInt i, ElInt j, double alpha );

void ElDistMatrixMakeReal_c( ElDistMatrix_c* A, ElInt i, ElInt j );
void ElDistMatrixMakeReal_z( ElDistMatrix_z* A, ElInt i, ElInt j );

void ElDistMatrixConjugate_c( ElDistMatrix_c* A, ElInt i, ElInt j );
void ElDistMatrixConjugate_z( ElDistMatrix_z* A, ElInt i, ElInt j );

ElDistMatrix_s* ElDistMatrixGetDiagonal_s( const ElDistMatrix_s* A, ElInt offset );
ElDistMatrix_d* ElDistMatrixGetDiagonal_d( const ElDistMatrix_d* A, ElInt offset );
ElDistMatrix_c* ElDistMatrixGetDiagonal_c( const ElDistMatrix_c* A, ElInt offset );
ElDistMatrix_z* ElDistMatrixGetDiagonal_z( const ElDistMatrix_z* A, ElInt offset );

ElDistMatrix_s* ElDistMatrixGetRealPartOfDiagonal_c
( const ElDistMatrix_c* A, ElInt offset );
ElDistMatrix_d* ElDistMatrixGetImagPartOfDiagonal_z
( const ElDistMatrix_z* A, ElInt offset );

void ElDistMatrixSetDiagonal
( ElDistMatrix_s* A, const ElDistMatrix_s* d, ElInt offset );
void ElDistMatrixSetDiagonal
( ElDistMatrix_d* A, const ElDistMatrix_d* d, ElInt offset );
void ElDistMatrixSetDiagonal
( ElDistMatrix_c* A, const ElDistMatrix_c* d, ElInt offset );
void ElDistMatrixSetDiagonal
( ElDistMatrix_z* A, const ElDistMatrix_z* d, ElInt offset );

void ElDistMatrixSetRealPartOfDiagonal
( ElDistMatrix_c* A, const ElDistMatrix_s* d, ElInt offset );
void ElDistMatrixSetRealPartOfDiagonal
( ElDistMatrix_z* A, const ElDistMatrix_d* d, ElInt offset );

void ElDistMatrixSetImagPartOfDiagonal
( ElDistMatrix_c* A, const ElDistMatrix_s* d, ElInt offset );
void ElDistMatrixSetImagPartOfDiagonal
( ElDistMatrix_z* A, const ElDistMatrix_d* d, ElInt offset );

void ElDistMatrixUpdateDiagonal_s
( ElDistMatrix_s* A, const ElDistMatrix_s* d, ElInt offset );
void ElDistMatrixUpdateDiagonal_d
( ElDistMatrix_d* A, const ElDistMatrix_d* d, ElInt offset );
void ElDistMatrixUpdateDiagonal_c
( ElDistMatrix_c* A, const ElDistMatrix_c* d, ElInt offset );
void ElDistMatrixUpdateDiagonal_z
( ElDistMatrix_z* A, const ElDistMatrix_z* d, ElInt offset );

void ElDistMatrixUpdateRealPartOfDiagonal_c
( ElDistMatrix_c* A, const ElDistMatrix_s* d, ElInt offset );
void ElDistMatrixUpdateRealPartOfDiagonal_z
( ElDistMatrix_z* A, const ElDistMatrix_d* d, ElInt offset );

void ElDistMatrixUpdateImagPartOfDiagonal_c
( ElDistMatrix_c* A, const ElDistMatrix_s* d, ElInt offset );
void ElDistMatrixUpdateImagPartOfDiagonal_z
( ElDistMatrix_z* A, const ElDistMatrix_d* d, ElInt offset );

void ElDistMatrixMakeDiagonalReal_c( ElDistMatrix_c* A, ElInt offset );
void ElDistMatrixMakeDiagonalReal_z( ElDistMatrix_z* A, ElInt offset );

void ElDistMatrixConjugateDiagonal_c( ElDistMatrix_c* A, ElInt offset );
void ElDistMatrixConjugateDiagonal_z( ElDistMatrix_z* A, ElInt offset );

ElDistMatrix_s* ElDistMatrixGetSubmatrix_s
( const ElDistMatrix_s* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );
ElDistMatrix_d* ElDistMatrixGetSubmatrix_d
( const ElDistMatrix_d* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );
ElDistMatrix_c* ElDistMatrixGetSubmatrix_c
( const ElDistMatrix_c* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );
ElDistMatrix_z* ElDistMatrixGetSubmatrix_z
( const ElDistMatrix_z* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );

ElDistMatrix_s* ElDistMatrixGetRealPartOfSubmatrix_c
( const ElDistMatrix_c* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );
ElDistMatrix_d* ElDistMatrixGetRealPartOfSubmatrix_z
( const ElDistMatrix_z* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );

ElDistMatrix_s* ElDistMatrixGetImagPartOfSubmatrix_c
( const ElDistMatrix_c* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );
ElDistMatrix_d* ElDistMatrixGetImagPartOfSubmatrix_z
( const ElDistMatrix_z* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );

void ElDistMatrixSetSubmatrix_s
( ElDistMatrix_s* A, ElInt* rowInds, ElInt* colInds, const ElDistMatrix_s* ASub );
void ElDistMatrixSetSubmatrix_d
( ElDistMatrix_d* A, ElInt* rowInds, ElInt* colInds, const ElDistMatrix_d* ASub );
void ElDistMatrixSetSubmatrix_c
( ElDistMatrix_c* A, ElInt* rowInds, ElInt* colInds, const ElDistMatrix_c* ASub );
void ElDistMatrixSetSubmatrix_z
( ElDistMatrix_z* A, ElInt* rowInds, ElInt* colInds, const ElDistMatrix_z* ASub );

void ElDistMatrixSetRealPartOfSubmatrix_c
( ElDistMatrix_c* A, ElInt* rowInds, ElInt* colInds, const ElDistMatrix_s* ASub );
void ElDistMatrixSetRealPartOfSubmatrix_z
( ElDistMatrix_z* A, ElInt* rowInds, ElInt* colInds, const ElDistMatrix_d* ASub );

void ElDistMatrixSetImagPartOfSubmatrix_c
( ElDistMatrix_c* A, ElInt* rowInds, ElInt* colInds, const ElDistMatrix_s* ASub );
void ElDistMatrixSetImagPartOfSubmatrix_z
( ElDistMatrix_z* A, ElInt* rowInds, ElInt* colInds, const ElDistMatrix_d* ASub );

void ElDistMatrixUpdateSubmatrix_s
( ElDistMatrix_s* A, ElInt* rowInds, ElInt* colInds, 
  float alpha, const ElDistMatrix_s* ASub );
void ElDistMatrixUpdateSubmatrix_d
( ElDistMatrix_d* A, ElInt* rowInds, ElInt* colInds, 
  double alpha, const ElDistMatrix_d* ASub );
void ElDistMatrixUpdateSubmatrix_c
( ElDistMatrix_c* A, ElInt* rowInds, ElInt* colInds, 
  void* alpha, const ElDistMatrix_c* ASub );
void ElDistMatrixUpdateSubmatrix_z
( ElDistMatrix_z* A, ElInt* rowInds, ElInt* colInds, 
  void* alpha, const ElDistMatrix_z* ASub );

void ElDistMatrixUpdateRealPartOfSubmatrix_c
( ElDistMatrix_c* A, ElInt* rowInds, ElInt* colInds, 
  void* alpha, const ElDistMatrix_s* ASub );
void ElDistMatrixUpdateRealPartOfSubmatrix_z
( ElDistMatrix_z* A, ElInt* rowInds, ElInt* colInds, 
  void* alpha, const ElDistMatrix_d* ASub );

void ElDistMatrixUpdateImagPartOfSubmatrix_c
( ElDistMatrix_c* A, ElInt* rowInds, ElInt* colInds, 
  void* alpha, const ElDistMatrix_s* ASub );
void ElDistMatrixUpdateImagPartOfSubmatrix_z
( ElDistMatrix_z* A, ElInt* rowInds, ElInt* colInds, 
  void* alpha, const ElDistMatrix_d* ASub );

void ElDistMatrixMakeSubmatrixReal_c
( ElDistMatrix_c* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );
void ElDistMatrixMakeSubmatrixReal_z
( ElDistMatrix_z* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );

void ElDistMatrixConjugateSubmatrix_c
( ElDistMatrix_c* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );
void ElDistMatrixConjugateSubmatrix_z
( ElDistMatrix_z* A, 
  ElInt numRowInds, ElInt* rowInds, ElInt numColInds, ElInt* colInds );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_DISTMATRIX_CINT_H */
