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

// An anonymous struct meant as a placeholder for DynamicDistMatrix<T>
// -------------------------------------------------------------------
struct ElDistMatrix_s; typedef struct ElDistMatrix_s ElDistMatrix_s;
struct ElDistMatrix_d; typedef struct ElDistMatrix_d ElDistMatrix_d;
struct ElDistMatrix_c; typedef struct ElDistMatrix_c ElDistMatrix_c;
struct ElDistMatrix_z; typedef struct ElDistMatrix_z ElDistMatrix_z;

// DistMatrix<T,MC,MR>::DistMatrix( const Grid& g )
// ------------------------------------------------
ElDistMatrix_s* ElDistMatrixCreate_s( const ElGrid* g );
ElDistMatrix_d* ElDistMatrixCreate_d( const ElGrid* g );
ElDistMatrix_c* ElDistMatrixCreate_c( const ElGrid* g );
ElDistMatrix_z* ElDistMatrixCreate_z( const ElGrid* g );

// DistMatrix<T,U,V>::DistMatrix( const Grid& g )
// ----------------------------------------------
ElDistMatrix_s* ElDistMatrixCreateSpecific_s
( ElDist U, ElDist V, const ElGrid* g );
ElDistMatrix_d* ElDistMatrixCreateSpecific_d
( ElDist U, ElDist V, const ElGrid* g );
ElDistMatrix_c* ElDistMatrixCreateSpecific_c
( ElDist U, ElDist V, const ElGrid* g );
ElDistMatrix_z* ElDistMatrixCreateSpecific_z
( ElDist U, ElDist V, const ElGrid* g );

// DistMatrix<T,U,V>::~DistMatrix()
// --------------------------------
void ElDistMatrixDestroy_s( const ElDistMatrix_s* A );
void ElDistMatrixDestroy_d( const ElDistMatrix_d* A );
void ElDistMatrixDestroy_c( const ElDistMatrix_c* A );
void ElDistMatrixDestroy_z( const ElDistMatrix_z* A );

// void DistMatrix<T,U,V>::Empty()
// -------------------------------
void ElDistMatrixEmpty_s( ElDistMatrix_s* A );
void ElDistMatrixEmpty_d( ElDistMatrix_d* A );
void ElDistMatrixEmpty_c( ElDistMatrix_c* A );
void ElDistMatrixEmpty_z( ElDistMatrix_z* A );

// void DistMatrix<T,U,V>::EmptyData()
// -----------------------------------
void ElDistMatrixEmptyData_s( ElDistMatrix_s* A );
void ElDistMatrixEmptyData_d( ElDistMatrix_d* A );
void ElDistMatrixEmptyData_c( ElDistMatrix_c* A );
void ElDistMatrixEmptyData_z( ElDistMatrix_z* A );

// void DistMatrix<T,U,V>::SetGrid( const Grid& g )
// ------------------------------------------------
void ElDistMatrixSetGrid_s( ElDistMatrix_s* AHandle, const ElGrid* gridHandle );
void ElDistMatrixSetGrid_d( ElDistMatrix_d* AHandle, const ElGrid* gridHandle );
void ElDistMatrixSetGrid_c( ElDistMatrix_c* AHandle, const ElGrid* gridHandle );
void ElDistMatrixSetGrid_z( ElDistMatrix_z* AHandle, const ElGrid* gridHandle );

// B = A
// -----
void ElDistMatrixCopy_s
( const ElDistMatrix_s* AHandle, ElDistMatrix_s* BHandle );
void ElDistMatrixCopy_d
( const ElDistMatrix_d* AHandle, ElDistMatrix_d* BHandle );
void ElDistMatrixCopy_c
( const ElDistMatrix_c* AHandle, ElDistMatrix_c* BHandle );
void ElDistMatrixCopy_z
( const ElDistMatrix_z* AHandle, ElDistMatrix_z* BHandle );

// void DistMatrix<T,U,V>::Resize( Int height, Int width )
// -------------------------------------------------------
void ElDistMatrixResize_s( ElDistMatrix_s* A, ElInt height, ElInt width );
void ElDistMatrixResize_d( ElDistMatrix_d* A, ElInt height, ElInt width );
void ElDistMatrixResize_c( ElDistMatrix_c* A, ElInt height, ElInt width );
void ElDistMatrixResize_z( ElDistMatrix_z* A, ElInt height, ElInt width );

// void DistMatrix<T,U,V>::Resize( Int height, Int width, Int ldim )
// -----------------------------------------------------------------
void ElDistMatrixResizeWithLDim_s
( ElDistMatrix_s* A, ElInt height, ElInt width, ElInt ldim );
void ElDistMatrixResizeWithLDim_d
( ElDistMatrix_d* A, ElInt height, ElInt width, ElInt ldim );
void ElDistMatrixResizeWithLDim_c
( ElDistMatrix_c* A, ElInt height, ElInt width, ElInt ldim );
void ElDistMatrixResizeWithLDim_z
( ElDistMatrix_z* A, ElInt height, ElInt width, ElInt ldim );

// void DistMatrix<T,U,V>::MakeConsistent( bool includeViewers )
// -------------------------------------------------------------
void ElDistMatrixMakeConsistent_s
( ElDistMatrix_s* AHandle, bool includeViewers );
void ElDistMatrixMakeConsistent_d
( ElDistMatrix_d* AHandle, bool includeViewers );
void ElDistMatrixMakeConsistent_c
( ElDistMatrix_c* AHandle, bool includeViewers );
void ElDistMatrixMakeConsistent_z
( ElDistMatrix_z* AHandle, bool includeViewers );

// void DistMatrix<T,U,V>::MakeSizeConsistent( bool includeViewers )
// -----------------------------------------------------------------
void ElDistMatrixMakeSizeConsistent_s
( ElDistMatrix_s* AHandle, bool includeViewers );
void ElDistMatrixMakeSizeConsistent_d
( ElDistMatrix_d* AHandle, bool includeViewers );
void ElDistMatrixMakeSizeConsistent_c
( ElDistMatrix_c* AHandle, bool includeViewers );
void ElDistMatrixMakeSizeConsistent_z
( ElDistMatrix_z* AHandle, bool includeViewers );

// TODO: Fill in a large number of routines here
// =============================================

// const Grid& DistMatrix<T,U,V>::Grid() const
// -------------------------------------------
const ElGrid* ElDistMatrixGrid_s( const ElDistMatrix_s* AHandle );
const ElGrid* ElDistMatrixGrid_d( const ElDistMatrix_d* AHandle );
const ElGrid* ElDistMatrixGrid_c( const ElDistMatrix_c* AHandle );
const ElGrid* ElDistMatrixGrid_z( const ElDistMatrix_z* AHandle );

// T DistMatrix<T,U,V>::Get( Int i, Int j ) const
// ----------------------------------------------
float  ElDistMatrixGet_s
( const ElDistMatrix_s* A, ElInt i, ElInt j );
double ElDistMatrixGet_d
( const ElDistMatrix_d* A, ElInt i, ElInt j );
void   ElDistMatrixGet_c
( const ElDistMatrix_c* A, ElInt i, ElInt j, void* alpha );
void   ElDistMatrixGet_z
( const ElDistMatrix_z* A, ElInt i, ElInt j, void* alpha );

// Base<T> DistMatrix<T,U,V>::GetRealPart( Int i, Int j ) const
// ------------------------------------------------------------
float  ElDistMatrixGetRealPart_c( const ElDistMatrix_c* A, ElInt i, ElInt j );
double ElDistMatrixGetRealPart_z( const ElDistMatrix_z* A, ElInt i, ElInt j );

// Base<T> DistMatrix<T,U,V>::GetImagPart( Int i, Int j ) const
// ------------------------------------------------------------
float  ElDistMatrixGetImagPart_c( const ElDistMatrix_c* A, ElInt i, ElInt j );
double ElDistMatrixGetImagPart_z( const ElDistMatrix_z* A, ElInt i, ElInt j );

// void DistMatrix<T,U,V>::Set( Int i, Int j, T alpha )
// ----------------------------------------------------
void ElDistMatrixSet_s( ElDistMatrix_s* A, ElInt i, ElInt j, float alpha );
void ElDistMatrixSet_d( ElDistMatrix_d* A, ElInt i, ElInt j, double alpha );
void ElDistMatrixSet_c( ElDistMatrix_c* A, ElInt i, ElInt j, void* alpha );
void ElDistMatrixSet_z( ElDistMatrix_z* A, ElInt i, ElInt j, void* alpha );

// void DistMatrix<T,U,V>::SetRealPart( Int i, Int j, Base<T> alpha )
// ------------------------------------------------------------------
void ElDistMatrixSetRealPart_c
( ElDistMatrix_c* A, ElInt i, ElInt j, float alpha );
void ElDistMatrixSetRealPart_z
( ElDistMatrix_z* A, ElInt i, ElInt j, double alpha );

// void DistMatrix<T,U,V>::SetImagPart( Int i, Int j, Base<T> alpha )
// ------------------------------------------------------------------
void ElDistMatrixSetImagPart_c
( ElDistMatrix_c* A, ElInt i, ElInt j, float alpha );
void ElDistMatrixSetImagPart_z
( ElDistMatrix_z* A, ElInt i, ElInt j, double alpha );

// void DistMatrix<T,U,V>::Update( Int i, Int j, T alpha )
// -------------------------------------------------------
void ElDistMatrixUpdate_s( ElDistMatrix_s* A, ElInt i, ElInt j, float alpha );
void ElDistMatrixUpdate_d( ElDistMatrix_d* A, ElInt i, ElInt j, double alpha );
void ElDistMatrixUpdate_c( ElDistMatrix_c* A, ElInt i, ElInt j, void* alpha );
void ElDistMatrixUpdate_z( ElDistMatrix_z* A, ElInt i, ElInt j, void* alpha );

// void DistMatrix<T,U,V>::UpdateRealPart( Int i, Int j, Base<T> alpha )
// ---------------------------------------------------------------------
void ElDistMatrixUpdateRealPart_c
( ElDistMatrix_c* A, ElInt i, ElInt j, float alpha );
void ElDistMatrixUpdateRealPart_z
( ElDistMatrix_z* A, ElInt i, ElInt j, double alpha );

// void DistMatrix<T,U,V>::UpdateImagPart( Int i, Int j, Base<T> alpha )
// ---------------------------------------------------------------------
void ElDistMatrixUpdateImagPart_c
( ElDistMatrix_c* A, ElInt i, ElInt j, float alpha );
void ElDistMatrixUpdateImagPart_z
( ElDistMatrix_z* A, ElInt i, ElInt j, double alpha );

// void DistMatrix<T,U,V>::MakeReal( Int i, Int j )
// ------------------------------------------------
void ElDistMatrixMakeReal_c( ElDistMatrix_c* AHandle, ElInt i, ElInt j );
void ElDistMatrixMakeReal_z( ElDistMatrix_z* AHandle, ElInt i, ElInt j );

// void DistMatrix<T,U,V>::Conjugate( Int i, Int j )
// ------------------------------------------------
void ElDistMatrixConjugate_c( ElDistMatrix_c* AHandle, ElInt i, ElInt j );
void ElDistMatrixConjugate_z( ElDistMatrix_z* AHandle, ElInt i, ElInt j );

// DistMatrix<T,UDiag,VDiag> DistMatrix<T,U,V>::GetDiagonal( Int offset ) const
// ----------------------------------------------------------------------------
ElDistMatrix_s* ElDistMatrixGetDiagonal_s
( const ElDistMatrix_s* AHandle, ElInt offset );
ElDistMatrix_d* ElDistMatrixGetDiagonal_d
( const ElDistMatrix_d* AHandle, ElInt offset );
ElDistMatrix_c* ElDistMatrixGetDiagonal_c
( const ElDistMatrix_c* AHandle, ElInt offset );
ElDistMatrix_z* ElDistMatrixGetDiagonal_z
( const ElDistMatrix_z* AHandle, ElInt offset );

// TODO: More diagonal manipulation
// ================================

// DistMatrix<T,STAR,STAR> DistMatrix<T,U,V>::GetSubmatrix
// ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
// --------------------------------------------------------------------------
ElDistMatrix_s* ElDistMatrixGetSubmatrix_s
( const ElDistMatrix_s* AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds );
ElDistMatrix_d* ElDistMatrixGetSubmatrix_d
( const ElDistMatrix_d* AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds );
ElDistMatrix_c* ElDistMatrixGetSubmatrix_c
( const ElDistMatrix_c* AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds );
ElDistMatrix_z* ElDistMatrixGetSubmatrix_z
( const ElDistMatrix_z* AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_DISTMATRIX_C_H */
