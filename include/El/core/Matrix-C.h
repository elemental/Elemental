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

/* An anonymous struct meant as a placeholder for Matrix<T>
   -------------------------------------------------------- */
typedef struct ElMatrix_sDummy* ElMatrix_s;
typedef struct ElMatrix_dDummy* ElMatrix_d;
typedef struct ElMatrix_cDummy* ElMatrix_c;
typedef struct ElMatrix_zDummy* ElMatrix_z;

typedef const struct ElMatrix_sDummy* ElConstMatrix_s;
typedef const struct ElMatrix_dDummy* ElConstMatrix_d;
typedef const struct ElMatrix_cDummy* ElConstMatrix_c;
typedef const struct ElMatrix_zDummy* ElConstMatrix_z;

/* Matrix<T>::Matrix()
   ------------------- */
ElMatrix_s ElMatrixCreate_s();
ElMatrix_d ElMatrixCreate_d();
ElMatrix_c ElMatrixCreate_c();
ElMatrix_z ElMatrixCreate_z();

/* Matrix<T>::~Matrix() 
   -------------------- */
void ElMatrixDestroy_s( ElConstMatrix_s A );
void ElMatrixDestroy_d( ElConstMatrix_d A );
void ElMatrixDestroy_c( ElConstMatrix_c A );
void ElMatrixDestroy_z( ElConstMatrix_z A );

/* void Matrix<T>::Empty()
   ----------------------- */
void ElMatrixEmpty_s( ElMatrix_s A );
void ElMatrixEmpty_d( ElMatrix_d A );
void ElMatrixEmpty_c( ElMatrix_c A );
void ElMatrixEmpty_z( ElMatrix_z A );

/* void Matrix<T>::Resize( Int height, Int width )
   ----------------------------------------------- */
void ElMatrixResize_s( ElMatrix_s A, ElInt height, ElInt width );
void ElMatrixResize_d( ElMatrix_d A, ElInt height, ElInt width );
void ElMatrixResize_c( ElMatrix_c A, ElInt height, ElInt width );
void ElMatrixResize_z( ElMatrix_z A, ElInt height, ElInt width );

/* void Matrix<T>::Resize( Int height, Int width, Int ldim )
   --------------------------------------------------------- */
void ElMatrixResizeWithLDim_s
( ElMatrix_s A, ElInt height, ElInt width, ElInt ldim );
void ElMatrixResizeWithLDim_d
( ElMatrix_d A, ElInt height, ElInt width, ElInt ldim );
void ElMatrixResizeWithLDim_c
( ElMatrix_c A, ElInt height, ElInt width, ElInt ldim );
void ElMatrixResizeWithLDim_z
( ElMatrix_z A, ElInt height, ElInt width, ElInt ldim );

/* void Matrix<T>::Attach( Int height, Int width, T* buffer, Int ldim )
   -------------------------------------------------------------------- */
void ElMatrixAttach_s
( ElMatrix_s A, ElInt height, ElInt width, float* buffer, ElInt ldim );
void ElMatrixAttach_d
( ElMatrix_d A, ElInt height, ElInt width, double* buffer, ElInt ldim );
void ElMatrixAttach_c
( ElMatrix_c A, ElInt height, ElInt width, 
  complex_float* buffer, ElInt ldim );
void ElMatrixAttach_z
( ElMatrix_z A, ElInt height, ElInt width, 
  complex_double* buffer, ElInt ldim );

/* void Matrix<T>::LockedAttach
   ( Int height, Int width, const T* buffer, Int ldim )
   ---------------------------------------------------- */
void ElMatrixLockedAttach_s
( ElMatrix_s A, ElInt height, ElInt width, const float* buffer, ElInt ldim );
void ElMatrixLockedAttach_d
( ElMatrix_d A, ElInt height, ElInt width, const double* buffer, ElInt ldim );
void ElMatrixLockedAttach_c
( ElMatrix_c A, ElInt height, ElInt width, 
  const complex_float* buffer, ElInt ldim );
void ElMatrixLockedAttach_z
( ElMatrix_z A, ElInt height, ElInt width, 
  const complex_double* buffer, ElInt ldim );

/* void Matrix<T>::Control( Int height, Int width, T* buffer, Int ldim )
   --------------------------------------------------------------------- */
void ElMatrixControl_s
( ElMatrix_s A, ElInt height, ElInt width, float* buffer, ElInt ldim );
void ElMatrixControl_d
( ElMatrix_d A, ElInt height, ElInt width, double* buffer, ElInt ldim );
void ElMatrixControl_c
( ElMatrix_c A, ElInt height, ElInt width, 
  complex_float* buffer, ElInt ldim );
void ElMatrixControl_z
( ElMatrix_z A, ElInt height, ElInt width, 
  complex_double* buffer, ElInt ldim );

/* B := A
   ------ */
void ElMatrixCopy_s( ElConstMatrix_s A, ElMatrix_s B );
void ElMatrixCopy_d( ElConstMatrix_d A, ElMatrix_d B );
void ElMatrixCopy_c( ElConstMatrix_c A, ElMatrix_c B );
void ElMatrixCopy_z( ElConstMatrix_z A, ElMatrix_z B );

/* Int Matrix<T>::Height() const
   ----------------------------- */
ElInt ElMatrixHeight_s( ElConstMatrix_s A );
ElInt ElMatrixHeight_d( ElConstMatrix_d A );
ElInt ElMatrixHeight_c( ElConstMatrix_c A );
ElInt ElMatrixHeight_z( ElConstMatrix_z A );

/* Int Matrix<T>::Width() const
   ---------------------------- */
ElInt ElMatrixWidth_s( ElConstMatrix_s A );
ElInt ElMatrixWidth_d( ElConstMatrix_d A );
ElInt ElMatrixWidth_c( ElConstMatrix_c A );
ElInt ElMatrixWidth_z( ElConstMatrix_z A );

/* Int Matrix<T>::LDim() const
   --------------------------- */
ElInt ElMatrixLDim_s( ElConstMatrix_s A );
ElInt ElMatrixLDim_d( ElConstMatrix_d A );
ElInt ElMatrixLDim_c( ElConstMatrix_c A );
ElInt ElMatrixLDim_z( ElConstMatrix_z A );

/* Int Matrix<T>::MemorySize() const
   --------------------------------- */
ElInt ElMatrixMemorySize_s( ElConstMatrix_s A );
ElInt ElMatrixMemorySize_d( ElConstMatrix_d A );
ElInt ElMatrixMemorySize_c( ElConstMatrix_c A );
ElInt ElMatrixMemorySize_z( ElConstMatrix_z A );

/* Int Matrix<T>::DiagonalLength( Int offset ) const
   ------------------------------------------------- */
ElInt ElMatrixDiagonalLength_s( ElConstMatrix_s A, ElInt offset );
ElInt ElMatrixDiagonalLength_d( ElConstMatrix_d A, ElInt offset );
ElInt ElMatrixDiagonalLength_c( ElConstMatrix_c A, ElInt offset );
ElInt ElMatrixDiagonalLength_z( ElConstMatrix_z A, ElInt offset );

/* T* Matrix<T>::Buffer()
   ---------------------- */
float*          ElMatrixBuffer_s( ElMatrix_s A );
double*         ElMatrixBuffer_d( ElMatrix_d A );
complex_float*  ElMatrixBuffer_c( ElMatrix_c A );
complex_double* ElMatrixBuffer_z( ElMatrix_z A );

/* const T* Matrix<T>::LockedBuffer() const
   ---------------------------------------- */
const float*          ElMatrixLockedBuffer_s( ElConstMatrix_s A );
const double*         ElMatrixLockedBuffer_d( ElConstMatrix_d A );
const complex_float*  ElMatrixLockedBuffer_c( ElConstMatrix_c A );
const complex_double* ElMatrixLockedBuffer_z( ElConstMatrix_z A );

/* bool Matrix<T>::Viewing() const
   ------------------------------- */
bool ElMatrixViewing_s( ElConstMatrix_s A );
bool ElMatrixViewing_d( ElConstMatrix_d A );
bool ElMatrixViewing_c( ElConstMatrix_c A );
bool ElMatrixViewing_z( ElConstMatrix_z A );

/* bool Matrix<T>::FixedSize() const
   --------------------------------- */
bool ElMatrixFixedSize_s( ElConstMatrix_s A );
bool ElMatrixFixedSize_d( ElConstMatrix_d A );
bool ElMatrixFixedSize_c( ElConstMatrix_c A );
bool ElMatrixFixedSize_z( ElConstMatrix_z A );

/* bool Matrix<T>::Locked() const
   ------------------------------ */
bool ElMatrixLocked_s( ElConstMatrix_s A );
bool ElMatrixLocked_d( ElConstMatrix_d A );
bool ElMatrixLocked_c( ElConstMatrix_c A );
bool ElMatrixLocked_z( ElConstMatrix_z A );

/* T Matrix<T>::Get( Int i, Int j ) const
   -------------------------------------- */
float          ElMatrixGet_s( ElConstMatrix_s A, ElInt i, ElInt j );
double         ElMatrixGet_d( ElConstMatrix_d A, ElInt i, ElInt j );
complex_float  ElMatrixGet_c( ElConstMatrix_c A, ElInt i, ElInt j );
complex_double ElMatrixGet_z( ElConstMatrix_z A, ElInt i, ElInt j );

/* Base<T> Matrix<T>::GetRealPart( Int i, Int j ) const
   ---------------------------------------------------- */
float  ElMatrixGetRealPart_c( ElConstMatrix_c A, ElInt i, ElInt j );
double ElMatrixGetRealPart_z( ElConstMatrix_z A, ElInt i, ElInt j );

/* Base<T> Matrix<T>::GetImagPart( Int i, Int j ) const
   ---------------------------------------------------- */
float  ElMatrixGetImagPart_c( ElConstMatrix_c A, ElInt i, ElInt j );
double ElMatrixGetImagPart_z( ElConstMatrix_z A, ElInt i, ElInt j );

/* void Matrix<T>::Set( Int i, Int j, T alpha )
   -------------------------------------------- */
void ElMatrixSet_s( ElMatrix_s A, ElInt i, ElInt j, float alpha );
void ElMatrixSet_d( ElMatrix_d A, ElInt i, ElInt j, double alpha );
void ElMatrixSet_c( ElMatrix_c A, ElInt i, ElInt j, complex_float alpha );
void ElMatrixSet_z( ElMatrix_z A, ElInt i, ElInt j, complex_double alpha );

/* void Matrix<T>::SetRealPart( Int i, Int j, Base<T> alpha )
   ---------------------------------------------------------- */
void ElMatrixSetRealPart_c( ElMatrix_c A, ElInt i, ElInt j, float alpha );
void ElMatrixSetRealPart_z( ElMatrix_z A, ElInt i, ElInt j, double alpha );

/* void Matrix<T>::SetImagPart( Int i, Int j, Base<T> alpha )
   ---------------------------------------------------------- */
void ElMatrixSetImagPart_c( ElMatrix_c A, ElInt i, ElInt j, float alpha );
void ElMatrixSetImagPart_z( ElMatrix_z A, ElInt i, ElInt j, double alpha );

/* void Matrix<T>::Update( Int i, Int j, T alpha )
   ----------------------------------------------- */
void ElMatrixUpdate_s( ElMatrix_s A, ElInt i, ElInt j, float alpha );
void ElMatrixUpdate_d( ElMatrix_d A, ElInt i, ElInt j, double alpha );
void ElMatrixUpdate_c( ElMatrix_c A, ElInt i, ElInt j, complex_float alpha );
void ElMatrixUpdate_z( ElMatrix_z A, ElInt i, ElInt j, complex_double alpha );

/* void Matrix<T>::UpdateRealPart( Int i, Int j, Base<T> alpha )
   ------------------------------------------------------------- */
void ElMatrixUpdateRealPart_c( ElMatrix_c A, ElInt i, ElInt j, float alpha );
void ElMatrixUpdateRealPart_z( ElMatrix_z A, ElInt i, ElInt j, double alpha );

/* void Matrix<T>::UpdateImagPart( Int i, Int j, Base<T> alpha )
   ------------------------------------------------------------- */
void ElMatrixUpdateImagPart_c( ElMatrix_c A, ElInt i, ElInt j, float alpha );
void ElMatrixUpdateImagPart_z( ElMatrix_z A, ElInt i, ElInt j, double alpha );

/* void Matrix<T>::MakeReal( Int i, Int j )
   ---------------------------------------- */
void ElMatrixMakeReal_c( ElMatrix_c A, ElInt i, ElInt j );
void ElMatrixMakeReal_z( ElMatrix_z A, ElInt i, ElInt j );

/* void Matrix<T>::Conjugate( Int i, Int j )
   ----------------------------------------- */
void ElMatrixConjugate_c( ElMatrix_c A, ElInt i, ElInt j );
void ElMatrixConjugate_z( ElMatrix_z A, ElInt i, ElInt j );

/* Matrix<T> Matrix<T>::GetDiagonal( Int offset ) const
   ---------------------------------------------------- */
ElMatrix_s ElMatrixGetDiagonal_s( ElConstMatrix_s A, ElInt offset );
ElMatrix_d ElMatrixGetDiagonal_d( ElConstMatrix_d A, ElInt offset );
ElMatrix_c ElMatrixGetDiagonal_c( ElConstMatrix_c A, ElInt offset );
ElMatrix_z ElMatrixGetDiagonal_z( ElConstMatrix_z A, ElInt offset );

/* Matrix<Base<T>> Matrix<T>::GetRealPartOfDiagonal( Int offset ) const
   -------------------------------------------------------------------- */
ElMatrix_s ElMatrixGetRealPartOfDiagonal_c( ElConstMatrix_c A, ElInt offset );
ElMatrix_d ElMatrixGetRealPartOfDiagonal_z( ElConstMatrix_z A, ElInt offset );

/* Matrix<Base<T>> Matrix<T>::GetImagPartOfDiagonal( Int offset ) const
   -------------------------------------------------------------------- */
ElMatrix_s ElMatrixGetImagPartOfDiagonal_c( ElConstMatrix_c A, ElInt offset );
ElMatrix_d ElMatrixGetImagPartOfDiagonal_z( ElConstMatrix_z A, ElInt offset );

/* void Matrix<T>::SetDiagonal( const Matrix<T>& d, Int offset )
   ------------------------------------------------------------- */
void ElMatrixSetDiagonal_s
( ElMatrix_s A, ElConstMatrix_s d, ElInt offset );
void ElMatrixSetDiagonal_d
( ElMatrix_d A, ElConstMatrix_d d, ElInt offset );
void ElMatrixSetDiagonal_c
( ElMatrix_c A, ElConstMatrix_c d, ElInt offset );
void ElMatrixSetDiagonal_z
( ElMatrix_z A, ElConstMatrix_z d, ElInt offset );

/* void Matrix<T>::SetRealPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
   -----------------------------------------------------------------------------
*/
void ElMatrixSetRealPartOfDiagonal_c
( ElMatrix_c A, ElConstMatrix_s d, ElInt offset );
void ElMatrixSetRealPartOfDiagonal_z
( ElMatrix_z A, ElConstMatrix_d d, ElInt offset );

/* void Matrix<T>::SetImagPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
   -----------------------------------------------------------------------------
*/
void ElMatrixSetImagPartOfDiagonal_c
( ElMatrix_c A, ElConstMatrix_s d, ElInt offset );
void ElMatrixSetImagPartOfDiagonal_z
( ElMatrix_z A, ElConstMatrix_d d, ElInt offset );

/* void Matrix<T>::UpdateDiagonal( const Matrix<T>& d, Int offset )
   ---------------------------------------------------------------- */
void ElMatrixUpdateDiagonal_s( ElMatrix_s A, ElConstMatrix_s d, ElInt offset );
void ElMatrixUpdateDiagonal_d( ElMatrix_d A, ElConstMatrix_d d, ElInt offset );
void ElMatrixUpdateDiagonal_c( ElMatrix_c A, ElConstMatrix_c d, ElInt offset );
void ElMatrixUpdateDiagonal_z( ElMatrix_z A, ElConstMatrix_z d, ElInt offset );

/* void Matrix<T>::UpdateRealPartOfDiagonal
   ( const Matrix<Base<T>>& d, Int offset )
   ---------------------------------------- */
void ElMatrixUpdateRealPartOfDiagonal_c
( ElMatrix_c A, ElConstMatrix_s d, ElInt offset );
void ElMatrixUpdateRealPartOfDiagonal_z
( ElMatrix_z A, ElConstMatrix_d d, ElInt offset );

/* void Matrix<T>::UpdateImagPartOfDiagonal
   ( const Matrix<Base<T>>& d, Int offset )
   ---------------------------------------- */
void ElMatrixUpdateImagPartOfDiagonal_c
( ElMatrix_c A, ElConstMatrix_s d, ElInt offset );
void ElMatrixUpdateImagPartOfDiagonal_z
( ElMatrix_z A, ElConstMatrix_d d, ElInt offset );

/* void Matrix<T>::MakeDiaogonalReal( Int offset )
   ----------------------------------------------- */
void ElMatrixMakeDiagonalReal_c( ElMatrix_c A, ElInt offset );
void ElMatrixMakeDiagonalReal_z( ElMatrix_z A, ElInt offset );

/* void Matrix<T>::ConjugateDiagonal Int offset )
   ---------------------------------------------- */
void ElMatrixConjugateDiagonal_c( ElMatrix_c A, ElInt offset );
void ElMatrixConjugateDiagonal_z( ElMatrix_z A, ElInt offset );

/* Matrix<T> Matrix<T>::GetSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
   -------------------------------------------------------------------------- */
ElMatrix_s ElMatrixGetSubmatrix_s
( ElConstMatrix_s A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );
ElMatrix_d ElMatrixGetSubmatrix_d
( ElConstMatrix_d A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );
ElMatrix_c ElMatrixGetSubmatrix_c
( ElConstMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );
ElMatrix_z ElMatrixGetSubmatrix_z
( ElConstMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );

/* Matrix<Base<T>> Matrix<T>::GetRealPartOfSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
   -------------------------------------------------------------------------- */
ElMatrix_s ElMatrixGetRealPartOfSubmatrix_c
( ElConstMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );
ElMatrix_d ElMatrixGetRealPartOfSubmatrix_z
( ElConstMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );

/* Matrix<Base<T>> Matrix<T>::GetImagPartOfSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
   -------------------------------------------------------------------------- */
ElMatrix_s ElMatrixGetImagPartOfSubmatrix_c
( ElConstMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );
ElMatrix_d ElMatrixGetImagPartOfSubmatrix_z
( ElConstMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );

/* void Matrix<T>::SetSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds, 
     const Matrix<T>& ASub )
   ------------------------------------------------------------------- */
void ElMatrixSetSubmatrix_s
( ElMatrix_s A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_s ASub );
void ElMatrixSetSubmatrix_d
( ElMatrix_d A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_d ASub );
void ElMatrixSetSubmatrix_c
( ElMatrix_c A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_c ASub );
void ElMatrixSetSubmatrix_z
( ElMatrix_z A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_z ASub );

/* void Matrix<T>::SetRealPartOfSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds, 
     const Matrix<Base<T>>& ASub )
   ------------------------------------------------------------------- */
void ElMatrixSetRealPartOfSubmatrix_c
( ElMatrix_c A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_s ASub );
void ElMatrixSetRealPartOfSubmatrix_z
( ElMatrix_z A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_d ASub );

/* void Matrix<T>::SetImagPartOfSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds, 
     const Matrix<Base<T>>& ASub )
   ------------------------------------------------------------------- */
void ElMatrixSetImagPartOfSubmatrix_c
( ElMatrix_c A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_s ASub );
void ElMatrixSetImagPartOfSubmatrix_z
( ElMatrix_z A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_d ASub );

/* void Matrix<T>::UpdateSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds, 
     const Matrix<T>& ASub )
   ------------------------------------------------------------------- */
void ElMatrixUpdateSubmatrix_s
( ElMatrix_s A, const ElInt* rowInds, const ElInt* colInds, 
  float alpha, ElConstMatrix_s ASub );
void ElMatrixUpdateSubmatrix_d
( ElMatrix_d A, const ElInt* rowInds, const ElInt* colInds, 
  double alpha, ElConstMatrix_d ASub );
void ElMatrixUpdateSubmatrix_c
( ElMatrix_c A, const ElInt* rowInds, const ElInt* colInds, 
  complex_float alpha, ElConstMatrix_c ASub );
void ElMatrixUpdateSubmatrix_z
( ElMatrix_z A, const ElInt* rowInds, const ElInt* colInds, 
  complex_double alpha, ElConstMatrix_z ASub );

/* void Matrix<T>::UpdateRealPartOfSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds, 
     const Matrix<Base<T>>& ASub )
   ------------------------------------------------------------------- */
void ElMatrixUpdateRealPartOfSubmatrix_c
( ElMatrix_c A, const ElInt* rowInds, const ElInt* colInds, 
  complex_float alpha, ElConstMatrix_s ASub );
void ElMatrixUpdateRealPartOfSubmatrix_z
( ElMatrix_z A, const ElInt* rowInds, const ElInt* colInds, 
  complex_double alpha, ElConstMatrix_d ASub );

/* void Matrix<T>::UpdateImagPartOfSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds, 
     const Matrix<Base<T>>& ASub )
   ------------------------------------------------------------------- */
void ElMatrixUpdateImagPartOfSubmatrix_c
( ElMatrix_c A, const ElInt* rowInds, const ElInt* colInds, 
  complex_float alpha, ElConstMatrix_s ASub );
void ElMatrixUpdateImagPartOfSubmatrix_z
( ElMatrix_z A, const ElInt* rowInds, const ElInt* colInds, 
  complex_double alpha, ElConstMatrix_d ASub );

/* void Matrix<T>::MakeSubmatrixReal
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds )
   -------------------------------------------------------------------- */
void ElMatrixMakeSubmatrixReal_c
( ElMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );
void ElMatrixMakeSubmatrixReal_z
( ElMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );

/* void Matrix<T>::ConjugateSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds )
   -------------------------------------------------------------------- */
void ElMatrixConjugateSubmatrix_c
( ElMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );
void ElMatrixConjugateSubmatrix_z
( ElMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_MATRIX_C_H */
