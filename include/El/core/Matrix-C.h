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
ElError ElMatrixCreate_s( ElMatrix_s* A );
ElError ElMatrixCreate_d( ElMatrix_d* A );
ElError ElMatrixCreate_c( ElMatrix_c* A );
ElError ElMatrixCreate_z( ElMatrix_z* A );

/* Matrix<T>::~Matrix() 
   -------------------- */
ElError ElMatrixDestroy_s( ElConstMatrix_s A );
ElError ElMatrixDestroy_d( ElConstMatrix_d A );
ElError ElMatrixDestroy_c( ElConstMatrix_c A );
ElError ElMatrixDestroy_z( ElConstMatrix_z A );

/* void Matrix<T>::Empty()
   ----------------------- */
ElError ElMatrixEmpty_s( ElMatrix_s A );
ElError ElMatrixEmpty_d( ElMatrix_d A );
ElError ElMatrixEmpty_c( ElMatrix_c A );
ElError ElMatrixEmpty_z( ElMatrix_z A );

/* void Matrix<T>::Resize( Int height, Int width )
   ----------------------------------------------- */
ElError ElMatrixResize_s( ElMatrix_s A, ElInt height, ElInt width );
ElError ElMatrixResize_d( ElMatrix_d A, ElInt height, ElInt width );
ElError ElMatrixResize_c( ElMatrix_c A, ElInt height, ElInt width );
ElError ElMatrixResize_z( ElMatrix_z A, ElInt height, ElInt width );

/* void Matrix<T>::Resize( Int height, Int width, Int ldim )
   --------------------------------------------------------- */
ElError ElMatrixResizeWithLDim_s
( ElMatrix_s A, ElInt height, ElInt width, ElInt ldim );
ElError ElMatrixResizeWithLDim_d
( ElMatrix_d A, ElInt height, ElInt width, ElInt ldim );
ElError ElMatrixResizeWithLDim_c
( ElMatrix_c A, ElInt height, ElInt width, ElInt ldim );
ElError ElMatrixResizeWithLDim_z
( ElMatrix_z A, ElInt height, ElInt width, ElInt ldim );

/* void Matrix<T>::Attach( Int height, Int width, T* buffer, Int ldim )
   -------------------------------------------------------------------- */
ElError ElMatrixAttach_s
( ElMatrix_s A, ElInt height, ElInt width, float* buffer, ElInt ldim );
ElError ElMatrixAttach_d
( ElMatrix_d A, ElInt height, ElInt width, double* buffer, ElInt ldim );
ElError ElMatrixAttach_c
( ElMatrix_c A, ElInt height, ElInt width, 
  complex_float* buffer, ElInt ldim );
ElError ElMatrixAttach_z
( ElMatrix_z A, ElInt height, ElInt width, 
  complex_double* buffer, ElInt ldim );

/* void Matrix<T>::LockedAttach
   ( Int height, Int width, const T* buffer, Int ldim )
   ---------------------------------------------------- */
ElError ElMatrixLockedAttach_s
( ElMatrix_s A, ElInt height, ElInt width, const float* buffer, ElInt ldim );
ElError ElMatrixLockedAttach_d
( ElMatrix_d A, ElInt height, ElInt width, const double* buffer, ElInt ldim );
ElError ElMatrixLockedAttach_c
( ElMatrix_c A, ElInt height, ElInt width, 
  const complex_float* buffer, ElInt ldim );
ElError ElMatrixLockedAttach_z
( ElMatrix_z A, ElInt height, ElInt width, 
  const complex_double* buffer, ElInt ldim );

/* void Matrix<T>::Control( Int height, Int width, T* buffer, Int ldim )
   --------------------------------------------------------------------- */
ElError ElMatrixControl_s
( ElMatrix_s A, ElInt height, ElInt width, float* buffer, ElInt ldim );
ElError ElMatrixControl_d
( ElMatrix_d A, ElInt height, ElInt width, double* buffer, ElInt ldim );
ElError ElMatrixControl_c
( ElMatrix_c A, ElInt height, ElInt width, 
  complex_float* buffer, ElInt ldim );
ElError ElMatrixControl_z
( ElMatrix_z A, ElInt height, ElInt width, 
  complex_double* buffer, ElInt ldim );

/* B := A
   ------ */
ElError ElMatrixCopy_s( ElConstMatrix_s A, ElMatrix_s B );
ElError ElMatrixCopy_d( ElConstMatrix_d A, ElMatrix_d B );
ElError ElMatrixCopy_c( ElConstMatrix_c A, ElMatrix_c B );
ElError ElMatrixCopy_z( ElConstMatrix_z A, ElMatrix_z B );

/* Int Matrix<T>::Height() const
   ----------------------------- */
ElError ElMatrixHeight_s( ElConstMatrix_s A, ElInt* height );
ElError ElMatrixHeight_d( ElConstMatrix_d A, ElInt* height );
ElError ElMatrixHeight_c( ElConstMatrix_c A, ElInt* height );
ElError ElMatrixHeight_z( ElConstMatrix_z A, ElInt* height );

/* Int Matrix<T>::Width() const
   ---------------------------- */
ElError ElMatrixWidth_s( ElConstMatrix_s A, ElInt* width );
ElError ElMatrixWidth_d( ElConstMatrix_d A, ElInt* width );
ElError ElMatrixWidth_c( ElConstMatrix_c A, ElInt* width );
ElError ElMatrixWidth_z( ElConstMatrix_z A, ElInt* width );

/* Int Matrix<T>::LDim() const
   --------------------------- */
ElError ElMatrixLDim_s( ElConstMatrix_s A, ElInt* ldim );
ElError ElMatrixLDim_d( ElConstMatrix_d A, ElInt* ldim );
ElError ElMatrixLDim_c( ElConstMatrix_c A, ElInt* ldim );
ElError ElMatrixLDim_z( ElConstMatrix_z A, ElInt* ldim );

/* Int Matrix<T>::MemorySize() const
   --------------------------------- */
ElError ElMatrixMemorySize_s( ElConstMatrix_s A, ElInt* memSize );
ElError ElMatrixMemorySize_d( ElConstMatrix_d A, ElInt* memSize );
ElError ElMatrixMemorySize_c( ElConstMatrix_c A, ElInt* memSize );
ElError ElMatrixMemorySize_z( ElConstMatrix_z A, ElInt* memSize );

/* Int Matrix<T>::DiagonalLength( Int offset ) const
   ------------------------------------------------- */
ElError ElMatrixDiagonalLength_s
( ElConstMatrix_s A, ElInt offset, ElInt* length );
ElError ElMatrixDiagonalLength_d
( ElConstMatrix_d A, ElInt offset, ElInt* length );
ElError ElMatrixDiagonalLength_c
( ElConstMatrix_c A, ElInt offset, ElInt* length );
ElError ElMatrixDiagonalLength_z
( ElConstMatrix_z A, ElInt offset, ElInt* length );

/* T* Matrix<T>::Buffer()
   ---------------------- */
ElError ElMatrixBuffer_s( ElMatrix_s A, float** buffer );
ElError ElMatrixBuffer_d( ElMatrix_d A, double** buffer );
ElError ElMatrixBuffer_c( ElMatrix_c A, complex_float** buffer );
ElError ElMatrixBuffer_z( ElMatrix_z A, complex_double** buffer );

/* const T* Matrix<T>::LockedBuffer() const
   ---------------------------------------- */
ElError ElMatrixLockedBuffer_s
( ElConstMatrix_s A, const float** buffer );
ElError ElMatrixLockedBuffer_d
( ElConstMatrix_d A, const double** buffer );
ElError ElMatrixLockedBuffer_c
( ElConstMatrix_c A, const complex_float** buffer );
ElError ElMatrixLockedBuffer_z
( ElConstMatrix_z A, const complex_double** buffer );

/* bool Matrix<T>::Viewing() const
   ------------------------------- */
ElError ElMatrixViewing_s( ElConstMatrix_s A, bool* viewing );
ElError ElMatrixViewing_d( ElConstMatrix_d A, bool* viewing );
ElError ElMatrixViewing_c( ElConstMatrix_c A, bool* viewing );
ElError ElMatrixViewing_z( ElConstMatrix_z A, bool* viewing );

/* bool Matrix<T>::FixedSize() const
   --------------------------------- */
ElError ElMatrixFixedSize_s( ElConstMatrix_s A, bool* fixedSize );
ElError ElMatrixFixedSize_d( ElConstMatrix_d A, bool* fixedSize );
ElError ElMatrixFixedSize_c( ElConstMatrix_c A, bool* fixedSize );
ElError ElMatrixFixedSize_z( ElConstMatrix_z A, bool* fixedSize );

/* bool Matrix<T>::Locked() const
   ------------------------------ */
ElError ElMatrixLocked_s( ElConstMatrix_s A, bool* locked );
ElError ElMatrixLocked_d( ElConstMatrix_d A, bool* locked );
ElError ElMatrixLocked_c( ElConstMatrix_c A, bool* locked );
ElError ElMatrixLocked_z( ElConstMatrix_z A, bool* locked );

/* T Matrix<T>::Get( Int i, Int j ) const
   -------------------------------------- */
ElError ElMatrixGet_s
( ElConstMatrix_s A, ElInt i, ElInt j, float* val );
ElError ElMatrixGet_d
( ElConstMatrix_d A, ElInt i, ElInt j, double* val );
ElError ElMatrixGet_c
( ElConstMatrix_c A, ElInt i, ElInt j, complex_float* val );
ElError ElMatrixGet_z
( ElConstMatrix_z A, ElInt i, ElInt j, complex_double* val );

/* Base<T> Matrix<T>::GetRealPart( Int i, Int j ) const
   ---------------------------------------------------- */
ElError ElMatrixGetRealPart_c
( ElConstMatrix_c A, ElInt i, ElInt j, float* val );
ElError ElMatrixGetRealPart_z
( ElConstMatrix_z A, ElInt i, ElInt j, double* val );

/* Base<T> Matrix<T>::GetImagPart( Int i, Int j ) const
   ---------------------------------------------------- */
ElError ElMatrixGetImagPart_c
( ElConstMatrix_c A, ElInt i, ElInt j, float* val );
ElError ElMatrixGetImagPart_z
( ElConstMatrix_z A, ElInt i, ElInt j, double* val );

/* void Matrix<T>::Set( Int i, Int j, T alpha )
   -------------------------------------------- */
ElError ElMatrixSet_s( ElMatrix_s A, ElInt i, ElInt j, float alpha );
ElError ElMatrixSet_d( ElMatrix_d A, ElInt i, ElInt j, double alpha );
ElError ElMatrixSet_c( ElMatrix_c A, ElInt i, ElInt j, complex_float alpha );
ElError ElMatrixSet_z( ElMatrix_z A, ElInt i, ElInt j, complex_double alpha );

/* void Matrix<T>::SetRealPart( Int i, Int j, Base<T> alpha )
   ---------------------------------------------------------- */
ElError ElMatrixSetRealPart_c( ElMatrix_c A, ElInt i, ElInt j, float alpha );
ElError ElMatrixSetRealPart_z( ElMatrix_z A, ElInt i, ElInt j, double alpha );

/* void Matrix<T>::SetImagPart( Int i, Int j, Base<T> alpha )
   ---------------------------------------------------------- */
ElError ElMatrixSetImagPart_c( ElMatrix_c A, ElInt i, ElInt j, float alpha );
ElError ElMatrixSetImagPart_z( ElMatrix_z A, ElInt i, ElInt j, double alpha );

/* void Matrix<T>::Update( Int i, Int j, T alpha )
   ----------------------------------------------- */
ElError ElMatrixUpdate_s
( ElMatrix_s A, ElInt i, ElInt j, float alpha );
ElError ElMatrixUpdate_d
( ElMatrix_d A, ElInt i, ElInt j, double alpha );
ElError ElMatrixUpdate_c
( ElMatrix_c A, ElInt i, ElInt j, complex_float alpha );
ElError ElMatrixUpdate_z
( ElMatrix_z A, ElInt i, ElInt j, complex_double alpha );

/* void Matrix<T>::UpdateRealPart( Int i, Int j, Base<T> alpha )
   ------------------------------------------------------------- */
ElError ElMatrixUpdateRealPart_c
( ElMatrix_c A, ElInt i, ElInt j, float alpha );
ElError ElMatrixUpdateRealPart_z
( ElMatrix_z A, ElInt i, ElInt j, double alpha );

/* void Matrix<T>::UpdateImagPart( Int i, Int j, Base<T> alpha )
   ------------------------------------------------------------- */
ElError ElMatrixUpdateImagPart_c
( ElMatrix_c A, ElInt i, ElInt j, float alpha );
ElError ElMatrixUpdateImagPart_z
( ElMatrix_z A, ElInt i, ElInt j, double alpha );

/* void Matrix<T>::MakeReal( Int i, Int j )
   ---------------------------------------- */
ElError ElMatrixMakeReal_c( ElMatrix_c A, ElInt i, ElInt j );
ElError ElMatrixMakeReal_z( ElMatrix_z A, ElInt i, ElInt j );

/* void Matrix<T>::Conjugate( Int i, Int j )
   ----------------------------------------- */
ElError ElMatrixConjugate_c( ElMatrix_c A, ElInt i, ElInt j );
ElError ElMatrixConjugate_z( ElMatrix_z A, ElInt i, ElInt j );

/* Matrix<T> Matrix<T>::GetDiagonal( Int offset ) const
   ---------------------------------------------------- */
ElError ElMatrixGetDiagonal_s( ElConstMatrix_s A, ElInt offset, ElMatrix_s* d );
ElError ElMatrixGetDiagonal_d( ElConstMatrix_d A, ElInt offset, ElMatrix_d* d );
ElError ElMatrixGetDiagonal_c( ElConstMatrix_c A, ElInt offset, ElMatrix_c* d );
ElError ElMatrixGetDiagonal_z( ElConstMatrix_z A, ElInt offset, ElMatrix_z* d );

/* Matrix<Base<T>> Matrix<T>::GetRealPartOfDiagonal( Int offset ) const
   -------------------------------------------------------------------- */
ElError ElMatrixGetRealPartOfDiagonal_c
( ElConstMatrix_c A, ElInt offset, ElMatrix_s* d );
ElError ElMatrixGetRealPartOfDiagonal_z
( ElConstMatrix_z A, ElInt offset, ElMatrix_d* d );

/* Matrix<Base<T>> Matrix<T>::GetImagPartOfDiagonal( Int offset ) const
   -------------------------------------------------------------------- */
ElError ElMatrixGetImagPartOfDiagonal_c
( ElConstMatrix_c A, ElInt offset, ElMatrix_s* d );
ElError ElMatrixGetImagPartOfDiagonal_z
( ElConstMatrix_z A, ElInt offset, ElMatrix_d* d );

/* void Matrix<T>::SetDiagonal( const Matrix<T>& d, Int offset )
   ------------------------------------------------------------- */
ElError ElMatrixSetDiagonal_s
( ElMatrix_s A, ElConstMatrix_s d, ElInt offset );
ElError ElMatrixSetDiagonal_d
( ElMatrix_d A, ElConstMatrix_d d, ElInt offset );
ElError ElMatrixSetDiagonal_c
( ElMatrix_c A, ElConstMatrix_c d, ElInt offset );
ElError ElMatrixSetDiagonal_z
( ElMatrix_z A, ElConstMatrix_z d, ElInt offset );

/* void Matrix<T>::SetRealPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
   -----------------------------------------------------------------------------
*/
ElError ElMatrixSetRealPartOfDiagonal_c
( ElMatrix_c A, ElConstMatrix_s d, ElInt offset );
ElError ElMatrixSetRealPartOfDiagonal_z
( ElMatrix_z A, ElConstMatrix_d d, ElInt offset );

/* void Matrix<T>::SetImagPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
   -----------------------------------------------------------------------------
*/
ElError ElMatrixSetImagPartOfDiagonal_c
( ElMatrix_c A, ElConstMatrix_s d, ElInt offset );
ElError ElMatrixSetImagPartOfDiagonal_z
( ElMatrix_z A, ElConstMatrix_d d, ElInt offset );

/* void Matrix<T>::UpdateDiagonal( const Matrix<T>& d, Int offset )
   ---------------------------------------------------------------- */
ElError ElMatrixUpdateDiagonal_s
( ElMatrix_s A, ElConstMatrix_s d, ElInt offset );
ElError ElMatrixUpdateDiagonal_d
( ElMatrix_d A, ElConstMatrix_d d, ElInt offset );
ElError ElMatrixUpdateDiagonal_c
( ElMatrix_c A, ElConstMatrix_c d, ElInt offset );
ElError ElMatrixUpdateDiagonal_z
( ElMatrix_z A, ElConstMatrix_z d, ElInt offset );

/* void Matrix<T>::UpdateRealPartOfDiagonal
   ( const Matrix<Base<T>>& d, Int offset )
   ---------------------------------------- */
ElError ElMatrixUpdateRealPartOfDiagonal_c
( ElMatrix_c A, ElConstMatrix_s d, ElInt offset );
ElError ElMatrixUpdateRealPartOfDiagonal_z
( ElMatrix_z A, ElConstMatrix_d d, ElInt offset );

/* void Matrix<T>::UpdateImagPartOfDiagonal
   ( const Matrix<Base<T>>& d, Int offset )
   ---------------------------------------- */
ElError ElMatrixUpdateImagPartOfDiagonal_c
( ElMatrix_c A, ElConstMatrix_s d, ElInt offset );
ElError ElMatrixUpdateImagPartOfDiagonal_z
( ElMatrix_z A, ElConstMatrix_d d, ElInt offset );

/* void Matrix<T>::MakeDiaogonalReal( Int offset )
   ----------------------------------------------- */
ElError ElMatrixMakeDiagonalReal_c( ElMatrix_c A, ElInt offset );
ElError ElMatrixMakeDiagonalReal_z( ElMatrix_z A, ElInt offset );

/* void Matrix<T>::ConjugateDiagonal Int offset )
   ---------------------------------------------- */
ElError ElMatrixConjugateDiagonal_c( ElMatrix_c A, ElInt offset );
ElError ElMatrixConjugateDiagonal_z( ElMatrix_z A, ElInt offset );

/* Matrix<T> Matrix<T>::GetSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
   -------------------------------------------------------------------------- */
ElError ElMatrixGetSubmatrix_s
( ElConstMatrix_s A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds, ElMatrix_s* ASub );
ElError ElMatrixGetSubmatrix_d
( ElConstMatrix_d A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds, ElMatrix_d* ASub );
ElError ElMatrixGetSubmatrix_c
( ElConstMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds, ElMatrix_c* ASub );
ElError ElMatrixGetSubmatrix_z
( ElConstMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds, ElMatrix_z* ASub );

/* Matrix<Base<T>> Matrix<T>::GetRealPartOfSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
   -------------------------------------------------------------------------- */
ElError ElMatrixGetRealPartOfSubmatrix_c
( ElConstMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds, ElMatrix_s* ASub );
ElError ElMatrixGetRealPartOfSubmatrix_z
( ElConstMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds, ElMatrix_d* ASub );

/* Matrix<Base<T>> Matrix<T>::GetImagPartOfSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
   -------------------------------------------------------------------------- */
ElError ElMatrixGetImagPartOfSubmatrix_c
( ElConstMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds, ElMatrix_s* ASub );
ElError ElMatrixGetImagPartOfSubmatrix_z
( ElConstMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds, ElMatrix_d* ASub );

/* void Matrix<T>::SetSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds, 
     const Matrix<T>& ASub )
   ------------------------------------------------------------------- */
ElError ElMatrixSetSubmatrix_s
( ElMatrix_s A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_s ASub );
ElError ElMatrixSetSubmatrix_d
( ElMatrix_d A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_d ASub );
ElError ElMatrixSetSubmatrix_c
( ElMatrix_c A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_c ASub );
ElError ElMatrixSetSubmatrix_z
( ElMatrix_z A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_z ASub );

/* void Matrix<T>::SetRealPartOfSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds, 
     const Matrix<Base<T>>& ASub )
   ------------------------------------------------------------------- */
ElError ElMatrixSetRealPartOfSubmatrix_c
( ElMatrix_c A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_s ASub );
ElError ElMatrixSetRealPartOfSubmatrix_z
( ElMatrix_z A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_d ASub );

/* void Matrix<T>::SetImagPartOfSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds, 
     const Matrix<Base<T>>& ASub )
   ------------------------------------------------------------------- */
ElError ElMatrixSetImagPartOfSubmatrix_c
( ElMatrix_c A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_s ASub );
ElError ElMatrixSetImagPartOfSubmatrix_z
( ElMatrix_z A, 
  const ElInt* rowInds, const ElInt* colInds, ElConstMatrix_d ASub );

/* void Matrix<T>::UpdateSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds, 
     const Matrix<T>& ASub )
   ------------------------------------------------------------------- */
ElError ElMatrixUpdateSubmatrix_s
( ElMatrix_s A, const ElInt* rowInds, const ElInt* colInds, 
  float alpha, ElConstMatrix_s ASub );
ElError ElMatrixUpdateSubmatrix_d
( ElMatrix_d A, const ElInt* rowInds, const ElInt* colInds, 
  double alpha, ElConstMatrix_d ASub );
ElError ElMatrixUpdateSubmatrix_c
( ElMatrix_c A, const ElInt* rowInds, const ElInt* colInds, 
  complex_float alpha, ElConstMatrix_c ASub );
ElError ElMatrixUpdateSubmatrix_z
( ElMatrix_z A, const ElInt* rowInds, const ElInt* colInds, 
  complex_double alpha, ElConstMatrix_z ASub );

/* void Matrix<T>::UpdateRealPartOfSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds, 
     const Matrix<Base<T>>& ASub )
   ------------------------------------------------------------------- */
ElError ElMatrixUpdateRealPartOfSubmatrix_c
( ElMatrix_c A, const ElInt* rowInds, const ElInt* colInds, 
  complex_float alpha, ElConstMatrix_s ASub );
ElError ElMatrixUpdateRealPartOfSubmatrix_z
( ElMatrix_z A, const ElInt* rowInds, const ElInt* colInds, 
  complex_double alpha, ElConstMatrix_d ASub );

/* void Matrix<T>::UpdateImagPartOfSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds, 
     const Matrix<Base<T>>& ASub )
   ------------------------------------------------------------------- */
ElError ElMatrixUpdateImagPartOfSubmatrix_c
( ElMatrix_c A, const ElInt* rowInds, const ElInt* colInds, 
  complex_float alpha, ElConstMatrix_s ASub );
ElError ElMatrixUpdateImagPartOfSubmatrix_z
( ElMatrix_z A, const ElInt* rowInds, const ElInt* colInds, 
  complex_double alpha, ElConstMatrix_d ASub );

/* void Matrix<T>::MakeSubmatrixReal
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds )
   -------------------------------------------------------------------- */
ElError ElMatrixMakeSubmatrixReal_c
( ElMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );
ElError ElMatrixMakeSubmatrixReal_z
( ElMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );

/* void Matrix<T>::ConjugateSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds )
   -------------------------------------------------------------------- */
ElError ElMatrixConjugateSubmatrix_c
( ElMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );
ElError ElMatrixConjugateSubmatrix_z
( ElMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_MATRIX_C_H */
