/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_MATRIX_C_H
#define EL_MATRIX_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* NOTE: "ElInt* I" and "ElInt* J" have been converted to 
         "ElInt* rowInd" and "ElInt* colInd" due to what appears to be a bug
         in /usr/bin/c++ in Ubuntu */

/* An anonymous struct meant as a placeholder for Matrix<T>
   -------------------------------------------------------- */
typedef struct ElMatrix_iDummy* ElMatrix_i;
typedef struct ElMatrix_sDummy* ElMatrix_s;
typedef struct ElMatrix_dDummy* ElMatrix_d;
typedef struct ElMatrix_cDummy* ElMatrix_c;
typedef struct ElMatrix_zDummy* ElMatrix_z;

typedef const struct ElMatrix_iDummy* ElConstMatrix_i;
typedef const struct ElMatrix_sDummy* ElConstMatrix_s;
typedef const struct ElMatrix_dDummy* ElConstMatrix_d;
typedef const struct ElMatrix_cDummy* ElConstMatrix_c;
typedef const struct ElMatrix_zDummy* ElConstMatrix_z;

/* Matrix<T>::Matrix()
   ------------------- */
EL_EXPORT ElError ElMatrixCreate_i( ElMatrix_i* A );
EL_EXPORT ElError ElMatrixCreate_s( ElMatrix_s* A );
EL_EXPORT ElError ElMatrixCreate_d( ElMatrix_d* A );
EL_EXPORT ElError ElMatrixCreate_c( ElMatrix_c* A );
EL_EXPORT ElError ElMatrixCreate_z( ElMatrix_z* A );

/* Matrix<T>::~Matrix() 
   -------------------- */
EL_EXPORT ElError ElMatrixDestroy_i( ElConstMatrix_i A );
EL_EXPORT ElError ElMatrixDestroy_s( ElConstMatrix_s A );
EL_EXPORT ElError ElMatrixDestroy_d( ElConstMatrix_d A );
EL_EXPORT ElError ElMatrixDestroy_c( ElConstMatrix_c A );
EL_EXPORT ElError ElMatrixDestroy_z( ElConstMatrix_z A );

/* void Matrix<T>::Empty()
   ----------------------- */
EL_EXPORT ElError ElMatrixEmpty_i( ElMatrix_i A );
EL_EXPORT ElError ElMatrixEmpty_s( ElMatrix_s A );
EL_EXPORT ElError ElMatrixEmpty_d( ElMatrix_d A );
EL_EXPORT ElError ElMatrixEmpty_c( ElMatrix_c A );
EL_EXPORT ElError ElMatrixEmpty_z( ElMatrix_z A );

/* void Matrix<T>::Resize( Int height, Int width )
   ----------------------------------------------- */
EL_EXPORT ElError ElMatrixResize_i( ElMatrix_i A, ElInt height, ElInt width );
EL_EXPORT ElError ElMatrixResize_s( ElMatrix_s A, ElInt height, ElInt width );
EL_EXPORT ElError ElMatrixResize_d( ElMatrix_d A, ElInt height, ElInt width );
EL_EXPORT ElError ElMatrixResize_c( ElMatrix_c A, ElInt height, ElInt width );
EL_EXPORT ElError ElMatrixResize_z( ElMatrix_z A, ElInt height, ElInt width );

/* void Matrix<T>::Resize( Int height, Int width, Int ldim )
   --------------------------------------------------------- */
EL_EXPORT ElError ElMatrixResizeWithLDim_i
( ElMatrix_i A, ElInt height, ElInt width, ElInt ldim );
EL_EXPORT ElError ElMatrixResizeWithLDim_s
( ElMatrix_s A, ElInt height, ElInt width, ElInt ldim );
EL_EXPORT ElError ElMatrixResizeWithLDim_d
( ElMatrix_d A, ElInt height, ElInt width, ElInt ldim );
EL_EXPORT ElError ElMatrixResizeWithLDim_c
( ElMatrix_c A, ElInt height, ElInt width, ElInt ldim );
EL_EXPORT ElError ElMatrixResizeWithLDim_z
( ElMatrix_z A, ElInt height, ElInt width, ElInt ldim );

/* void Matrix<T>::Attach( Int height, Int width, T* buffer, Int ldim )
   -------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixAttach_i
( ElMatrix_i A, ElInt height, ElInt width, ElInt* buffer, ElInt ldim );
EL_EXPORT ElError ElMatrixAttach_s
( ElMatrix_s A, ElInt height, ElInt width, float* buffer, ElInt ldim );
EL_EXPORT ElError ElMatrixAttach_d
( ElMatrix_d A, ElInt height, ElInt width, double* buffer, ElInt ldim );
EL_EXPORT ElError ElMatrixAttach_c
( ElMatrix_c A, ElInt height, ElInt width, 
  complex_float* buffer, ElInt ldim );
EL_EXPORT ElError ElMatrixAttach_z
( ElMatrix_z A, ElInt height, ElInt width, 
  complex_double* buffer, ElInt ldim );

/* void Matrix<T>::LockedAttach
   ( Int height, Int width, const T* buffer, Int ldim )
   ---------------------------------------------------- */
EL_EXPORT ElError ElMatrixLockedAttach_i
( ElMatrix_i A, ElInt height, ElInt width, const ElInt* buffer, ElInt ldim );
EL_EXPORT ElError ElMatrixLockedAttach_s
( ElMatrix_s A, ElInt height, ElInt width, const float* buffer, ElInt ldim );
EL_EXPORT ElError ElMatrixLockedAttach_d
( ElMatrix_d A, ElInt height, ElInt width, const double* buffer, ElInt ldim );
EL_EXPORT ElError ElMatrixLockedAttach_c
( ElMatrix_c A, ElInt height, ElInt width, 
  const complex_float* buffer, ElInt ldim );
EL_EXPORT ElError ElMatrixLockedAttach_z
( ElMatrix_z A, ElInt height, ElInt width, 
  const complex_double* buffer, ElInt ldim );

/* void Matrix<T>::Control( Int height, Int width, T* buffer, Int ldim )
   --------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixControl_i
( ElMatrix_i A, ElInt height, ElInt width, ElInt* buffer, ElInt ldim );
EL_EXPORT ElError ElMatrixControl_s
( ElMatrix_s A, ElInt height, ElInt width, float* buffer, ElInt ldim );
EL_EXPORT ElError ElMatrixControl_d
( ElMatrix_d A, ElInt height, ElInt width, double* buffer, ElInt ldim );
EL_EXPORT ElError ElMatrixControl_c
( ElMatrix_c A, ElInt height, ElInt width, 
  complex_float* buffer, ElInt ldim );
EL_EXPORT ElError ElMatrixControl_z
( ElMatrix_z A, ElInt height, ElInt width, 
  complex_double* buffer, ElInt ldim );

/* B := A
   ------ */
EL_EXPORT ElError ElMatrixCopy_i( ElConstMatrix_i A, ElMatrix_i B );
EL_EXPORT ElError ElMatrixCopy_s( ElConstMatrix_s A, ElMatrix_s B );
EL_EXPORT ElError ElMatrixCopy_d( ElConstMatrix_d A, ElMatrix_d B );
EL_EXPORT ElError ElMatrixCopy_c( ElConstMatrix_c A, ElMatrix_c B );
EL_EXPORT ElError ElMatrixCopy_z( ElConstMatrix_z A, ElMatrix_z B );

/* Int Matrix<T>::Height() const
   ----------------------------- */
EL_EXPORT ElError ElMatrixHeight_i( ElConstMatrix_i A, ElInt* height );
EL_EXPORT ElError ElMatrixHeight_s( ElConstMatrix_s A, ElInt* height );
EL_EXPORT ElError ElMatrixHeight_d( ElConstMatrix_d A, ElInt* height );
EL_EXPORT ElError ElMatrixHeight_c( ElConstMatrix_c A, ElInt* height );
EL_EXPORT ElError ElMatrixHeight_z( ElConstMatrix_z A, ElInt* height );

/* Int Matrix<T>::Width() const
   ---------------------------- */
EL_EXPORT ElError ElMatrixWidth_i( ElConstMatrix_i A, ElInt* width );
EL_EXPORT ElError ElMatrixWidth_s( ElConstMatrix_s A, ElInt* width );
EL_EXPORT ElError ElMatrixWidth_d( ElConstMatrix_d A, ElInt* width );
EL_EXPORT ElError ElMatrixWidth_c( ElConstMatrix_c A, ElInt* width );
EL_EXPORT ElError ElMatrixWidth_z( ElConstMatrix_z A, ElInt* width );

/* Int Matrix<T>::LDim() const
   --------------------------- */
EL_EXPORT ElError ElMatrixLDim_i( ElConstMatrix_i A, ElInt* ldim );
EL_EXPORT ElError ElMatrixLDim_s( ElConstMatrix_s A, ElInt* ldim );
EL_EXPORT ElError ElMatrixLDim_d( ElConstMatrix_d A, ElInt* ldim );
EL_EXPORT ElError ElMatrixLDim_c( ElConstMatrix_c A, ElInt* ldim );
EL_EXPORT ElError ElMatrixLDim_z( ElConstMatrix_z A, ElInt* ldim );

/* Int Matrix<T>::MemorySize() const
   --------------------------------- */
EL_EXPORT ElError ElMatrixMemorySize_i( ElConstMatrix_i A, ElInt* memSize );
EL_EXPORT ElError ElMatrixMemorySize_s( ElConstMatrix_s A, ElInt* memSize );
EL_EXPORT ElError ElMatrixMemorySize_d( ElConstMatrix_d A, ElInt* memSize );
EL_EXPORT ElError ElMatrixMemorySize_c( ElConstMatrix_c A, ElInt* memSize );
EL_EXPORT ElError ElMatrixMemorySize_z( ElConstMatrix_z A, ElInt* memSize );

/* Int Matrix<T>::DiagonalLength( Int offset ) const
   ------------------------------------------------- */
EL_EXPORT ElError ElMatrixDiagonalLength_i
( ElConstMatrix_i A, ElInt offset, ElInt* length );
EL_EXPORT ElError ElMatrixDiagonalLength_s
( ElConstMatrix_s A, ElInt offset, ElInt* length );
EL_EXPORT ElError ElMatrixDiagonalLength_d
( ElConstMatrix_d A, ElInt offset, ElInt* length );
EL_EXPORT ElError ElMatrixDiagonalLength_c
( ElConstMatrix_c A, ElInt offset, ElInt* length );
EL_EXPORT ElError ElMatrixDiagonalLength_z
( ElConstMatrix_z A, ElInt offset, ElInt* length );

/* T* Matrix<T>::Buffer()
   ---------------------- */
EL_EXPORT ElError ElMatrixBuffer_i( ElMatrix_i A, ElInt** buffer );
EL_EXPORT ElError ElMatrixBuffer_s( ElMatrix_s A, float** buffer );
EL_EXPORT ElError ElMatrixBuffer_d( ElMatrix_d A, double** buffer );
EL_EXPORT ElError ElMatrixBuffer_c( ElMatrix_c A, complex_float** buffer );
EL_EXPORT ElError ElMatrixBuffer_z( ElMatrix_z A, complex_double** buffer );

/* const T* Matrix<T>::LockedBuffer() const
   ---------------------------------------- */
EL_EXPORT ElError ElMatrixLockedBuffer_i
( ElConstMatrix_i A, const ElInt** buffer );
EL_EXPORT ElError ElMatrixLockedBuffer_s
( ElConstMatrix_s A, const float** buffer );
EL_EXPORT ElError ElMatrixLockedBuffer_d
( ElConstMatrix_d A, const double** buffer );
EL_EXPORT ElError ElMatrixLockedBuffer_c
( ElConstMatrix_c A, const complex_float** buffer );
EL_EXPORT ElError ElMatrixLockedBuffer_z
( ElConstMatrix_z A, const complex_double** buffer );

/* bool Matrix<T>::Viewing() const
   ------------------------------- */
EL_EXPORT ElError ElMatrixViewing_i( ElConstMatrix_i A, bool* viewing );
EL_EXPORT ElError ElMatrixViewing_s( ElConstMatrix_s A, bool* viewing );
EL_EXPORT ElError ElMatrixViewing_d( ElConstMatrix_d A, bool* viewing );
EL_EXPORT ElError ElMatrixViewing_c( ElConstMatrix_c A, bool* viewing );
EL_EXPORT ElError ElMatrixViewing_z( ElConstMatrix_z A, bool* viewing );

/* bool Matrix<T>::FixedSize() const
   --------------------------------- */
EL_EXPORT ElError ElMatrixFixedSize_i( ElConstMatrix_i A, bool* fixedSize );
EL_EXPORT ElError ElMatrixFixedSize_s( ElConstMatrix_s A, bool* fixedSize );
EL_EXPORT ElError ElMatrixFixedSize_d( ElConstMatrix_d A, bool* fixedSize );
EL_EXPORT ElError ElMatrixFixedSize_c( ElConstMatrix_c A, bool* fixedSize );
EL_EXPORT ElError ElMatrixFixedSize_z( ElConstMatrix_z A, bool* fixedSize );

/* bool Matrix<T>::Locked() const
   ------------------------------ */
EL_EXPORT ElError ElMatrixLocked_i( ElConstMatrix_i A, bool* locked );
EL_EXPORT ElError ElMatrixLocked_s( ElConstMatrix_s A, bool* locked );
EL_EXPORT ElError ElMatrixLocked_d( ElConstMatrix_d A, bool* locked );
EL_EXPORT ElError ElMatrixLocked_c( ElConstMatrix_c A, bool* locked );
EL_EXPORT ElError ElMatrixLocked_z( ElConstMatrix_z A, bool* locked );

/* T Matrix<T>::Get( Int i, Int j ) const
   -------------------------------------- */
EL_EXPORT ElError ElMatrixGet_i
( ElConstMatrix_i A, ElInt i, ElInt j, ElInt* val );
EL_EXPORT ElError ElMatrixGet_s
( ElConstMatrix_s A, ElInt i, ElInt j, float* val );
EL_EXPORT ElError ElMatrixGet_d
( ElConstMatrix_d A, ElInt i, ElInt j, double* val );
EL_EXPORT ElError ElMatrixGet_c
( ElConstMatrix_c A, ElInt i, ElInt j, complex_float* val );
EL_EXPORT ElError ElMatrixGet_z
( ElConstMatrix_z A, ElInt i, ElInt j, complex_double* val );

/* Base<T> Matrix<T>::GetRealPart( Int i, Int j ) const
   ---------------------------------------------------- */
EL_EXPORT ElError ElMatrixGetRealPart_c
( ElConstMatrix_c A, ElInt i, ElInt j, float* val );
EL_EXPORT ElError ElMatrixGetRealPart_z
( ElConstMatrix_z A, ElInt i, ElInt j, double* val );

/* Base<T> Matrix<T>::GetImagPart( Int i, Int j ) const
   ---------------------------------------------------- */
EL_EXPORT ElError ElMatrixGetImagPart_c
( ElConstMatrix_c A, ElInt i, ElInt j, float* val );
EL_EXPORT ElError ElMatrixGetImagPart_z
( ElConstMatrix_z A, ElInt i, ElInt j, double* val );

/* void Matrix<T>::Set( Int i, Int j, T alpha )
   -------------------------------------------- */
EL_EXPORT ElError ElMatrixSet_s
( ElMatrix_s A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElMatrixSet_s
( ElMatrix_s A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElMatrixSet_d
( ElMatrix_d A, ElInt i, ElInt j, double alpha );
EL_EXPORT ElError ElMatrixSet_c
( ElMatrix_c A, ElInt i, ElInt j, complex_float alpha );
EL_EXPORT ElError ElMatrixSet_z
( ElMatrix_z A, ElInt i, ElInt j, complex_double alpha );

/* void Matrix<T>::SetRealPart( Int i, Int j, Base<T> alpha )
   ---------------------------------------------------------- */
EL_EXPORT ElError ElMatrixSetRealPart_c
( ElMatrix_c A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElMatrixSetRealPart_z
( ElMatrix_z A, ElInt i, ElInt j, double alpha );

/* void Matrix<T>::SetImagPart( Int i, Int j, Base<T> alpha )
   ---------------------------------------------------------- */
EL_EXPORT ElError ElMatrixSetImagPart_c
( ElMatrix_c A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElMatrixSetImagPart_z
( ElMatrix_z A, ElInt i, ElInt j, double alpha );

/* void Matrix<T>::Update( Int i, Int j, T alpha )
   ----------------------------------------------- */
EL_EXPORT ElError ElMatrixUpdate_i
( ElMatrix_i A, ElInt i, ElInt j, ElInt alpha );
EL_EXPORT ElError ElMatrixUpdate_s
( ElMatrix_s A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElMatrixUpdate_d
( ElMatrix_d A, ElInt i, ElInt j, double alpha );
EL_EXPORT ElError ElMatrixUpdate_c
( ElMatrix_c A, ElInt i, ElInt j, complex_float alpha );
EL_EXPORT ElError ElMatrixUpdate_z
( ElMatrix_z A, ElInt i, ElInt j, complex_double alpha );

/* void Matrix<T>::UpdateRealPart( Int i, Int j, Base<T> alpha )
   ------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixUpdateRealPart_c
( ElMatrix_c A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElMatrixUpdateRealPart_z
( ElMatrix_z A, ElInt i, ElInt j, double alpha );

/* void Matrix<T>::UpdateImagPart( Int i, Int j, Base<T> alpha )
   ------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixUpdateImagPart_c
( ElMatrix_c A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElMatrixUpdateImagPart_z
( ElMatrix_z A, ElInt i, ElInt j, double alpha );

/* void Matrix<T>::MakeReal( Int i, Int j )
   ---------------------------------------- */
EL_EXPORT ElError ElMatrixMakeReal_c( ElMatrix_c A, ElInt i, ElInt j );
EL_EXPORT ElError ElMatrixMakeReal_z( ElMatrix_z A, ElInt i, ElInt j );

/* void Matrix<T>::Conjugate( Int i, Int j )
   ----------------------------------------- */
EL_EXPORT ElError ElMatrixConjugate_c( ElMatrix_c A, ElInt i, ElInt j );
EL_EXPORT ElError ElMatrixConjugate_z( ElMatrix_z A, ElInt i, ElInt j );

/* Matrix<T> Matrix<T>::GetDiagonal( Int offset ) const
   ---------------------------------------------------- */
EL_EXPORT ElError ElMatrixGetDiagonal_i
( ElConstMatrix_i A, ElInt offset, ElMatrix_i* d );
EL_EXPORT ElError ElMatrixGetDiagonal_s
( ElConstMatrix_s A, ElInt offset, ElMatrix_s* d );
EL_EXPORT ElError ElMatrixGetDiagonal_d
( ElConstMatrix_d A, ElInt offset, ElMatrix_d* d );
EL_EXPORT ElError ElMatrixGetDiagonal_c
( ElConstMatrix_c A, ElInt offset, ElMatrix_c* d );
EL_EXPORT ElError ElMatrixGetDiagonal_z
( ElConstMatrix_z A, ElInt offset, ElMatrix_z* d );

/* Matrix<Base<T>> Matrix<T>::GetRealPartOfDiagonal( Int offset ) const
   -------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixGetRealPartOfDiagonal_c
( ElConstMatrix_c A, ElInt offset, ElMatrix_s* d );
EL_EXPORT ElError ElMatrixGetRealPartOfDiagonal_z
( ElConstMatrix_z A, ElInt offset, ElMatrix_d* d );

/* Matrix<Base<T>> Matrix<T>::GetImagPartOfDiagonal( Int offset ) const
   -------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixGetImagPartOfDiagonal_c
( ElConstMatrix_c A, ElInt offset, ElMatrix_s* d );
EL_EXPORT ElError ElMatrixGetImagPartOfDiagonal_z
( ElConstMatrix_z A, ElInt offset, ElMatrix_d* d );

/* void Matrix<T>::SetDiagonal( const Matrix<T>& d, Int offset )
   ------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixSetDiagonal_i
( ElMatrix_i A, ElConstMatrix_i d, ElInt offset );
EL_EXPORT ElError ElMatrixSetDiagonal_s
( ElMatrix_s A, ElConstMatrix_s d, ElInt offset );
EL_EXPORT ElError ElMatrixSetDiagonal_d
( ElMatrix_d A, ElConstMatrix_d d, ElInt offset );
EL_EXPORT ElError ElMatrixSetDiagonal_c
( ElMatrix_c A, ElConstMatrix_c d, ElInt offset );
EL_EXPORT ElError ElMatrixSetDiagonal_z
( ElMatrix_z A, ElConstMatrix_z d, ElInt offset );

/* void Matrix<T>::SetRealPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
   -----------------------------------------------------------------------------
*/
EL_EXPORT ElError ElMatrixSetRealPartOfDiagonal_c
( ElMatrix_c A, ElConstMatrix_s d, ElInt offset );
EL_EXPORT ElError ElMatrixSetRealPartOfDiagonal_z
( ElMatrix_z A, ElConstMatrix_d d, ElInt offset );

/* void Matrix<T>::SetImagPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
   -----------------------------------------------------------------------------
*/
EL_EXPORT ElError ElMatrixSetImagPartOfDiagonal_c
( ElMatrix_c A, ElConstMatrix_s d, ElInt offset );
EL_EXPORT ElError ElMatrixSetImagPartOfDiagonal_z
( ElMatrix_z A, ElConstMatrix_d d, ElInt offset );

/* void Matrix<T>::UpdateDiagonal( const Matrix<T>& d, Int offset )
   ---------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixUpdateDiagonal_i
( ElMatrix_i A, ElConstMatrix_i d, ElInt offset );
EL_EXPORT ElError ElMatrixUpdateDiagonal_s
( ElMatrix_s A, ElConstMatrix_s d, ElInt offset );
EL_EXPORT ElError ElMatrixUpdateDiagonal_d
( ElMatrix_d A, ElConstMatrix_d d, ElInt offset );
EL_EXPORT ElError ElMatrixUpdateDiagonal_c
( ElMatrix_c A, ElConstMatrix_c d, ElInt offset );
EL_EXPORT ElError ElMatrixUpdateDiagonal_z
( ElMatrix_z A, ElConstMatrix_z d, ElInt offset );

/* void Matrix<T>::UpdateRealPartOfDiagonal
   ( const Matrix<Base<T>>& d, Int offset )
   ---------------------------------------- */
EL_EXPORT ElError ElMatrixUpdateRealPartOfDiagonal_c
( ElMatrix_c A, ElConstMatrix_s d, ElInt offset );
EL_EXPORT ElError ElMatrixUpdateRealPartOfDiagonal_z
( ElMatrix_z A, ElConstMatrix_d d, ElInt offset );

/* void Matrix<T>::UpdateImagPartOfDiagonal
   ( const Matrix<Base<T>>& d, Int offset )
   ---------------------------------------- */
EL_EXPORT ElError ElMatrixUpdateImagPartOfDiagonal_c
( ElMatrix_c A, ElConstMatrix_s d, ElInt offset );
EL_EXPORT ElError ElMatrixUpdateImagPartOfDiagonal_z
( ElMatrix_z A, ElConstMatrix_d d, ElInt offset );

/* void Matrix<T>::MakeDiaogonalReal( Int offset )
   ----------------------------------------------- */
EL_EXPORT ElError ElMatrixMakeDiagonalReal_c( ElMatrix_c A, ElInt offset );
EL_EXPORT ElError ElMatrixMakeDiagonalReal_z( ElMatrix_z A, ElInt offset );

/* void Matrix<T>::ConjugateDiagonal Int offset )
   ---------------------------------------------- */
EL_EXPORT ElError ElMatrixConjugateDiagonal_c( ElMatrix_c A, ElInt offset );
EL_EXPORT ElError ElMatrixConjugateDiagonal_z( ElMatrix_z A, ElInt offset );

/* Matrix<T> Matrix<T>::GetSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J ) const
   -------------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixGetSubmatrix_i
( ElConstMatrix_i A, 
  ElInt numRowInds, const ElInt* rowInd, 
  ElInt numColInds, const ElInt* colInd, ElMatrix_i* ASub );
EL_EXPORT ElError ElMatrixGetSubmatrix_s
( ElConstMatrix_s A, 
  ElInt numRowInds, const ElInt* rowInd, 
  ElInt numColInds, const ElInt* colInd, ElMatrix_s* ASub );
EL_EXPORT ElError ElMatrixGetSubmatrix_d
( ElConstMatrix_d A, 
  ElInt numRowInds, const ElInt* rowInd, 
  ElInt numColInds, const ElInt* colInd, ElMatrix_d* ASub );
EL_EXPORT ElError ElMatrixGetSubmatrix_c
( ElConstMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInd, 
  ElInt numColInds, const ElInt* colInd, ElMatrix_c* ASub );
EL_EXPORT ElError ElMatrixGetSubmatrix_z
( ElConstMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInd, 
  ElInt numColInds, const ElInt* colInd, ElMatrix_z* ASub );

/* Matrix<Base<T>> Matrix<T>::GetRealPartOfSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J ) const
   -------------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixGetRealPartOfSubmatrix_c
( ElConstMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInd, 
  ElInt numColInds, const ElInt* colInd, ElMatrix_s* ASub );
EL_EXPORT ElError ElMatrixGetRealPartOfSubmatrix_z
( ElConstMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInd, 
  ElInt numColInds, const ElInt* colInd, ElMatrix_d* ASub );

/* Matrix<Base<T>> Matrix<T>::GetImagPartOfSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J ) const
   -------------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixGetImagPartOfSubmatrix_c
( ElConstMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInd, 
  ElInt numColInds, const ElInt* colInd, ElMatrix_s* ASub );
EL_EXPORT ElError ElMatrixGetImagPartOfSubmatrix_z
( ElConstMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInd, 
  ElInt numColInds, const ElInt* colInd, ElMatrix_d* ASub );

/* void Matrix<T>::SetSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J, 
     const Matrix<T>& ASub )
   ------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixSetSubmatrix_i
( ElMatrix_i A, const ElInt* rowInd, const ElInt* colInd, ElConstMatrix_i ASub );
EL_EXPORT ElError ElMatrixSetSubmatrix_s
( ElMatrix_s A, const ElInt* rowInd, const ElInt* colInd, ElConstMatrix_s ASub );
EL_EXPORT ElError ElMatrixSetSubmatrix_d
( ElMatrix_d A, const ElInt* rowInd, const ElInt* colInd, ElConstMatrix_d ASub );
EL_EXPORT ElError ElMatrixSetSubmatrix_c
( ElMatrix_c A, const ElInt* rowInd, const ElInt* colInd, ElConstMatrix_c ASub );
EL_EXPORT ElError ElMatrixSetSubmatrix_z
( ElMatrix_z A, const ElInt* rowInd, const ElInt* colInd, ElConstMatrix_z ASub );

/* void Matrix<T>::SetRealPartOfSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J, 
     const Matrix<Base<T>>& ASub )
   ------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixSetRealPartOfSubmatrix_c
( ElMatrix_c A, const ElInt* rowInd, const ElInt* colInd, ElConstMatrix_s ASub );
EL_EXPORT ElError ElMatrixSetRealPartOfSubmatrix_z
( ElMatrix_z A, const ElInt* rowInd, const ElInt* colInd, ElConstMatrix_d ASub );

/* void Matrix<T>::SetImagPartOfSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J, 
     const Matrix<Base<T>>& ASub )
   ------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixSetImagPartOfSubmatrix_c
( ElMatrix_c A, const ElInt* rowInd, const ElInt* colInd, ElConstMatrix_s ASub );
EL_EXPORT ElError ElMatrixSetImagPartOfSubmatrix_z
( ElMatrix_z A, const ElInt* rowInd, const ElInt* colInd, ElConstMatrix_d ASub );

/* void Matrix<T>::UpdateSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J, 
     const Matrix<T>& ASub )
   ------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixUpdateSubmatrix_i
( ElMatrix_i A, const ElInt* rowInd, const ElInt* colInd, 
  ElInt alpha, ElConstMatrix_i ASub );
EL_EXPORT ElError ElMatrixUpdateSubmatrix_s
( ElMatrix_s A, const ElInt* rowInd, const ElInt* colInd, 
  float alpha, ElConstMatrix_s ASub );
EL_EXPORT ElError ElMatrixUpdateSubmatrix_d
( ElMatrix_d A, const ElInt* rowInd, const ElInt* colInd, 
  double alpha, ElConstMatrix_d ASub );
EL_EXPORT ElError ElMatrixUpdateSubmatrix_c
( ElMatrix_c A, const ElInt* rowInd, const ElInt* colInd, 
  complex_float alpha, ElConstMatrix_c ASub );
EL_EXPORT ElError ElMatrixUpdateSubmatrix_z
( ElMatrix_z A, const ElInt* rowInd, const ElInt* colInd, 
  complex_double alpha, ElConstMatrix_z ASub );

/* void Matrix<T>::UpdateRealPartOfSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J, 
     Base<T> alpha, const Matrix<Base<T>>& ASub )
   ------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixUpdateRealPartOfSubmatrix_c
( ElMatrix_c A, const ElInt* rowInd, const ElInt* colInd, 
  float alpha, ElConstMatrix_s ASub );
EL_EXPORT ElError ElMatrixUpdateRealPartOfSubmatrix_z
( ElMatrix_z A, const ElInt* rowInd, const ElInt* colInd, 
  double alpha, ElConstMatrix_d ASub );

/* void Matrix<T>::UpdateImagPartOfSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J, 
     Base<T> alpha, const Matrix<Base<T>>& ASub )
   ------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixUpdateImagPartOfSubmatrix_c
( ElMatrix_c A, const ElInt* rowInd, const ElInt* colInd, 
  float alpha, ElConstMatrix_s ASub );
EL_EXPORT ElError ElMatrixUpdateImagPartOfSubmatrix_z
( ElMatrix_z A, const ElInt* rowInd, const ElInt* colInd, 
  double alpha, ElConstMatrix_d ASub );

/* void Matrix<T>::MakeSubmatrixReal
   ( const std::vector<Int>& I, const std::vector<Int>& J )
   -------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixMakeSubmatrixReal_c
( ElMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInd, ElInt numColInds, const ElInt* colInd );
EL_EXPORT ElError ElMatrixMakeSubmatrixReal_z
( ElMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInd, ElInt numColInds, const ElInt* colInd );

/* void Matrix<T>::ConjugateSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J )
   -------------------------------------------------------------------- */
EL_EXPORT ElError ElMatrixConjugateSubmatrix_c
( ElMatrix_c A, 
  ElInt numRowInds, const ElInt* rowInd, 
  ElInt numColInds, const ElInt* colInd );
EL_EXPORT ElError ElMatrixConjugateSubmatrix_z
( ElMatrix_z A, 
  ElInt numRowInds, const ElInt* rowInd, 
  ElInt numColInds, const ElInt* colInd );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_MATRIX_C_H */
