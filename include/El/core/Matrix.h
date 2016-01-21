/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
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

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_MATRIX_C_H */
