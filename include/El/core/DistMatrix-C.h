/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISTMATRIX_C_H
#define EL_DISTMATRIX_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* NOTE: "ElInt* rowInd" and "ElInt* colInd" have been converted to 
         "ElInt* rowInd" and "ElInt* colInd" due to what appears to be a bug
         in /usr/bin/c++ in Ubuntu */

typedef struct 
{
    ElDist colDist, rowDist;
    ElInt colAlign, rowAlign;
    ElInt root;
    ElConstGrid grid;
} ElDistData;

/* An anonymous struct meant as a placeholder for AbstractDistMatrix<T>
   -------------------------------------------------------------------- */
typedef struct ElDistMatrix_iDummy* ElDistMatrix_i;
typedef struct ElDistMatrix_sDummy* ElDistMatrix_s;
typedef struct ElDistMatrix_dDummy* ElDistMatrix_d;
typedef struct ElDistMatrix_cDummy* ElDistMatrix_c;
typedef struct ElDistMatrix_zDummy* ElDistMatrix_z;

typedef const struct ElDistMatrix_iDummy* ElConstDistMatrix_i;
typedef const struct ElDistMatrix_sDummy* ElConstDistMatrix_s;
typedef const struct ElDistMatrix_dDummy* ElConstDistMatrix_d;
typedef const struct ElDistMatrix_cDummy* ElConstDistMatrix_c;
typedef const struct ElDistMatrix_zDummy* ElConstDistMatrix_z;

/* DistMatrix<T,MC,MR>::DistMatrix( const Grid& g )
   ------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixCreate_i( ElConstGrid g, ElDistMatrix_i* A );
EL_EXPORT ElError ElDistMatrixCreate_s( ElConstGrid g, ElDistMatrix_s* A );
EL_EXPORT ElError ElDistMatrixCreate_d( ElConstGrid g, ElDistMatrix_d* A );
EL_EXPORT ElError ElDistMatrixCreate_c( ElConstGrid g, ElDistMatrix_c* A );
EL_EXPORT ElError ElDistMatrixCreate_z( ElConstGrid g, ElDistMatrix_z* A );

/* DistMatrix<T,U,V>::DistMatrix( const Grid& g )
   ---------------------------------------------- */
EL_EXPORT ElError ElDistMatrixCreateSpecific_i
( ElDist U, ElDist V, ElConstGrid g, ElDistMatrix_i* A );
EL_EXPORT ElError ElDistMatrixCreateSpecific_s
( ElDist U, ElDist V, ElConstGrid g, ElDistMatrix_s* A );
EL_EXPORT ElError ElDistMatrixCreateSpecific_d
( ElDist U, ElDist V, ElConstGrid g, ElDistMatrix_d* A );
EL_EXPORT ElError ElDistMatrixCreateSpecific_c
( ElDist U, ElDist V, ElConstGrid g, ElDistMatrix_c* A );
EL_EXPORT ElError ElDistMatrixCreateSpecific_z
( ElDist U, ElDist V, ElConstGrid g, ElDistMatrix_z* A );

/* DistMatrix<T,U,V>::~DistMatrix()
   -------------------------------- */
EL_EXPORT ElError ElDistMatrixDestroy_i( ElConstDistMatrix_i A );
EL_EXPORT ElError ElDistMatrixDestroy_s( ElConstDistMatrix_s A );
EL_EXPORT ElError ElDistMatrixDestroy_d( ElConstDistMatrix_d A );
EL_EXPORT ElError ElDistMatrixDestroy_c( ElConstDistMatrix_c A );
EL_EXPORT ElError ElDistMatrixDestroy_z( ElConstDistMatrix_z A );

/* void DistMatrix<T,U,V>::Empty()
   ------------------------------- */
EL_EXPORT ElError ElDistMatrixEmpty_i( ElDistMatrix_i A );
EL_EXPORT ElError ElDistMatrixEmpty_s( ElDistMatrix_s A );
EL_EXPORT ElError ElDistMatrixEmpty_d( ElDistMatrix_d A );
EL_EXPORT ElError ElDistMatrixEmpty_c( ElDistMatrix_c A );
EL_EXPORT ElError ElDistMatrixEmpty_z( ElDistMatrix_z A );

/* void DistMatrix<T,U,V>::EmptyData()
   ----------------------------------- */
EL_EXPORT ElError ElDistMatrixEmptyData_i( ElDistMatrix_i A );
EL_EXPORT ElError ElDistMatrixEmptyData_s( ElDistMatrix_s A );
EL_EXPORT ElError ElDistMatrixEmptyData_d( ElDistMatrix_d A );
EL_EXPORT ElError ElDistMatrixEmptyData_c( ElDistMatrix_c A );
EL_EXPORT ElError ElDistMatrixEmptyData_z( ElDistMatrix_z A );

/* void DistMatrix<T,U,V>::SetGrid( const Grid& g )
   ------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixSetGrid_i( ElDistMatrix_i A, ElConstGrid grid );
EL_EXPORT ElError ElDistMatrixSetGrid_s( ElDistMatrix_s A, ElConstGrid grid );
EL_EXPORT ElError ElDistMatrixSetGrid_d( ElDistMatrix_d A, ElConstGrid grid );
EL_EXPORT ElError ElDistMatrixSetGrid_c( ElDistMatrix_c A, ElConstGrid grid );
EL_EXPORT ElError ElDistMatrixSetGrid_z( ElDistMatrix_z A, ElConstGrid grid );

/* void DistMatrix<T,U,V>::Resize( Int height, Int width )
   ------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixResize_i
( ElDistMatrix_i A, ElInt height, ElInt width );
EL_EXPORT ElError ElDistMatrixResize_s
( ElDistMatrix_s A, ElInt height, ElInt width );
EL_EXPORT ElError ElDistMatrixResize_d
( ElDistMatrix_d A, ElInt height, ElInt width );
EL_EXPORT ElError ElDistMatrixResize_c
( ElDistMatrix_c A, ElInt height, ElInt width );
EL_EXPORT ElError ElDistMatrixResize_z
( ElDistMatrix_z A, ElInt height, ElInt width );

/* void DistMatrix<T,U,V>::Resize( Int height, Int width, Int ldim )
   ----------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixResizeWithLDim_i
( ElDistMatrix_i A, ElInt height, ElInt width, ElInt ldim );
EL_EXPORT ElError ElDistMatrixResizeWithLDim_s
( ElDistMatrix_s A, ElInt height, ElInt width, ElInt ldim );
EL_EXPORT ElError ElDistMatrixResizeWithLDim_d
( ElDistMatrix_d A, ElInt height, ElInt width, ElInt ldim );
EL_EXPORT ElError ElDistMatrixResizeWithLDim_c
( ElDistMatrix_c A, ElInt height, ElInt width, ElInt ldim );
EL_EXPORT ElError ElDistMatrixResizeWithLDim_z
( ElDistMatrix_z A, ElInt height, ElInt width, ElInt ldim );

/* void DistMatrix<T,U,V>::MakeConsistent( bool includeViewers )
   ------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixMakeConsistent_i
( ElDistMatrix_i A, bool includeViewers );
EL_EXPORT ElError ElDistMatrixMakeConsistent_s
( ElDistMatrix_s A, bool includeViewers );
EL_EXPORT ElError ElDistMatrixMakeConsistent_d
( ElDistMatrix_d A, bool includeViewers );
EL_EXPORT ElError ElDistMatrixMakeConsistent_c
( ElDistMatrix_c A, bool includeViewers );
EL_EXPORT ElError ElDistMatrixMakeConsistent_z
( ElDistMatrix_z A, bool includeViewers );

/* void DistMatrix<T,U,V>::MakeSizeConsistent( bool includeViewers )
   ----------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixMakeSizeConsistent_i
( ElDistMatrix_i A, bool includeViewers );
EL_EXPORT ElError ElDistMatrixMakeSizeConsistent_s
( ElDistMatrix_s A, bool includeViewers );
EL_EXPORT ElError ElDistMatrixMakeSizeConsistent_d
( ElDistMatrix_d A, bool includeViewers );
EL_EXPORT ElError ElDistMatrixMakeSizeConsistent_c
( ElDistMatrix_c A, bool includeViewers );
EL_EXPORT ElError ElDistMatrixMakeSizeConsistent_z
( ElDistMatrix_z A, bool includeViewers );

/* void DistMatrix<T,U,V>::Align( Int colAlign, Int rowAlign, bool constrain )
   --------------------------------------------------------------------------- 
*/
EL_EXPORT ElError ElDistMatrixAlign_i
( ElDistMatrix_i A, ElInt colAlign, ElInt rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlign_s
( ElDistMatrix_s A, ElInt colAlign, ElInt rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlign_d
( ElDistMatrix_d A, ElInt colAlign, ElInt rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlign_c
( ElDistMatrix_c A, ElInt colAlign, ElInt rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlign_z
( ElDistMatrix_z A, ElInt colAlign, ElInt rowAlign, bool constrain );

/* void DistMatrix<T,U,V>::AlignCols( Int colAlign, bool constrain )
   ----------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixAlignCols_i
( ElDistMatrix_i A, ElInt colAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignCols_s
( ElDistMatrix_s A, ElInt colAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignCols_d
( ElDistMatrix_d A, ElInt colAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignCols_c
( ElDistMatrix_c A, ElInt colAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignCols_z
( ElDistMatrix_z A, ElInt colAlign, bool constrain );

/* void DistMatrix<T,U,V>::AlignRows( Int rowAlign, bool constrain )
   ----------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixAlignRows_i
( ElDistMatrix_i A, ElInt rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRows_s
( ElDistMatrix_s A, ElInt rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRows_d
( ElDistMatrix_d A, ElInt rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRows_c
( ElDistMatrix_c A, ElInt rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRows_z
( ElDistMatrix_z A, ElInt rowAlign, bool constrain );

/* void DistMatrix<T,U,V>::FreeAlignments()
   ---------------------------------------- */
EL_EXPORT ElError ElDistMatrixFreeAlignments_i( ElDistMatrix_i A );
EL_EXPORT ElError ElDistMatrixFreeAlignments_s( ElDistMatrix_s A );
EL_EXPORT ElError ElDistMatrixFreeAlignments_d( ElDistMatrix_d A );
EL_EXPORT ElError ElDistMatrixFreeAlignments_c( ElDistMatrix_c A );
EL_EXPORT ElError ElDistMatrixFreeAlignments_z( ElDistMatrix_z A );

/* void DistMatrix<T,U,V>::SetRoot( Int root )
   ------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSetRoot_i( ElDistMatrix_i A, ElInt root );
EL_EXPORT ElError ElDistMatrixSetRoot_s( ElDistMatrix_s A, ElInt root );
EL_EXPORT ElError ElDistMatrixSetRoot_d( ElDistMatrix_d A, ElInt root );
EL_EXPORT ElError ElDistMatrixSetRoot_c( ElDistMatrix_c A, ElInt root );
EL_EXPORT ElError ElDistMatrixSetRoot_z( ElDistMatrix_z A, ElInt root );

/* void DistMatrix<T,U,V>::AlignWith( const DistData& data, bool constrain )
   ------------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixAlignWith_i
( ElDistMatrix_i A, ElDistData distData, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignWith_s
( ElDistMatrix_s A, ElDistData distData, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignWith_d
( ElDistMatrix_d A, ElDistData distData, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignWith_c
( ElDistMatrix_c A, ElDistData distData, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignWith_z
( ElDistMatrix_z A, ElDistData distData, bool constrain );

/* void DistMatrix<T,U,V>::AlignColsWith( const DistData& data, bool constrain )
   -----------------------------------------------------------------------------
*/
EL_EXPORT ElError ElDistMatrixAlignColsWith_i
( ElDistMatrix_i A, ElDistData distData, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignColsWith_s
( ElDistMatrix_s A, ElDistData distData, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignColsWith_d
( ElDistMatrix_d A, ElDistData distData, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignColsWith_c
( ElDistMatrix_c A, ElDistData distData, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignColsWith_z
( ElDistMatrix_z A, ElDistData distData, bool constrain );

/* void DistMatrix<T,U,V>::AlignRowsWith( const DistData& data, bool constrain )
   -----------------------------------------------------------------------------
*/
EL_EXPORT ElError ElDistMatrixAlignRowsWith_i
( ElDistMatrix_i A, ElDistData distData, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRowsWith_s
( ElDistMatrix_s A, ElDistData distData, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRowsWith_d
( ElDistMatrix_d A, ElDistData distData, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRowsWith_c
( ElDistMatrix_c A, ElDistData distData, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRowsWith_z
( ElDistMatrix_z A, ElDistData distData, bool constrain );

/* void DistMatrix<T,U,V>::AlignAndResize
  ( Int colAlign, Int rowAlign, Int height, Int width, 
    bool force, bool constrain )
   --------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixAlignAndResize_i
( ElDistMatrix_i A, ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignAndResize_s
( ElDistMatrix_s A, ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignAndResize_d
( ElDistMatrix_d A, ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignAndResize_c
( ElDistMatrix_c A, ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignAndResize_z
( ElDistMatrix_z A, ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );

/* void DistMatrix<T,U,V>::AlignColsAndResize
   ( Int colAlign, Int height, Int width, bool force, bool constrain )
   ------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixAlignColsAndResize_i
( ElDistMatrix_i A, ElInt colAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignColsAndResize_s
( ElDistMatrix_s A, ElInt colAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignColsAndResize_d
( ElDistMatrix_d A, ElInt colAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignColsAndResize_c
( ElDistMatrix_c A, ElInt colAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignColsAndResize_z
( ElDistMatrix_z A, ElInt colAlign, ElInt height, ElInt width, 
  bool force, bool constrain );

/* void DistMatrix<T,U,V>::AlignRowsAndResize
  ( Int rowAlign, Int height, Int width, bool force, bool constrain )
   ------------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixAlignRowsAndResize_i
( ElDistMatrix_i A, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRowsAndResize_s
( ElDistMatrix_s A, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRowsAndResize_d
( ElDistMatrix_d A, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRowsAndResize_c
( ElDistMatrix_c A, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRowsAndResize_z
( ElDistMatrix_z A, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );

/* void DistMatrix<T,U,V>::Attach
   ( Int height, Int width, const Grid& grid, Int colAlign, Int rowAlign, 
     T* buffer, Int ldim, Int root )
   ---------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixAttach_i
( ElDistMatrix_i A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, ElInt* buffer, ElInt ldim, 
  ElInt root );
EL_EXPORT ElError ElDistMatrixAttach_s
( ElDistMatrix_s A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, float* buffer, ElInt ldim, 
  ElInt root );
EL_EXPORT ElError ElDistMatrixAttach_d
( ElDistMatrix_d A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, double* buffer, ElInt ldim, 
  ElInt root );
EL_EXPORT ElError ElDistMatrixAttach_c
( ElDistMatrix_c A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, complex_float* buffer, ElInt ldim, 
  ElInt root );
EL_EXPORT ElError ElDistMatrixAttach_z
( ElDistMatrix_z A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, complex_double* buffer, ElInt ldim, 
  ElInt root );

/* void DistMatrix<T,U,V>::LockedAttach
   ( Int height, Int width, const Grid& grid, Int colAlign, Int rowAlign, 
     const T* buffer, Int ldim, Int root )
   ---------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixLockedAttach_i
( ElDistMatrix_i A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, const ElInt* buffer, 
  ElInt ldim, ElInt root );
EL_EXPORT ElError ElDistMatrixLockedAttach_s
( ElDistMatrix_s A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, const float* buffer, 
  ElInt ldim, ElInt root );
EL_EXPORT ElError ElDistMatrixLockedAttach_d
( ElDistMatrix_d A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, const double* buffer, 
  ElInt ldim, ElInt root );
EL_EXPORT ElError ElDistMatrixLockedAttach_c
( ElDistMatrix_c A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, const complex_float* buffer, 
  ElInt ldim, ElInt root );
EL_EXPORT ElError ElDistMatrixLockedAttach_z
( ElDistMatrix_z A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, const complex_double* buffer, 
  ElInt ldim, ElInt root );

/* Int DistMatrix<T,U,V>::Height() const
   ------------------------------------- */
EL_EXPORT ElError ElDistMatrixHeight_i( ElConstDistMatrix_i A, ElInt* height );
EL_EXPORT ElError ElDistMatrixHeight_s( ElConstDistMatrix_s A, ElInt* height );
EL_EXPORT ElError ElDistMatrixHeight_d( ElConstDistMatrix_d A, ElInt* height );
EL_EXPORT ElError ElDistMatrixHeight_c( ElConstDistMatrix_c A, ElInt* height );
EL_EXPORT ElError ElDistMatrixHeight_z( ElConstDistMatrix_z A, ElInt* height );

/* Int DistMatrix<T,U,V>::Width() const
   ------------------------------------ */
EL_EXPORT ElError ElDistMatrixWidth_i( ElConstDistMatrix_i A, ElInt* width );
EL_EXPORT ElError ElDistMatrixWidth_s( ElConstDistMatrix_s A, ElInt* width );
EL_EXPORT ElError ElDistMatrixWidth_d( ElConstDistMatrix_d A, ElInt* width );
EL_EXPORT ElError ElDistMatrixWidth_c( ElConstDistMatrix_c A, ElInt* width );
EL_EXPORT ElError ElDistMatrixWidth_z( ElConstDistMatrix_z A, ElInt* width );

/* Int DistMatrix<T,U,V>::DiagonalLength( Int offset ) const
   --------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixDiagonalLength_i
( ElConstDistMatrix_i A, ElInt offset, ElInt* length );
EL_EXPORT ElError ElDistMatrixDiagonalLength_s
( ElConstDistMatrix_s A, ElInt offset, ElInt* length );
EL_EXPORT ElError ElDistMatrixDiagonalLength_d
( ElConstDistMatrix_d A, ElInt offset, ElInt* length );
EL_EXPORT ElError ElDistMatrixDiagonalLength_c
( ElConstDistMatrix_c A, ElInt offset, ElInt* length );
EL_EXPORT ElError ElDistMatrixDiagonalLength_z
( ElConstDistMatrix_z A, ElInt offset, ElInt* length );

/* bool DistMatrix<T,U,V>::Viewing() const
   --------------------------------------- */
EL_EXPORT ElError ElDistMatrixViewing_i( ElConstDistMatrix_i A, bool* viewing );
EL_EXPORT ElError ElDistMatrixViewing_s( ElConstDistMatrix_s A, bool* viewing );
EL_EXPORT ElError ElDistMatrixViewing_d( ElConstDistMatrix_d A, bool* viewing );
EL_EXPORT ElError ElDistMatrixViewing_c( ElConstDistMatrix_c A, bool* viewing );
EL_EXPORT ElError ElDistMatrixViewing_z( ElConstDistMatrix_z A, bool* viewing );

/* bool DistMatrix<T,U,V>::Locked() const
   --------------------------------------- */
EL_EXPORT ElError ElDistMatrixLocked_i( ElConstDistMatrix_i A, bool* locked );
EL_EXPORT ElError ElDistMatrixLocked_s( ElConstDistMatrix_s A, bool* locked );
EL_EXPORT ElError ElDistMatrixLocked_d( ElConstDistMatrix_d A, bool* locked );
EL_EXPORT ElError ElDistMatrixLocked_c( ElConstDistMatrix_c A, bool* locked );
EL_EXPORT ElError ElDistMatrixLocked_z( ElConstDistMatrix_z A, bool* locked );

/* Int DistMatrix<T,U,V>::LocalHeight() const
   ------------------------------------------ */
EL_EXPORT ElError ElDistMatrixLocalHeight_i
( ElConstDistMatrix_i A, ElInt* localHeight );
EL_EXPORT ElError ElDistMatrixLocalHeight_s
( ElConstDistMatrix_s A, ElInt* localHeight );
EL_EXPORT ElError ElDistMatrixLocalHeight_d
( ElConstDistMatrix_d A, ElInt* localHeight );
EL_EXPORT ElError ElDistMatrixLocalHeight_c
( ElConstDistMatrix_c A, ElInt* localHeight );
EL_EXPORT ElError ElDistMatrixLocalHeight_z
( ElConstDistMatrix_z A, ElInt* localHeight );

/* Int DistMatrix<T,U,V>::LocalWidth() const
   ----------------------------------------- */
EL_EXPORT ElError ElDistMatrixLocalWidth_i
( ElConstDistMatrix_i A, ElInt* localWidth );
EL_EXPORT ElError ElDistMatrixLocalWidth_s
( ElConstDistMatrix_s A, ElInt* localWidth );
EL_EXPORT ElError ElDistMatrixLocalWidth_d
( ElConstDistMatrix_d A, ElInt* localWidth );
EL_EXPORT ElError ElDistMatrixLocalWidth_c
( ElConstDistMatrix_c A, ElInt* localWidth );
EL_EXPORT ElError ElDistMatrixLocalWidth_z
( ElConstDistMatrix_z A, ElInt* localWidth );

/* Int DistMatrix<T,U,V>::LDim() const
   ----------------------------------- */
EL_EXPORT ElError ElDistMatrixLDim_i( ElConstDistMatrix_i A, ElInt* ldim );
EL_EXPORT ElError ElDistMatrixLDim_s( ElConstDistMatrix_s A, ElInt* ldim );
EL_EXPORT ElError ElDistMatrixLDim_d( ElConstDistMatrix_d A, ElInt* ldim );
EL_EXPORT ElError ElDistMatrixLDim_c( ElConstDistMatrix_c A, ElInt* ldim );
EL_EXPORT ElError ElDistMatrixLDim_z( ElConstDistMatrix_z A, ElInt* ldim );

/* Matrix<T>& DistMatrix<T,U,V>::Matrix() 
   -------------------------------------- */
EL_EXPORT ElError ElDistMatrixMatrix_i( ElDistMatrix_i A, ElMatrix_i* ALoc );
EL_EXPORT ElError ElDistMatrixMatrix_s( ElDistMatrix_s A, ElMatrix_s* ALoc );
EL_EXPORT ElError ElDistMatrixMatrix_d( ElDistMatrix_d A, ElMatrix_d* ALoc );
EL_EXPORT ElError ElDistMatrixMatrix_c( ElDistMatrix_c A, ElMatrix_c* ALoc );
EL_EXPORT ElError ElDistMatrixMatrix_z( ElDistMatrix_z A, ElMatrix_z* ALoc );

/* const Matrix<T>& DistMatrix<T,U,V>::LockedMatrix() const
   -------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixLockedMatrix_i
( ElConstDistMatrix_i A, ElConstMatrix_i* ALoc );
EL_EXPORT ElError ElDistMatrixLockedMatrix_s
( ElConstDistMatrix_s A, ElConstMatrix_s* ALoc );
EL_EXPORT ElError ElDistMatrixLockedMatrix_d
( ElConstDistMatrix_d A, ElConstMatrix_d* ALoc );
EL_EXPORT ElError ElDistMatrixLockedMatrix_c
( ElConstDistMatrix_c A, ElConstMatrix_c* ALoc );
EL_EXPORT ElError ElDistMatrixLockedMatrix_z
( ElConstDistMatrix_z A, ElConstMatrix_z* ALoc );

/* size_t DistMatrix<T,U,V>::AllocatedMemory() const
   ------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixAllocatedMemory_i
( ElConstDistMatrix_i A, size_t* mem );
EL_EXPORT ElError ElDistMatrixAllocatedMemory_s
( ElConstDistMatrix_s A, size_t* mem );
EL_EXPORT ElError ElDistMatrixAllocatedMemory_d
( ElConstDistMatrix_d A, size_t* mem );
EL_EXPORT ElError ElDistMatrixAllocatedMemory_c
( ElConstDistMatrix_c A, size_t* mem );
EL_EXPORT ElError ElDistMatrixAllocatedMemory_z
( ElConstDistMatrix_z A, size_t* mem );

/* T* DistMatrix<T,U,V>::Buffer()
   ------------------------------ */
EL_EXPORT ElError ElDistMatrixBuffer_i
( ElDistMatrix_i A, ElInt** buffer );
EL_EXPORT ElError ElDistMatrixBuffer_s
( ElDistMatrix_s A, float** buffer );
EL_EXPORT ElError ElDistMatrixBuffer_d
( ElDistMatrix_d A, double** buffer );
EL_EXPORT ElError ElDistMatrixBuffer_c
( ElDistMatrix_c A, complex_float** buffer );
EL_EXPORT ElError ElDistMatrixBuffer_z
( ElDistMatrix_z A, complex_double** buffer );

/* const T* DistMatrix<T,U,V>::LockedBuffer() const
   ------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixLockedBuffer_i
( ElConstDistMatrix_i A, const ElInt** buffer );
EL_EXPORT ElError ElDistMatrixLockedBuffer_s
( ElConstDistMatrix_s A, const float** buffer );
EL_EXPORT ElError ElDistMatrixLockedBuffer_d
( ElConstDistMatrix_d A, const double** buffer );
EL_EXPORT ElError ElDistMatrixLockedBuffer_c
( ElConstDistMatrix_c A, const complex_float** buffer );
EL_EXPORT ElError ElDistMatrixLockedBuffer_z
( ElConstDistMatrix_z A, const complex_double** buffer );

/* const Grid& DistMatrix<T,U,V>::Grid() const
   ------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGrid_i
( ElConstDistMatrix_i A, ElConstGrid* grid );
EL_EXPORT ElError ElDistMatrixGrid_s
( ElConstDistMatrix_s A, ElConstGrid* grid );
EL_EXPORT ElError ElDistMatrixGrid_d
( ElConstDistMatrix_d A, ElConstGrid* grid );
EL_EXPORT ElError ElDistMatrixGrid_c
( ElConstDistMatrix_c A, ElConstGrid* grid );
EL_EXPORT ElError ElDistMatrixGrid_z
( ElConstDistMatrix_z A, ElConstGrid* grid );

/* bool DistMatrix<T,U,V>::ColConstrained() const
   ---------------------------------------------- */
EL_EXPORT ElError ElDistMatrixColConstrained_i
( ElConstDistMatrix_i A, bool* colConst );
EL_EXPORT ElError ElDistMatrixColConstrained_s
( ElConstDistMatrix_s A, bool* colConst );
EL_EXPORT ElError ElDistMatrixColConstrained_d
( ElConstDistMatrix_d A, bool* colConst );
EL_EXPORT ElError ElDistMatrixColConstrained_c
( ElConstDistMatrix_c A, bool* colConst );
EL_EXPORT ElError ElDistMatrixColConstrained_z
( ElConstDistMatrix_z A, bool* colConst );

/* bool DistMatrix<T,U,V>::RowConstrained() const
   ---------------------------------------------- */
EL_EXPORT ElError ElDistMatrixRowConstrained_i
( ElConstDistMatrix_i A, bool* rowConst );
EL_EXPORT ElError ElDistMatrixRowConstrained_s
( ElConstDistMatrix_s A, bool* rowConst );
EL_EXPORT ElError ElDistMatrixRowConstrained_d
( ElConstDistMatrix_d A, bool* rowConst );
EL_EXPORT ElError ElDistMatrixRowConstrained_c
( ElConstDistMatrix_c A, bool* rowConst );
EL_EXPORT ElError ElDistMatrixRowConstrained_z
( ElConstDistMatrix_z A, bool* rowConst );

/* bool DistMatrix<T,U,V>::RootConstrained() const
   ----------------------------------------------- */
EL_EXPORT ElError ElDistMatrixRootConstrained_i
( ElConstDistMatrix_i A, bool* rootConst );
EL_EXPORT ElError ElDistMatrixRootConstrained_s
( ElConstDistMatrix_s A, bool* rootConst );
EL_EXPORT ElError ElDistMatrixRootConstrained_d
( ElConstDistMatrix_d A, bool* rootConst );
EL_EXPORT ElError ElDistMatrixRootConstrained_c
( ElConstDistMatrix_c A, bool* rootConst );
EL_EXPORT ElError ElDistMatrixRootConstrained_z
( ElConstDistMatrix_z A, bool* rootConst );

/* Int DistMatrix<T,U,V>::ColAlign() const
   --------------------------------------- */
EL_EXPORT ElError ElDistMatrixColAlign_i
( ElConstDistMatrix_i A, ElInt* colAlign );
EL_EXPORT ElError ElDistMatrixColAlign_s
( ElConstDistMatrix_s A, ElInt* colAlign );
EL_EXPORT ElError ElDistMatrixColAlign_d
( ElConstDistMatrix_d A, ElInt* colAlign );
EL_EXPORT ElError ElDistMatrixColAlign_c
( ElConstDistMatrix_c A, ElInt* colAlign );
EL_EXPORT ElError ElDistMatrixColAlign_z
( ElConstDistMatrix_z A, ElInt* colAlign );

/* Int DistMatrix<T,U,V>::RowAlign() const
   --------------------------------------- */
EL_EXPORT ElError ElDistMatrixRowAlign_i
( ElConstDistMatrix_i A, ElInt* rowAlign );
EL_EXPORT ElError ElDistMatrixRowAlign_s
( ElConstDistMatrix_s A, ElInt* rowAlign );
EL_EXPORT ElError ElDistMatrixRowAlign_d
( ElConstDistMatrix_d A, ElInt* rowAlign );
EL_EXPORT ElError ElDistMatrixRowAlign_c
( ElConstDistMatrix_c A, ElInt* rowAlign );
EL_EXPORT ElError ElDistMatrixRowAlign_z
( ElConstDistMatrix_z A, ElInt* rowAlign );

/* Int DistMatrix<T,U,V>::ColShift() const
   --------------------------------------- */
EL_EXPORT ElError ElDistMatrixColShift_i
( ElConstDistMatrix_i A, ElInt* colShift );
EL_EXPORT ElError ElDistMatrixColShift_s
( ElConstDistMatrix_s A, ElInt* colShift );
EL_EXPORT ElError ElDistMatrixColShift_d
( ElConstDistMatrix_d A, ElInt* colShift );
EL_EXPORT ElError ElDistMatrixColShift_c
( ElConstDistMatrix_c A, ElInt* colShift );
EL_EXPORT ElError ElDistMatrixColShift_z
( ElConstDistMatrix_z A, ElInt* colShift );

/* Int DistMatrix<T,U,V>::RowShift() const
   --------------------------------------- */
EL_EXPORT ElError ElDistMatrixRowShift_i
( ElConstDistMatrix_i A, ElInt* rowShift );
EL_EXPORT ElError ElDistMatrixRowShift_s
( ElConstDistMatrix_s A, ElInt* rowShift );
EL_EXPORT ElError ElDistMatrixRowShift_d
( ElConstDistMatrix_d A, ElInt* rowShift );
EL_EXPORT ElError ElDistMatrixRowShift_c
( ElConstDistMatrix_c A, ElInt* rowShift );
EL_EXPORT ElError ElDistMatrixRowShift_z
( ElConstDistMatrix_z A, ElInt* rowShift );

/* Int DistMatrix<T,U,V>::ColRank() const
   -------------------------------------- */
EL_EXPORT ElError ElDistMatrixColRank_i( ElConstDistMatrix_i A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixColRank_s( ElConstDistMatrix_s A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixColRank_d( ElConstDistMatrix_d A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixColRank_c( ElConstDistMatrix_c A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixColRank_z( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::RowRank() const
   -------------------------------------- */
EL_EXPORT ElError ElDistMatrixRowRank_i( ElConstDistMatrix_i A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixRowRank_s( ElConstDistMatrix_s A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixRowRank_d( ElConstDistMatrix_d A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixRowRank_c( ElConstDistMatrix_c A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixRowRank_z( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::PartialColRank() const
   --------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialColRank_i
( ElConstDistMatrix_i A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialColRank_s
( ElConstDistMatrix_s A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialColRank_d
( ElConstDistMatrix_d A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialColRank_c
( ElConstDistMatrix_c A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialColRank_z
( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::PartialRowRank() const
   --------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialRowRank_i
( ElConstDistMatrix_i A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialRowRank_s
( ElConstDistMatrix_s A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialRowRank_d
( ElConstDistMatrix_d A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialRowRank_c
( ElConstDistMatrix_c A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialRowRank_z
( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::PartialUnionColRank() const
   -------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialUnionColRank_i
( ElConstDistMatrix_i A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionColRank_s
( ElConstDistMatrix_s A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionColRank_d
( ElConstDistMatrix_d A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionColRank_c
( ElConstDistMatrix_c A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionColRank_z
( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::PartialUnionRowRank() const
   -------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialUnionRowRank_i
( ElConstDistMatrix_i A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionRowRank_s
( ElConstDistMatrix_s A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionRowRank_d
( ElConstDistMatrix_d A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionRowRank_c
( ElConstDistMatrix_c A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionRowRank_z
( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::DistRank() const
   --------------------------------------- */
EL_EXPORT ElError ElDistMatrixDistRank_i
( ElConstDistMatrix_i A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixDistRank_s
( ElConstDistMatrix_s A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixDistRank_d
( ElConstDistMatrix_d A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixDistRank_c
( ElConstDistMatrix_c A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixDistRank_z
( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::CrossRank() const
   ---------------------------------------- */
EL_EXPORT ElError ElDistMatrixCrossRank_i( ElConstDistMatrix_i A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixCrossRank_s( ElConstDistMatrix_s A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixCrossRank_d( ElConstDistMatrix_d A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixCrossRank_c( ElConstDistMatrix_c A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixCrossRank_z( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::RedundantRank() const
   -------------------------------------------- */
EL_EXPORT ElError ElDistMatrixRedundantRank_i
( ElConstDistMatrix_i A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixRedundantRank_s
( ElConstDistMatrix_s A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixRedundantRank_d
( ElConstDistMatrix_d A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixRedundantRank_c
( ElConstDistMatrix_c A, ElInt* rank );
EL_EXPORT ElError ElDistMatrixRedundantRank_z
( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::Root() const
   ----------------------------------- */
EL_EXPORT ElError ElDistMatrixRoot_i( ElConstDistMatrix_i A, ElInt* root );
EL_EXPORT ElError ElDistMatrixRoot_s( ElConstDistMatrix_s A, ElInt* root );
EL_EXPORT ElError ElDistMatrixRoot_d( ElConstDistMatrix_d A, ElInt* root );
EL_EXPORT ElError ElDistMatrixRoot_c( ElConstDistMatrix_c A, ElInt* root );
EL_EXPORT ElError ElDistMatrixRoot_z( ElConstDistMatrix_z A, ElInt* root );

/* bool DistMatrix<T,U,V>::Participating() const
   --------------------------------------------- */
EL_EXPORT ElError ElDistMatrixParticipating_i
( ElConstDistMatrix_i A, bool* participating );
EL_EXPORT ElError ElDistMatrixParticipating_s
( ElConstDistMatrix_s A, bool* participating );
EL_EXPORT ElError ElDistMatrixParticipating_d
( ElConstDistMatrix_d A, bool* participating );
EL_EXPORT ElError ElDistMatrixParticipating_c
( ElConstDistMatrix_c A, bool* participating );
EL_EXPORT ElError ElDistMatrixParticipating_z
( ElConstDistMatrix_z A, bool* participating );

/* Int DistMatrix<T,U,V>::RowOwner( Int i ) const
   ---------------------------------------------- */
EL_EXPORT ElError ElDistMatrixRowOwner_i
( ElConstDistMatrix_i A, ElInt i, ElInt* rowOwner );
EL_EXPORT ElError ElDistMatrixRowOwner_s
( ElConstDistMatrix_s A, ElInt i, ElInt* rowOwner );
EL_EXPORT ElError ElDistMatrixRowOwner_d
( ElConstDistMatrix_d A, ElInt i, ElInt* rowOwner );
EL_EXPORT ElError ElDistMatrixRowOwner_c
( ElConstDistMatrix_c A, ElInt i, ElInt* rowOwner );
EL_EXPORT ElError ElDistMatrixRowOwner_z
( ElConstDistMatrix_z A, ElInt i, ElInt* rowOwner );

/* Int DistMatrix<T,U,V>::ColOwner( Int j ) const
   ---------------------------------------------- */
EL_EXPORT ElError ElDistMatrixColOwner_i
( ElConstDistMatrix_i A, ElInt j, ElInt* colOwner );
EL_EXPORT ElError ElDistMatrixColOwner_s
( ElConstDistMatrix_s A, ElInt j, ElInt* colOwner );
EL_EXPORT ElError ElDistMatrixColOwner_d
( ElConstDistMatrix_d A, ElInt j, ElInt* colOwner );
EL_EXPORT ElError ElDistMatrixColOwner_c
( ElConstDistMatrix_c A, ElInt j, ElInt* colOwner );
EL_EXPORT ElError ElDistMatrixColOwner_z
( ElConstDistMatrix_z A, ElInt j, ElInt* colOwner );

/* Int DistMatrix<T,U,V>::Owner( Int i, Int j ) const
   -------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixOwner_i
( ElConstDistMatrix_i A, ElInt i, ElInt j, ElInt* owner );
EL_EXPORT ElError ElDistMatrixOwner_s
( ElConstDistMatrix_s A, ElInt i, ElInt j, ElInt* owner );
EL_EXPORT ElError ElDistMatrixOwner_d
( ElConstDistMatrix_d A, ElInt i, ElInt j, ElInt* owner );
EL_EXPORT ElError ElDistMatrixOwner_c
( ElConstDistMatrix_c A, ElInt i, ElInt j, ElInt* owner );
EL_EXPORT ElError ElDistMatrixOwner_z
( ElConstDistMatrix_z A, ElInt i, ElInt j, ElInt* owner );

/* Int DistMatrix<T,U,V>::LocalRow( Int i ) const
   ---------------------------------------------- */
EL_EXPORT ElError ElDistMatrixLocalRow_i
( ElConstDistMatrix_i A, ElInt i, ElInt* iLoc );
EL_EXPORT ElError ElDistMatrixLocalRow_s
( ElConstDistMatrix_s A, ElInt i, ElInt* iLoc );
EL_EXPORT ElError ElDistMatrixLocalRow_d
( ElConstDistMatrix_d A, ElInt i, ElInt* iLoc );
EL_EXPORT ElError ElDistMatrixLocalRow_c
( ElConstDistMatrix_c A, ElInt i, ElInt* iLoc );
EL_EXPORT ElError ElDistMatrixLocalRow_z
( ElConstDistMatrix_z A, ElInt i, ElInt* iLoc );

/* Int DistMatrix<T,U,V>::LocalCol( Int j ) const
   ---------------------------------------------- */
EL_EXPORT ElError ElDistMatrixLocalCol_i
( ElConstDistMatrix_i A, ElInt j, ElInt* jLoc );
EL_EXPORT ElError ElDistMatrixLocalCol_s
( ElConstDistMatrix_s A, ElInt j, ElInt* jLoc );
EL_EXPORT ElError ElDistMatrixLocalCol_d
( ElConstDistMatrix_d A, ElInt j, ElInt* jLoc );
EL_EXPORT ElError ElDistMatrixLocalCol_c
( ElConstDistMatrix_c A, ElInt j, ElInt* jLoc );
EL_EXPORT ElError ElDistMatrixLocalCol_z
( ElConstDistMatrix_z A, ElInt j, ElInt* jLoc );

/* Int DistMatrix<T,U,V>::LocalRowOffset( Int i ) const
   ---------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixLocalRowOffset_i
( ElConstDistMatrix_i A, ElInt i, ElInt* iLoc );
EL_EXPORT ElError ElDistMatrixLocalRowOffset_s
( ElConstDistMatrix_s A, ElInt i, ElInt* iLoc );
EL_EXPORT ElError ElDistMatrixLocalRowOffset_d
( ElConstDistMatrix_d A, ElInt i, ElInt* iLoc );
EL_EXPORT ElError ElDistMatrixLocalRowOffset_c
( ElConstDistMatrix_c A, ElInt i, ElInt* iLoc );
EL_EXPORT ElError ElDistMatrixLocalRowOffset_z
( ElConstDistMatrix_z A, ElInt i, ElInt* iLoc );

/* Int DistMatrix<T,U,V>::LocalColOffset( Int j ) const
   ---------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixLocalColOffset_i
( ElConstDistMatrix_i A, ElInt j, ElInt* jLoc );
EL_EXPORT ElError ElDistMatrixLocalColOffset_s
( ElConstDistMatrix_s A, ElInt j, ElInt* jLoc );
EL_EXPORT ElError ElDistMatrixLocalColOffset_d
( ElConstDistMatrix_d A, ElInt j, ElInt* jLoc );
EL_EXPORT ElError ElDistMatrixLocalColOffset_c
( ElConstDistMatrix_c A, ElInt j, ElInt* jLoc );
EL_EXPORT ElError ElDistMatrixLocalColOffset_z
( ElConstDistMatrix_z A, ElInt j, ElInt* jLoc );

/* Int DistMatrix<T,U,V>::GlobalRow( Int iLoc ) const
   -------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGlobalRow_i
( ElConstDistMatrix_i A, ElInt iLoc, ElInt* i );
EL_EXPORT ElError ElDistMatrixGlobalRow_s
( ElConstDistMatrix_s A, ElInt iLoc, ElInt* i );
EL_EXPORT ElError ElDistMatrixGlobalRow_d
( ElConstDistMatrix_d A, ElInt iLoc, ElInt* i );
EL_EXPORT ElError ElDistMatrixGlobalRow_c
( ElConstDistMatrix_c A, ElInt iLoc, ElInt* i );
EL_EXPORT ElError ElDistMatrixGlobalRow_z
( ElConstDistMatrix_z A, ElInt iLoc, ElInt* i );

/* Int DistMatrix<T,U,V>::GlobalCol( Int jLoc ) const
   -------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGlobalCol_i
( ElConstDistMatrix_i A, ElInt jLoc, ElInt* j );
EL_EXPORT ElError ElDistMatrixGlobalCol_s
( ElConstDistMatrix_s A, ElInt jLoc, ElInt* j );
EL_EXPORT ElError ElDistMatrixGlobalCol_d
( ElConstDistMatrix_d A, ElInt jLoc, ElInt* j );
EL_EXPORT ElError ElDistMatrixGlobalCol_c
( ElConstDistMatrix_c A, ElInt jLoc, ElInt* j );
EL_EXPORT ElError ElDistMatrixGlobalCol_z
( ElConstDistMatrix_z A, ElInt jLoc, ElInt* j );

/* bool DistMatrix<T,U,V>::IsLocalRow( Int i ) const
   ------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixIsLocalRow_i
( ElConstDistMatrix_i A, ElInt i, bool* isLocal );
EL_EXPORT ElError ElDistMatrixIsLocalRow_s
( ElConstDistMatrix_s A, ElInt i, bool* isLocal );
EL_EXPORT ElError ElDistMatrixIsLocalRow_d
( ElConstDistMatrix_d A, ElInt i, bool* isLocal );
EL_EXPORT ElError ElDistMatrixIsLocalRow_c
( ElConstDistMatrix_c A, ElInt i, bool* isLocal );
EL_EXPORT ElError ElDistMatrixIsLocalRow_z
( ElConstDistMatrix_z A, ElInt i, bool* isLocal );

/* bool DistMatrix<T,U,V>::IsLocalCol( Int j ) const
   ------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixIsLocalCol_i
( ElConstDistMatrix_i A, ElInt j, bool* isLocal );
EL_EXPORT ElError ElDistMatrixIsLocalCol_s
( ElConstDistMatrix_s A, ElInt j, bool* isLocal );
EL_EXPORT ElError ElDistMatrixIsLocalCol_d
( ElConstDistMatrix_d A, ElInt j, bool* isLocal );
EL_EXPORT ElError ElDistMatrixIsLocalCol_c
( ElConstDistMatrix_c A, ElInt j, bool* isLocal );
EL_EXPORT ElError ElDistMatrixIsLocalCol_z
( ElConstDistMatrix_z A, ElInt j, bool* isLocal );

/* bool DistMatrix<T,U,V>::IsLocal( Int i, Int j ) const
   ----------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixIsLocal_i
( ElConstDistMatrix_i A, ElInt i, ElInt j, bool* isLocal );
EL_EXPORT ElError ElDistMatrixIsLocal_s
( ElConstDistMatrix_s A, ElInt i, ElInt j, bool* isLocal );
EL_EXPORT ElError ElDistMatrixIsLocal_d
( ElConstDistMatrix_d A, ElInt i, ElInt j, bool* isLocal );
EL_EXPORT ElError ElDistMatrixIsLocal_c
( ElConstDistMatrix_c A, ElInt i, ElInt j, bool* isLocal );
EL_EXPORT ElError ElDistMatrixIsLocal_z
( ElConstDistMatrix_z A, ElInt i, ElInt j, bool* isLocal );

/* DistData DistMatrix<T,U,V>::DistData() const
   -------------------------------------------- */
EL_EXPORT ElError ElDistMatrixDistData_i
( ElConstDistMatrix_i A, ElDistData* distData );
EL_EXPORT ElError ElDistMatrixDistData_s
( ElConstDistMatrix_s A, ElDistData* distData );
EL_EXPORT ElError ElDistMatrixDistData_d
( ElConstDistMatrix_d A, ElDistData* distData );
EL_EXPORT ElError ElDistMatrixDistData_c
( ElConstDistMatrix_c A, ElDistData* distData );
EL_EXPORT ElError ElDistMatrixDistData_z
( ElConstDistMatrix_z A, ElDistData* distData );

/* mpi::Comm DistMatrix<T,U,V>::DistComm() const
   --------------------------------------------- */
EL_EXPORT ElError ElDistMatrixDistComm_i
( ElConstDistMatrix_i A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixDistComm_s
( ElConstDistMatrix_s A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixDistComm_d
( ElConstDistMatrix_d A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixDistComm_c
( ElConstDistMatrix_c A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixDistComm_z
( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::CrossComm() const
   ---------------------------------------------- */
EL_EXPORT ElError ElDistMatrixCrossComm_i
( ElConstDistMatrix_i A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixCrossComm_s
( ElConstDistMatrix_s A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixCrossComm_d
( ElConstDistMatrix_d A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixCrossComm_c
( ElConstDistMatrix_c A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixCrossComm_z
( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::RedundantComm() const
   -------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixRedundantComm_i
( ElConstDistMatrix_i A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixRedundantComm_s
( ElConstDistMatrix_s A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixRedundantComm_d
( ElConstDistMatrix_d A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixRedundantComm_c
( ElConstDistMatrix_c A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixRedundantComm_z
( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::ColComm() const
   -------------------------------------------- */
EL_EXPORT ElError ElDistMatrixColComm_i
( ElConstDistMatrix_i A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixColComm_s
( ElConstDistMatrix_s A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixColComm_d
( ElConstDistMatrix_d A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixColComm_c
( ElConstDistMatrix_c A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixColComm_z
( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::RowComm() const
   -------------------------------------------- */
EL_EXPORT ElError ElDistMatrixRowComm_i
( ElConstDistMatrix_i A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixRowComm_s
( ElConstDistMatrix_s A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixRowComm_d
( ElConstDistMatrix_d A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixRowComm_c
( ElConstDistMatrix_c A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixRowComm_z
( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::PartialColComm() const
   --------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialColComm_i
( ElConstDistMatrix_i A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialColComm_s
( ElConstDistMatrix_s A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialColComm_d
( ElConstDistMatrix_d A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialColComm_c
( ElConstDistMatrix_c A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialColComm_z
( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::PartialRowComm() const
   --------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialRowComm_i
( ElConstDistMatrix_i A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialRowComm_s
( ElConstDistMatrix_s A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialRowComm_d
( ElConstDistMatrix_d A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialRowComm_c
( ElConstDistMatrix_c A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialRowComm_z
( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::PartialUnionColComm() const
   -------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialUnionColComm_i
( ElConstDistMatrix_i A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialUnionColComm_s
( ElConstDistMatrix_s A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialUnionColComm_d
( ElConstDistMatrix_d A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialUnionColComm_c
( ElConstDistMatrix_c A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialUnionColComm_z
( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::PartialUnionRowComm() const
   -------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialUnionRowComm_i
( ElConstDistMatrix_i A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialUnionRowComm_s
( ElConstDistMatrix_s A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialUnionRowComm_d
( ElConstDistMatrix_d A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialUnionRowComm_c
( ElConstDistMatrix_c A, MPI_Comm* comm );
EL_EXPORT ElError ElDistMatrixPartialUnionRowComm_z
( ElConstDistMatrix_z A, MPI_Comm* comm );

/* Int DistMatrix<T,U,V>::ColStride() const
   ---------------------------------------- */
EL_EXPORT ElError ElDistMatrixColStride_i
( ElConstDistMatrix_i A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixColStride_s
( ElConstDistMatrix_s A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixColStride_d
( ElConstDistMatrix_d A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixColStride_c
( ElConstDistMatrix_c A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixColStride_z
( ElConstDistMatrix_z A, ElInt* stride );

/* Int DistMatrix<T,U,V>::RowStride() const
   ---------------------------------------- */
EL_EXPORT ElError ElDistMatrixRowStride_i
( ElConstDistMatrix_i A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixRowStride_s
( ElConstDistMatrix_s A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixRowStride_d
( ElConstDistMatrix_d A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixRowStride_c
( ElConstDistMatrix_c A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixRowStride_z
( ElConstDistMatrix_z A, ElInt* stride );

/* Int DistMatrix<T,U,V>::PartialColStride() const
   ----------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialColStride_i
( ElConstDistMatrix_i A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialColStride_s
( ElConstDistMatrix_s A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialColStride_d
( ElConstDistMatrix_d A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialColStride_c
( ElConstDistMatrix_c A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialColStride_z
( ElConstDistMatrix_z A, ElInt* stride );

/* Int DistMatrix<T,U,V>::PartialRowStride() const
   ----------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialRowStride_i
( ElConstDistMatrix_i A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialRowStride_s
( ElConstDistMatrix_s A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialRowStride_d
( ElConstDistMatrix_d A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialRowStride_c
( ElConstDistMatrix_c A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialRowStride_z
( ElConstDistMatrix_z A, ElInt* stride );

/* Int DistMatrix<T,U,V>::PartialUnionColStride() const
   ---------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialUnionColStride_i
( ElConstDistMatrix_i A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionColStride_s
( ElConstDistMatrix_s A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionColStride_d
( ElConstDistMatrix_d A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionColStride_c
( ElConstDistMatrix_c A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionColStride_z
( ElConstDistMatrix_z A, ElInt* stride );

/* Int DistMatrix<T,U,V>::PartialUnionRowStride() const
   ---------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialUnionRowStride_i
( ElConstDistMatrix_i A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionRowStride_s
( ElConstDistMatrix_s A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionRowStride_d
( ElConstDistMatrix_d A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionRowStride_c
( ElConstDistMatrix_c A, ElInt* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionRowStride_z
( ElConstDistMatrix_z A, ElInt* stride );

/* Int DistMatrix<T,U,V>::DistSize() const
   --------------------------------------- */
EL_EXPORT ElError ElDistMatrixDistSize_i
( ElConstDistMatrix_i A, ElInt* commSize );
EL_EXPORT ElError ElDistMatrixDistSize_s
( ElConstDistMatrix_s A, ElInt* commSize );
EL_EXPORT ElError ElDistMatrixDistSize_d
( ElConstDistMatrix_d A, ElInt* commSize );
EL_EXPORT ElError ElDistMatrixDistSize_c
( ElConstDistMatrix_c A, ElInt* commSize );
EL_EXPORT ElError ElDistMatrixDistSize_z
( ElConstDistMatrix_z A, ElInt* commSize );

/* Int DistMatrix<T,U,V>::CrossSize() const
   ---------------------------------------- */
EL_EXPORT ElError ElDistMatrixCrossSize_i
( ElConstDistMatrix_i A, ElInt* commSize );
EL_EXPORT ElError ElDistMatrixCrossSize_s
( ElConstDistMatrix_s A, ElInt* commSize );
EL_EXPORT ElError ElDistMatrixCrossSize_d
( ElConstDistMatrix_d A, ElInt* commSize );
EL_EXPORT ElError ElDistMatrixCrossSize_c
( ElConstDistMatrix_c A, ElInt* commSize );
EL_EXPORT ElError ElDistMatrixCrossSize_z
( ElConstDistMatrix_z A, ElInt* commSize );

/* Int DistMatrix<T,U,V>::RedundantSize() const
   -------------------------------------------- */
EL_EXPORT ElError ElDistMatrixRedundantSize_i
( ElConstDistMatrix_i A, ElInt* commSize );
EL_EXPORT ElError ElDistMatrixRedundantSize_s
( ElConstDistMatrix_s A, ElInt* commSize );
EL_EXPORT ElError ElDistMatrixRedundantSize_d
( ElConstDistMatrix_d A, ElInt* commSize );
EL_EXPORT ElError ElDistMatrixRedundantSize_c
( ElConstDistMatrix_c A, ElInt* commSize );
EL_EXPORT ElError ElDistMatrixRedundantSize_z
( ElConstDistMatrix_z A, ElInt* commSize );

/* T DistMatrix<T,U,V>::Get( Int i, Int j ) const
   ---------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGet_i
( ElConstDistMatrix_i A, ElInt i, ElInt j, ElInt* val );
EL_EXPORT ElError ElDistMatrixGet_s
( ElConstDistMatrix_s A, ElInt i, ElInt j, float* val );
EL_EXPORT ElError ElDistMatrixGet_d
( ElConstDistMatrix_d A, ElInt i, ElInt j, double* val );
EL_EXPORT ElError ElDistMatrixGet_c
( ElConstDistMatrix_c A, ElInt i, ElInt j, complex_float* val );
EL_EXPORT ElError ElDistMatrixGet_z
( ElConstDistMatrix_z A, ElInt i, ElInt j, complex_double* val );

/* Base<T> DistMatrix<T,U,V>::GetRealPart( Int i, Int j ) const
   ------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixGetRealPart_c
( ElConstDistMatrix_c A, ElInt i, ElInt j, float* val );
EL_EXPORT ElError ElDistMatrixGetRealPart_z
( ElConstDistMatrix_z A, ElInt i, ElInt j, double* val );

/* Base<T> DistMatrix<T,U,V>::GetImagPart( Int i, Int j ) const
   ------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixGetImagPart_c
( ElConstDistMatrix_c A, ElInt i, ElInt j, float* val );
EL_EXPORT ElError ElDistMatrixGetImagPart_z
( ElConstDistMatrix_z A, ElInt i, ElInt j, double* val );

/* void DistMatrix<T,U,V>::Set( Int i, Int j, T alpha )
   ---------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSet_i
( ElDistMatrix_i A, ElInt i, ElInt j, ElInt alpha );
EL_EXPORT ElError ElDistMatrixSet_s
( ElDistMatrix_s A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElDistMatrixSet_d
( ElDistMatrix_d A, ElInt i, ElInt j, double alpha );
EL_EXPORT ElError ElDistMatrixSet_c
( ElDistMatrix_c A, ElInt i, ElInt j, complex_float alpha );
EL_EXPORT ElError ElDistMatrixSet_z
( ElDistMatrix_z A, ElInt i, ElInt j, complex_double alpha );

/* void DistMatrix<T,U,V>::SetRealPart( Int i, Int j, Base<T> alpha )
   ------------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixSetRealPart_c
( ElDistMatrix_c A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElDistMatrixSetRealPart_z
( ElDistMatrix_z A, ElInt i, ElInt j, double alpha );

/* void DistMatrix<T,U,V>::SetImagPart( Int i, Int j, Base<T> alpha )
   ------------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixSetImagPart_c
( ElDistMatrix_c A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElDistMatrixSetImagPart_z
( ElDistMatrix_z A, ElInt i, ElInt j, double alpha );

/* void DistMatrix<T,U,V>::Update( Int i, Int j, T alpha )
   ------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixUpdate_i
( ElDistMatrix_i A, ElInt i, ElInt j, ElInt alpha );
EL_EXPORT ElError ElDistMatrixUpdate_s
( ElDistMatrix_s A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElDistMatrixUpdate_d
( ElDistMatrix_d A, ElInt i, ElInt j, double alpha );
EL_EXPORT ElError ElDistMatrixUpdate_c
( ElDistMatrix_c A, ElInt i, ElInt j, complex_float alpha );
EL_EXPORT ElError ElDistMatrixUpdate_z
( ElDistMatrix_z A, ElInt i, ElInt j, complex_double alpha );

/* void DistMatrix<T,U,V>::UpdateRealPart( Int i, Int j, Base<T> alpha )
   --------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixUpdateRealPart_c
( ElDistMatrix_c A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElDistMatrixUpdateRealPart_z
( ElDistMatrix_z A, ElInt i, ElInt j, double alpha );

/* void DistMatrix<T,U,V>::UpdateImagPart( Int i, Int j, Base<T> alpha )
   --------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixUpdateImagPart_c
( ElDistMatrix_c A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElDistMatrixUpdateImagPart_z
( ElDistMatrix_z A, ElInt i, ElInt j, double alpha );

/* void DistMatrix<T,U,V>::MakeReal( Int i, Int j )
   ------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixMakeReal_c( ElDistMatrix_c A, ElInt i, ElInt j );
EL_EXPORT ElError ElDistMatrixMakeReal_z( ElDistMatrix_z A, ElInt i, ElInt j );

/* void DistMatrix<T,U,V>::Conjugate( Int i, Int j )
   ------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixConjugate_c( ElDistMatrix_c A, ElInt i, ElInt j );
EL_EXPORT ElError ElDistMatrixConjugate_z( ElDistMatrix_z A, ElInt i, ElInt j );

/* T DistMatrix<T,U,V>::GetLocal( Int iLoc, Int jLoc ) const
   --------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGetLocal_i
( ElConstDistMatrix_i A, ElInt iLoc, ElInt jLoc, ElInt* val );
EL_EXPORT ElError ElDistMatrixGetLocal_s
( ElConstDistMatrix_s A, ElInt iLoc, ElInt jLoc, float* val );
EL_EXPORT ElError ElDistMatrixGetLocal_d
( ElConstDistMatrix_d A, ElInt iLoc, ElInt jLoc, double* val );
EL_EXPORT ElError ElDistMatrixGetLocal_c
( ElConstDistMatrix_c A, ElInt iLoc, ElInt jLoc, complex_float* val );
EL_EXPORT ElError ElDistMatrixGetLocal_z
( ElConstDistMatrix_z A, ElInt iLoc, ElInt jLoc, complex_double* val );

/* Base<T> DistMatrix<T,U,V>::GetLocalRealPart( Int iLoc, Int jLoc ) const
   ----------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGetLocalRealPart_c
( ElConstDistMatrix_c A, ElInt iLoc, ElInt jLoc, float* val );
EL_EXPORT ElError ElDistMatrixGetLocalRealPart_z
( ElConstDistMatrix_z A, ElInt iLoc, ElInt jLoc, double* val );

/* Base<T> DistMatrix<T,U,V>::GetLocalImagPart( Int iLoc, Int jLoc ) const
   ----------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGetLocalImagPart_c
( ElConstDistMatrix_c A, ElInt iLoc, ElInt jLoc, float* val );
EL_EXPORT ElError ElDistMatrixGetLocalImagPart_z
( ElConstDistMatrix_z A, ElInt iLoc, ElInt jLoc, double* val );

/* void DistMatrix<T,U,V>::SetLocal( Int iLoc, Int jLoc, T val )
   ------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSetLocal_i
( ElDistMatrix_i A, ElInt iLoc, ElInt jLoc, ElInt val );
EL_EXPORT ElError ElDistMatrixSetLocal_s
( ElDistMatrix_s A, ElInt iLoc, ElInt jLoc, float val );
EL_EXPORT ElError ElDistMatrixSetLocal_d
( ElDistMatrix_d A, ElInt iLoc, ElInt jLoc, double val );
EL_EXPORT ElError ElDistMatrixSetLocal_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, complex_float val );
EL_EXPORT ElError ElDistMatrixSetLocal_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, complex_double val );

/* void DistMatrix<T,U,V>::SetLocalRealPart( Int iLoc, Int jLoc, Base<T> val )
   --------------------------------------------------------------------------- 
*/
EL_EXPORT ElError ElDistMatrixSetLocalRealPart_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, float val );
EL_EXPORT ElError ElDistMatrixSetLocalRealPart_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, double val );

/* void DistMatrix<T,U,V>::SetLocalImagPart( Int iLoc, Int jLoc, Base<T> val )
   --------------------------------------------------------------------------- 
*/
EL_EXPORT ElError ElDistMatrixSetLocalImagPart_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, float val );
EL_EXPORT ElError ElDistMatrixSetLocalImagPart_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, double val );

/* void DistMatrix<T,U,V>::UpdateLocal( Int iLoc, Int jLoc, T val )
   ---------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixUpdateLocal_i
( ElDistMatrix_i A, ElInt iLoc, ElInt jLoc, ElInt val );
EL_EXPORT ElError ElDistMatrixUpdateLocal_s
( ElDistMatrix_s A, ElInt iLoc, ElInt jLoc, float val );
EL_EXPORT ElError ElDistMatrixUpdateLocal_d
( ElDistMatrix_d A, ElInt iLoc, ElInt jLoc, double val );
EL_EXPORT ElError ElDistMatrixUpdateLocal_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, complex_float val );
EL_EXPORT ElError ElDistMatrixUpdateLocal_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, complex_double val );

/* void DistMatrix<T,U,V>::UpdateLocalRealPart
   ( Int iLoc, Int jLoc, Base<T> val )
   ------------------------------------------- 
*/
EL_EXPORT ElError ElDistMatrixUpdateLocalRealPart_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, float val );
EL_EXPORT ElError ElDistMatrixUpdateLocalRealPart_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, double val );

/* void DistMatrix<T,U,V>::UpdateLocalImagPart
   ( Int iLoc, Int jLoc, Base<T> val )
   ------------------------------------------- 
*/
EL_EXPORT ElError ElDistMatrixUpdateLocalImagPart_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, float val );
EL_EXPORT ElError ElDistMatrixUpdateLocalImagPart_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, double val );

/* void DistMatrix<T,U,V>::MakeDiagonalReal( Int offset )
   ------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixMakeDiagonalReal_c
( ElDistMatrix_c A, ElInt offset );
EL_EXPORT ElError ElDistMatrixMakeDiagonalReal_z
( ElDistMatrix_z A, ElInt offset );

/* void DistMatrix<T,U,V>::ConjugateDiagonal( Int offset )
   ------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixConjugateDiagonal_c
( ElDistMatrix_c A, ElInt offset );
EL_EXPORT ElError ElDistMatrixConjugateDiagonal_z
( ElDistMatrix_z A, ElInt offset );

/* bool DistMatrix<T,U,V>::DiagonalAlignedWith
   ( const DistData& data, Int offset ) const
   ------------------------------------------- */
EL_EXPORT ElError ElDistMatrixAlignedWith_i
( ElConstDistMatrix_i A, ElDistData distData, ElInt offset, bool* aligned );
EL_EXPORT ElError ElDistMatrixAlignedWith_s
( ElConstDistMatrix_s A, ElDistData distData, ElInt offset, bool* aligned );
EL_EXPORT ElError ElDistMatrixAlignedWith_d
( ElConstDistMatrix_d A, ElDistData distData, ElInt offset, bool* aligned );
EL_EXPORT ElError ElDistMatrixAlignedWith_c
( ElConstDistMatrix_c A, ElDistData distData, ElInt offset, bool* aligned );
EL_EXPORT ElError ElDistMatrixAlignedWith_z
( ElConstDistMatrix_z A, ElDistData distData, ElInt offset, bool* aligned );

/* Int DistMatrix<T,U,V>::DiagonalRoot( Int offset ) const
   ------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixDiagonalRoot_i
( ElConstDistMatrix_i A, ElInt offset, ElInt* root );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_s
( ElConstDistMatrix_s A, ElInt offset, ElInt* root );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_d
( ElConstDistMatrix_d A, ElInt offset, ElInt* root );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_c
( ElConstDistMatrix_c A, ElInt offset, ElInt* root );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_z
( ElConstDistMatrix_z A, ElInt offset, ElInt* root );

/* Int DistMatrix<T,U,V>::DiagonalAlign( Int offset ) const
   -------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixDiagonalRoot_i
( ElConstDistMatrix_i A, ElInt offset, ElInt* align );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_s
( ElConstDistMatrix_s A, ElInt offset, ElInt* align );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_d
( ElConstDistMatrix_d A, ElInt offset, ElInt* align );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_c
( ElConstDistMatrix_c A, ElInt offset, ElInt* align );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_z
( ElConstDistMatrix_z A, ElInt offset, ElInt* align );

/* DistMatrix<T,UDiag,VDiag> 
   DistMatrix<T,U,V>::GetDiagonal( Int offset ) const
   -------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGetDiagonal_i
( ElConstDistMatrix_i A, ElInt offset, ElDistMatrix_i* d );
EL_EXPORT ElError ElDistMatrixGetDiagonal_s
( ElConstDistMatrix_s A, ElInt offset, ElDistMatrix_s* d );
EL_EXPORT ElError ElDistMatrixGetDiagonal_d
( ElConstDistMatrix_d A, ElInt offset, ElDistMatrix_d* d );
EL_EXPORT ElError ElDistMatrixGetDiagonal_c
( ElConstDistMatrix_c A, ElInt offset, ElDistMatrix_c* d );
EL_EXPORT ElError ElDistMatrixGetDiagonal_z
( ElConstDistMatrix_z A, ElInt offset, ElDistMatrix_z* d );

/* DistMatrix<Base<T>,UDiag,VDiag> 
   DistMatrix<T,U,V>::GetRealPartOfDiagonal( Int offset ) const
   ------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixGetRealPartOfDiagonal_c
( ElConstDistMatrix_c A, ElInt offset, ElDistMatrix_s* d );
EL_EXPORT ElError ElDistMatrixGetRealPartOfDiagonal_z
( ElConstDistMatrix_z A, ElInt offset, ElDistMatrix_d* d );

/* DistMatrix<Base<T>,UDiag,VDiag> 
   DistMatrix<T,U,V>::GetImagPartOfDiagonal( Int offset ) const
   ------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixGetImagPartOfDiagonal_c
( ElConstDistMatrix_c A, ElInt offset, ElDistMatrix_s* d );
EL_EXPORT ElError ElDistMatrixGetImagPartOfDiagonal_z
( ElConstDistMatrix_z A, ElInt offset, ElDistMatrix_d* d );

/* TODO: More diagonal manipulation
   ================================ */

/* DistMatrix<T,STAR,STAR> DistMatrix<T,U,V>::GetSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J ) const
   -------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGetSubmatrix_i
( ElConstDistMatrix_i A,
  ElInt numRowInds, const ElInt* rowInd,
  ElInt numColInds, const ElInt* colInd, ElDistMatrix_i* ASub );
EL_EXPORT ElError ElDistMatrixGetSubmatrix_s
( ElConstDistMatrix_s A,
  ElInt numRowInds, const ElInt* rowInd,
  ElInt numColInds, const ElInt* colInd, ElDistMatrix_s* ASub );
EL_EXPORT ElError ElDistMatrixGetSubmatrix_d
( ElConstDistMatrix_d A,
  ElInt numRowInds, const ElInt* rowInd,
  ElInt numColInds, const ElInt* colInd, ElDistMatrix_d* ASub );
EL_EXPORT ElError ElDistMatrixGetSubmatrix_c
( ElConstDistMatrix_c A,
  ElInt numRowInds, const ElInt* rowInd,
  ElInt numColInds, const ElInt* colInd, ElDistMatrix_c* ASub );
EL_EXPORT ElError ElDistMatrixGetSubmatrix_z
( ElConstDistMatrix_z A,
  ElInt numRowInds, const ElInt* rowInd,
  ElInt numColInds, const ElInt* colInd, ElDistMatrix_z* ASub );

/* DistMatrix<Base<T>,STAR,STAR> DistMatrix<T,U,V>::GetRealPartOfSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J ) const
   ----------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGetRealPartOfSubmatrix_c
( ElConstDistMatrix_c A,
  ElInt numRowInds, const ElInt* rowInd,
  ElInt numColInds, const ElInt* colInd, ElDistMatrix_s* ASub );
EL_EXPORT ElError ElDistMatrixGetRealPartOfSubmatrix_z
( ElConstDistMatrix_z A,
  ElInt numRowInds, const ElInt* rowInd,
  ElInt numColInds, const ElInt* colInd, ElDistMatrix_d* ASub );

/* DistMatrix<Base<T>,STAR,STAR> DistMatrix<T,U,V>::GetImagPartOfSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J ) const
   ----------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGetImagPartOfSubmatrix_c
( ElConstDistMatrix_c A,
  ElInt numRowInds, const ElInt* rowInd,
  ElInt numColInds, const ElInt* colInd, ElDistMatrix_s* ASub );
EL_EXPORT ElError ElDistMatrixGetImagPartOfSubmatrix_z
( ElConstDistMatrix_z A,
  ElInt numRowInds, const ElInt* rowInd,
  ElInt numColInds, const ElInt* colInd, ElDistMatrix_d* ASub );

/* void DistMatrix<T,U,V>::SetSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J, 
     const DistMatrix<T,STAR,STAR>& ASub )
   ------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSetSubmatrix_i
( ElDistMatrix_i A, const ElInt* rowInd, const ElInt* colInd, 
  ElConstDistMatrix_i ASub );
EL_EXPORT ElError ElDistMatrixSetSubmatrix_s
( ElDistMatrix_s A, const ElInt* rowInd, const ElInt* colInd, 
  ElConstDistMatrix_s ASub );
EL_EXPORT ElError ElDistMatrixSetSubmatrix_d
( ElDistMatrix_d A, const ElInt* rowInd, const ElInt* colInd, 
  ElConstDistMatrix_d ASub );
EL_EXPORT ElError ElDistMatrixSetSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowInd, const ElInt* colInd, 
  ElConstDistMatrix_c ASub );
EL_EXPORT ElError ElDistMatrixSetSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowInd, const ElInt* colInd, 
  ElConstDistMatrix_z ASub );

/* void DistMatrix<T,U,V>::SetRealPartOfSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J, 
     const DistMatrix<Base<T>,STAR,STAR>& ASub )
   ------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSetRealPartOfSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowInd, const ElInt* colInd, 
  ElConstDistMatrix_s ASub );
EL_EXPORT ElError ElDistMatrixSetRealPartOfSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowInd, const ElInt* colInd, 
  ElConstDistMatrix_d ASub );

/* void DistMatrix<T,U,V>::SetImagPartOfSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J, 
     const DistMatrix<Base<T>,STAR,STAR>& ASub )
   ------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSetImagPartOfSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowInd, const ElInt* colInd, 
  ElConstDistMatrix_s ASub );
EL_EXPORT ElError ElDistMatrixSetImagPartOfSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowInd, const ElInt* colInd, 
  ElConstDistMatrix_d ASub );

/* void DistMatrix<T,U,V>::UpdateSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J, 
     T alpha, const DistMatrix<T,STAR,STAR>& ASub )
   ------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixUpdateSubmatrix_i
( ElDistMatrix_i A, const ElInt* rowInd, const ElInt* colInd, 
  ElInt alpha, ElConstDistMatrix_i ASub );
EL_EXPORT ElError ElDistMatrixUpdateSubmatrix_s
( ElDistMatrix_s A, const ElInt* rowInd, const ElInt* colInd, 
  float alpha, ElConstDistMatrix_s ASub );
EL_EXPORT ElError ElDistMatrixUpdateSubmatrix_d
( ElDistMatrix_d A, const ElInt* rowInd, const ElInt* colInd, 
  double alpha, ElConstDistMatrix_d ASub );
EL_EXPORT ElError ElDistMatrixUpdateSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowInd, const ElInt* colInd, 
  complex_float alpha, ElConstDistMatrix_c ASub );
EL_EXPORT ElError ElDistMatrixUpdateSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowInd, const ElInt* colInd, 
  complex_double alpha, ElConstDistMatrix_z ASub );

/* void DistMatrix<T,U,V>::UpdateRealPartOfSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J, 
     Base<T> alpha, const DistMatrix<Base<T>,STAR,STAR>& ASub )
   ------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixUpdateRealPartOfSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowInd, const ElInt* colInd, 
  float alpha, ElConstDistMatrix_s ASub );
EL_EXPORT ElError ElDistMatrixUpdateRealPartOfSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowInd, const ElInt* colInd, 
  double alpha, ElConstDistMatrix_d ASub );

/* void DistMatrix<T,U,V>::UpdateImagPartOfSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J, 
     Base<T> alpha, const DistMatrix<Base<T>,STAR,STAR>& ASub )
   ------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixUpdateImagPartOfSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowInd, const ElInt* colInd, 
  float alpha, ElConstDistMatrix_s ASub );
EL_EXPORT ElError ElDistMatrixUpdateImagPartOfSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowInd, const ElInt* colInd, 
  double alpha, ElConstDistMatrix_d ASub );

/* void DistMatrix<T,U,V>::MakeSubmatrixReal
   ( const std::vector<Int>& I, const std::vector<Int>& J )
   -------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixMakeSubmatrixReal_c
( ElDistMatrix_c A, ElInt numRowInds, const ElInt* rowInd, 
                    ElInt numColInds, const ElInt* colInd );
EL_EXPORT ElError ElDistMatrixMakeSubmatrixReal_z
( ElDistMatrix_z A, ElInt numRowInds, const ElInt* rowInd, 
                    ElInt numColInds, const ElInt* colInd );

/* void DistMatrix<T,U,V>::ConjugateSubmatrix
   ( const std::vector<Int>& I, const std::vector<Int>& J )
   -------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixConjugateSubmatrix_c
( ElDistMatrix_c A, ElInt numRowInds, const ElInt* rowInd, 
                    ElInt numColInds, const ElInt* colInd );
EL_EXPORT ElError ElDistMatrixConjugateSubmatrix_z
( ElDistMatrix_z A, ElInt numRowInds, const ElInt* rowInd, 
                    ElInt numColInds, const ElInt* colInd );

/* Matrix<T> DistMatrix<T,U,V>::GetLocalSubmatrix
   ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc ) const
   -------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGetLocalSubmatrix_i
( ElConstDistMatrix_i A,
  ElInt numRowInds, const ElInt* rowIndLoc,
  ElInt numColInds, const ElInt* colIndLoc, ElMatrix_i* ASub );
EL_EXPORT ElError ElDistMatrixGetLocalSubmatrix_s
( ElConstDistMatrix_s A,
  ElInt numRowInds, const ElInt* rowIndLoc,
  ElInt numColInds, const ElInt* colIndLoc, ElMatrix_s* ASub );
EL_EXPORT ElError ElDistMatrixGetLocalSubmatrix_d
( ElConstDistMatrix_d A,
  ElInt numRowInds, const ElInt* rowIndLoc,
  ElInt numColInds, const ElInt* colIndLoc, ElMatrix_d* ASub );
EL_EXPORT ElError ElDistMatrixGetLocalSubmatrix_c
( ElConstDistMatrix_c A,
  ElInt numRowInds, const ElInt* rowIndLoc,
  ElInt numColInds, const ElInt* colIndLoc, ElMatrix_c* ASub );
EL_EXPORT ElError ElDistMatrixGetLocalSubmatrix_z
( ElConstDistMatrix_z A,
  ElInt numRowInds, const ElInt* rowIndLoc,
  ElInt numColInds, const ElInt* colIndLoc, ElMatrix_z* ASub );

/* Matrix<Base<T>> DistMatrix<T,U,V>::GetRealPartOfLocalSubmatrix
   ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc ) const
   -------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGetRealPartOfLocalSubmatrix_c
( ElConstDistMatrix_c A,
  ElInt numRowInds, const ElInt* rowIndLoc,
  ElInt numColInds, const ElInt* colIndLoc, ElMatrix_s* ASub );
EL_EXPORT ElError ElDistMatrixGetRealPartOfLocalSubmatrix_z
( ElConstDistMatrix_z A,
  ElInt numRowInds, const ElInt* rowIndLoc,
  ElInt numColInds, const ElInt* colIndLoc, ElMatrix_d* ASub );

/* Matrix<Base<T>> DistMatrix<T,U,V>::GetImagPartOfLocalSubmatrix
   ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc ) const
   -------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGetImagPartOfLocalSubmatrix_c
( ElConstDistMatrix_c A,
  ElInt numRowInds, const ElInt* rowIndLoc,
  ElInt numColInds, const ElInt* colIndLoc, ElMatrix_s* ASub );
EL_EXPORT ElError ElDistMatrixGetImagPartOfLocalSubmatrix_z
( ElConstDistMatrix_z A,
  ElInt numRowInds, const ElInt* rowIndLoc,
  ElInt numColInds, const ElInt* colIndLoc, ElMatrix_d* ASub );

/* void DistMatrix<T,U,V>::SetLocalSubmatrix
   ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc, 
     const Matrix<T>& ASub )
   ---------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSetLocalSubmatrix_i
( ElDistMatrix_i A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  ElConstMatrix_i ASub );
EL_EXPORT ElError ElDistMatrixSetLocalSubmatrix_s
( ElDistMatrix_s A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  ElConstMatrix_s ASub );
EL_EXPORT ElError ElDistMatrixSetLocalSubmatrix_d
( ElDistMatrix_d A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  ElConstMatrix_d ASub );
EL_EXPORT ElError ElDistMatrixSetLocalSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  ElConstMatrix_c ASub );
EL_EXPORT ElError ElDistMatrixSetLocalSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  ElConstMatrix_z ASub );

/* void DistMatrix<T,U,V>::SetRealPartOfLocalSubmatrix
   ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc, 
     const Matrix<Base<T>>& ASub )
   ---------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSetRealPartOfLocalSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  ElConstMatrix_s ASub );
EL_EXPORT ElError ElDistMatrixSetRealPartOfLocalSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  ElConstMatrix_d ASub );

/* void DistMatrix<T,U,V>::SetImagPartOfLocalSubmatrix
   ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc, 
     const Matrix<Base<T>>& ASub )
   ---------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSetImagPartOfLocalSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  ElConstMatrix_s ASub );
EL_EXPORT ElError ElDistMatrixSetImagPartOfLocalSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  ElConstMatrix_d ASub );

/* void DistMatrix<T,U,V>::UpdateLocalSubmatrix
   ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc, 
     T alpha, const Matrix<T>& ASub )
   ----------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixUpdateLocalSubmatrix_i
( ElDistMatrix_i A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  ElInt alpha, ElConstMatrix_i ASub );
EL_EXPORT ElError ElDistMatrixUpdateLocalSubmatrix_s
( ElDistMatrix_s A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  float alpha, ElConstMatrix_s ASub );
EL_EXPORT ElError ElDistMatrixUpdateLocalSubmatrix_d
( ElDistMatrix_d A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  double alpha, ElConstMatrix_d ASub );
EL_EXPORT ElError ElDistMatrixUpdateLocalSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  complex_float alpha, ElConstMatrix_c ASub );
EL_EXPORT ElError ElDistMatrixUpdateLocalSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  complex_double alpha, ElConstMatrix_z ASub );

/* void DistMatrix<T,U,V>::UpdateRealPartOfLocalSubmatrix
   ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc, 
     Base<T> alpha, const Matrix<Base<T>>& ASub )
   ----------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixUpdateRealPartOfLocalSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  float alpha, ElConstMatrix_s ASub );
EL_EXPORT ElError ElDistMatrixUpdateRealPartOfLocalSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  double alpha, ElConstMatrix_d ASub );

/* void DistMatrix<T,U,V>::UpdateImagPartOfLocalSubmatrix
   ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc, 
     Base<T> alpha, const Matrix<Base<T>>& ASub )
   ----------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixUpdateImagPartOfLocalSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  float alpha, ElConstMatrix_s ASub );
EL_EXPORT ElError ElDistMatrixUpdateImagPartOfLocalSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowIndLoc, const ElInt* colIndLoc, 
  double alpha, ElConstMatrix_d ASub );

/* void DistMatrix<T,U,V>::MakeLocalSubmatrixReal
   ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc )
   ------------------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixMakeLocalSubmatrixReal_c
( ElDistMatrix_c A, ElInt numRowInds, const ElInt* rowIndLoc, 
                    ElInt numColInds, const ElInt* colIndLoc );
EL_EXPORT ElError ElDistMatrixMakeLocalSubmatrixReal_z
( ElDistMatrix_z A, ElInt numRowInds, const ElInt* rowIndLoc, 
                    ElInt numColInds, const ElInt* colIndLoc );

/* void DistMatrix<T,U,V>::ConjugateLocalSubmatrix
   ( const std::vector<Int>& ILoc, const std::vector<Int>& JLoc )
   ------------------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixConjugateLocalSubmatrix_c
( ElDistMatrix_c A, ElInt numRowInds, const ElInt* rowIndLoc, 
                    ElInt numColInds, const ElInt* colIndLoc );
EL_EXPORT ElError ElDistMatrixConjugateLocalSubmatrix_z
( ElDistMatrix_z A, ElInt numRowInds, const ElInt* rowIndLoc, 
                    ElInt numColInds, const ElInt* colIndLoc );

/* void DistMatrix<T,U,V>::SumOver( mpi::Comm comm )
   ------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSumOver_i( ElDistMatrix_i A, MPI_Comm comm );
EL_EXPORT ElError ElDistMatrixSumOver_s( ElDistMatrix_s A, MPI_Comm comm );
EL_EXPORT ElError ElDistMatrixSumOver_d( ElDistMatrix_d A, MPI_Comm comm );
EL_EXPORT ElError ElDistMatrixSumOver_c( ElDistMatrix_c A, MPI_Comm comm );
EL_EXPORT ElError ElDistMatrixSumOver_z( ElDistMatrix_z A, MPI_Comm comm );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_DISTMATRIX_C_H */
