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

typedef struct 
{
    ElDist colDist, rowDist;
    ElInt colAlign, rowAlign;
    ElInt root;
    ElConstGrid grid;
} ElDistData;

/* An anonymous struct meant as a placeholder for AbstractDistMatrix<T>
   -------------------------------------------------------------------- */
typedef struct ElDistMatrix_sDummy* ElDistMatrix_s;
typedef struct ElDistMatrix_dDummy* ElDistMatrix_d;
typedef struct ElDistMatrix_cDummy* ElDistMatrix_c;
typedef struct ElDistMatrix_zDummy* ElDistMatrix_z;

typedef const struct ElDistMatrix_sDummy* ElConstDistMatrix_s;
typedef const struct ElDistMatrix_dDummy* ElConstDistMatrix_d;
typedef const struct ElDistMatrix_cDummy* ElConstDistMatrix_c;
typedef const struct ElDistMatrix_zDummy* ElConstDistMatrix_z;

/* DistMatrix<T,MC,MR>::DistMatrix( const Grid& g )
   ------------------------------------------------ */
ElError ElDistMatrixCreate_s( ElConstGrid g, ElDistMatrix_s* A );
ElError ElDistMatrixCreate_d( ElConstGrid g, ElDistMatrix_d* A );
ElError ElDistMatrixCreate_c( ElConstGrid g, ElDistMatrix_c* A );
ElError ElDistMatrixCreate_z( ElConstGrid g, ElDistMatrix_z* A );

/* DistMatrix<T,U,V>::DistMatrix( const Grid& g )
   ---------------------------------------------- */
ElError ElDistMatrixCreateSpecific_s
( ElDist U, ElDist V, ElConstGrid g, ElDistMatrix_s* A );
ElError ElDistMatrixCreateSpecific_d
( ElDist U, ElDist V, ElConstGrid g, ElDistMatrix_d* A );
ElError ElDistMatrixCreateSpecific_c
( ElDist U, ElDist V, ElConstGrid g, ElDistMatrix_c* A );
ElError ElDistMatrixCreateSpecific_z
( ElDist U, ElDist V, ElConstGrid g, ElDistMatrix_z* A );

/* DistMatrix<T,U,V>::~DistMatrix()
   -------------------------------- */
ElError ElDistMatrixDestroy_s( ElConstDistMatrix_s A );
ElError ElDistMatrixDestroy_d( ElConstDistMatrix_d A );
ElError ElDistMatrixDestroy_c( ElConstDistMatrix_c A );
ElError ElDistMatrixDestroy_z( ElConstDistMatrix_z A );

/* void DistMatrix<T,U,V>::Empty()
   ------------------------------- */
ElError ElDistMatrixEmpty_s( ElDistMatrix_s A );
ElError ElDistMatrixEmpty_d( ElDistMatrix_d A );
ElError ElDistMatrixEmpty_c( ElDistMatrix_c A );
ElError ElDistMatrixEmpty_z( ElDistMatrix_z A );

/* void DistMatrix<T,U,V>::EmptyData()
   ----------------------------------- */
ElError ElDistMatrixEmptyData_s( ElDistMatrix_s A );
ElError ElDistMatrixEmptyData_d( ElDistMatrix_d A );
ElError ElDistMatrixEmptyData_c( ElDistMatrix_c A );
ElError ElDistMatrixEmptyData_z( ElDistMatrix_z A );

/* void DistMatrix<T,U,V>::SetGrid( const Grid& g )
   ------------------------------------------------ */
ElError ElDistMatrixSetGrid_s( ElDistMatrix_s A, ElConstGrid grid );
ElError ElDistMatrixSetGrid_d( ElDistMatrix_d A, ElConstGrid grid );
ElError ElDistMatrixSetGrid_c( ElDistMatrix_c A, ElConstGrid grid );
ElError ElDistMatrixSetGrid_z( ElDistMatrix_z A, ElConstGrid grid );

/* B = A
   ----- */
ElError ElDistMatrixCopy_s( ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElDistMatrixCopy_d( ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElDistMatrixCopy_c( ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElDistMatrixCopy_z( ElConstDistMatrix_z A, ElDistMatrix_z B );

/* void DistMatrix<T,U,V>::Resize( Int height, Int width )
   ------------------------------------------------------- */
ElError ElDistMatrixResize_s( ElDistMatrix_s A, ElInt height, ElInt width );
ElError ElDistMatrixResize_d( ElDistMatrix_d A, ElInt height, ElInt width );
ElError ElDistMatrixResize_c( ElDistMatrix_c A, ElInt height, ElInt width );
ElError ElDistMatrixResize_z( ElDistMatrix_z A, ElInt height, ElInt width );

/* void DistMatrix<T,U,V>::Resize( Int height, Int width, Int ldim )
   ----------------------------------------------------------------- */
ElError ElDistMatrixResizeWithLDim_s
( ElDistMatrix_s A, ElInt height, ElInt width, ElInt ldim );
ElError ElDistMatrixResizeWithLDim_d
( ElDistMatrix_d A, ElInt height, ElInt width, ElInt ldim );
ElError ElDistMatrixResizeWithLDim_c
( ElDistMatrix_c A, ElInt height, ElInt width, ElInt ldim );
ElError ElDistMatrixResizeWithLDim_z
( ElDistMatrix_z A, ElInt height, ElInt width, ElInt ldim );

/* void DistMatrix<T,U,V>::MakeConsistent( bool includeViewers )
   ------------------------------------------------------------- */
ElError ElDistMatrixMakeConsistent_s( ElDistMatrix_s A, bool includeViewers );
ElError ElDistMatrixMakeConsistent_d( ElDistMatrix_d A, bool includeViewers );
ElError ElDistMatrixMakeConsistent_c( ElDistMatrix_c A, bool includeViewers );
ElError ElDistMatrixMakeConsistent_z( ElDistMatrix_z A, bool includeViewers );

/* void DistMatrix<T,U,V>::MakeSizeConsistent( bool includeViewers )
   ----------------------------------------------------------------- */
ElError ElDistMatrixMakeSizeConsistent_s
( ElDistMatrix_s A, bool includeViewers );
ElError ElDistMatrixMakeSizeConsistent_d
( ElDistMatrix_d A, bool includeViewers );
ElError ElDistMatrixMakeSizeConsistent_c
( ElDistMatrix_c A, bool includeViewers );
ElError ElDistMatrixMakeSizeConsistent_z
( ElDistMatrix_z A, bool includeViewers );

/* void DistMatrix<T,U,V>::Align( Int colAlign, Int rowAlign, bool constrain )
   --------------------------------------------------------------------------- 
*/
ElError ElDistMatrixAlign_s
( ElDistMatrix_s A, ElInt colAlign, ElInt rowAlign, bool constrain );
ElError ElDistMatrixAlign_d
( ElDistMatrix_d A, ElInt colAlign, ElInt rowAlign, bool constrain );
ElError ElDistMatrixAlign_c
( ElDistMatrix_c A, ElInt colAlign, ElInt rowAlign, bool constrain );
ElError ElDistMatrixAlign_z
( ElDistMatrix_z A, ElInt colAlign, ElInt rowAlign, bool constrain );

/* void DistMatrix<T,U,V>::AlignCols( Int colAlign, bool constrain )
   ----------------------------------------------------------------- */
ElError ElDistMatrixAlignCols_s
( ElDistMatrix_s A, ElInt colAlign, bool constrain );
ElError ElDistMatrixAlignCols_d
( ElDistMatrix_d A, ElInt colAlign, bool constrain );
ElError ElDistMatrixAlignCols_c
( ElDistMatrix_c A, ElInt colAlign, bool constrain );
ElError ElDistMatrixAlignCols_z
( ElDistMatrix_z A, ElInt colAlign, bool constrain );

/* void DistMatrix<T,U,V>::AlignRows( Int rowAlign, bool constrain )
   ----------------------------------------------------------------- */
ElError ElDistMatrixAlignRows_s
( ElDistMatrix_s A, ElInt rowAlign, bool constrain );
ElError ElDistMatrixAlignRows_d
( ElDistMatrix_d A, ElInt rowAlign, bool constrain );
ElError ElDistMatrixAlignRows_c
( ElDistMatrix_c A, ElInt rowAlign, bool constrain );
ElError ElDistMatrixAlignRows_z
( ElDistMatrix_z A, ElInt rowAlign, bool constrain );

/* void DistMatrix<T,U,V>::FreeAlignments()
   ---------------------------------------- */
ElError ElDistMatrixFreeAlignments_s( ElDistMatrix_s A );
ElError ElDistMatrixFreeAlignments_d( ElDistMatrix_d A );
ElError ElDistMatrixFreeAlignments_c( ElDistMatrix_c A );
ElError ElDistMatrixFreeAlignments_z( ElDistMatrix_z A );

/* void DistMatrix<T,U,V>::SetRoot( Int root )
   ------------------------------------------- */
ElError ElDistMatrixSetRoot_s( ElDistMatrix_s A, ElInt root );
ElError ElDistMatrixSetRoot_d( ElDistMatrix_d A, ElInt root );
ElError ElDistMatrixSetRoot_c( ElDistMatrix_c A, ElInt root );
ElError ElDistMatrixSetRoot_z( ElDistMatrix_z A, ElInt root );

/* void DistMatrix<T,U,V>::AlignWith( const DistData& data, bool constrain )
   ------------------------------------------------------------------------- */
ElError ElDistMatrixAlignWith_s
( ElDistMatrix_s A, ElDistData distData, bool constrain );
ElError ElDistMatrixAlignWith_d
( ElDistMatrix_d A, ElDistData distData, bool constrain );
ElError ElDistMatrixAlignWith_c
( ElDistMatrix_c A, ElDistData distData, bool constrain );
ElError ElDistMatrixAlignWith_z
( ElDistMatrix_z A, ElDistData distData, bool constrain );

/* void DistMatrix<T,U,V>::AlignColsWith( const DistData& data, bool constrain )
   -----------------------------------------------------------------------------
*/
ElError ElDistMatrixAlignColsWith_s
( ElDistMatrix_s A, ElDistData distData, bool constrain );
ElError ElDistMatrixAlignColsWith_d
( ElDistMatrix_d A, ElDistData distData, bool constrain );
ElError ElDistMatrixAlignColsWith_c
( ElDistMatrix_c A, ElDistData distData, bool constrain );
ElError ElDistMatrixAlignColsWith_z
( ElDistMatrix_z A, ElDistData distData, bool constrain );

/* void DistMatrix<T,U,V>::AlignRowsWith( const DistData& data, bool constrain )
   -----------------------------------------------------------------------------
*/
ElError ElDistMatrixAlignRowsWith_s
( ElDistMatrix_s A, ElDistData distData, bool constrain );
ElError ElDistMatrixAlignRowsWith_d
( ElDistMatrix_d A, ElDistData distData, bool constrain );
ElError ElDistMatrixAlignRowsWith_c
( ElDistMatrix_c A, ElDistData distData, bool constrain );
ElError ElDistMatrixAlignRowsWith_z
( ElDistMatrix_z A, ElDistData distData, bool constrain );

/* void DistMatrix<T,U,V>::AlignAndResize
  ( Int colAlign, Int rowAlign, Int height, Int width, 
    bool force, bool constrain )
   --------------------------------------------------- */
ElError ElDistMatrixAlignAndResize_s
( ElDistMatrix_s A, ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
ElError ElDistMatrixAlignAndResize_d
( ElDistMatrix_d A, ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
ElError ElDistMatrixAlignAndResize_c
( ElDistMatrix_c A, ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
ElError ElDistMatrixAlignAndResize_z
( ElDistMatrix_z A, ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );

/* void DistMatrix<T,U,V>::AlignColsAndResize
   ( Int colAlign, Int height, Int width, bool force, bool constrain )
   ------------------------------------------------------------------- */
ElError ElDistMatrixAlignColsAndResize_s
( ElDistMatrix_s A, ElInt colAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
ElError ElDistMatrixAlignColsAndResize_d
( ElDistMatrix_d A, ElInt colAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
ElError ElDistMatrixAlignColsAndResize_c
( ElDistMatrix_c A, ElInt colAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
ElError ElDistMatrixAlignColsAndResize_z
( ElDistMatrix_z A, ElInt colAlign, ElInt height, ElInt width, 
  bool force, bool constrain );

/* void DistMatrix<T,U,V>::AlignRowsAndResize
  ( Int rowAlign, Int height, Int width, bool force, bool constrain )
   ------------------------------------------------------------------ */
ElError ElDistMatrixAlignRowsAndResize_s
( ElDistMatrix_s A, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
ElError ElDistMatrixAlignRowsAndResize_d
( ElDistMatrix_d A, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
ElError ElDistMatrixAlignRowsAndResize_c
( ElDistMatrix_c A, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
ElError ElDistMatrixAlignRowsAndResize_z
( ElDistMatrix_z A, ElInt rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );

/* void DistMatrix<T,U,V>::Attach
   ( Int height, Int width, const Grid& grid, Int colAlign, Int rowAlign, 
     T* buffer, Int ldim, Int root )
   ---------------------------------------------------------------------- */
ElError ElDistMatrixAttach_s
( ElDistMatrix_s A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, float* buffer, ElInt ldim, 
  ElInt root );
ElError ElDistMatrixAttach_d
( ElDistMatrix_d A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, double* buffer, ElInt ldim, 
  ElInt root );
ElError ElDistMatrixAttach_c
( ElDistMatrix_c A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, complex_float* buffer, ElInt ldim, 
  ElInt root );
ElError ElDistMatrixAttach_z
( ElDistMatrix_z A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, complex_double* buffer, ElInt ldim, 
  ElInt root );

/* void DistMatrix<T,U,V>::LockedAttach
   ( Int height, Int width, const Grid& grid, Int colAlign, Int rowAlign, 
     const T* buffer, Int ldim, Int root )
   ---------------------------------------------------------------------- */
ElError ElDistMatrixLockedAttach_s
( ElDistMatrix_s A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, const float* buffer, 
  ElInt ldim, ElInt root );
ElError ElDistMatrixLockedAttach_d
( ElDistMatrix_d A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, const double* buffer, 
  ElInt ldim, ElInt root );
ElError ElDistMatrixLockedAttach_c
( ElDistMatrix_c A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, const complex_float* buffer, 
  ElInt ldim, ElInt root );
ElError ElDistMatrixLockedAttach_z
( ElDistMatrix_z A, ElInt height, ElInt width, ElConstGrid grid,
  ElInt colAlign, ElInt rowAlign, const complex_double* buffer, 
  ElInt ldim, ElInt root );

/* Int DistMatrix<T,U,V>::Height() const
   ------------------------------------- */
ElError ElDistMatrixHeight_s( ElConstDistMatrix_s A, ElInt* height );
ElError ElDistMatrixHeight_d( ElConstDistMatrix_d A, ElInt* height );
ElError ElDistMatrixHeight_c( ElConstDistMatrix_c A, ElInt* height );
ElError ElDistMatrixHeight_z( ElConstDistMatrix_z A, ElInt* height );

/* Int DistMatrix<T,U,V>::Width() const
   ------------------------------------ */
ElError ElDistMatrixWidth_s( ElConstDistMatrix_s A, ElInt* width );
ElError ElDistMatrixWidth_d( ElConstDistMatrix_d A, ElInt* width );
ElError ElDistMatrixWidth_c( ElConstDistMatrix_c A, ElInt* width );
ElError ElDistMatrixWidth_z( ElConstDistMatrix_z A, ElInt* width );

/* Int DistMatrix<T,U,V>::DiagonalLength( Int offset ) const
   --------------------------------------------------------- */
ElError ElDistMatrixDiagonalLength_s
( ElConstDistMatrix_s A, ElInt offset, ElInt* length );
ElError ElDistMatrixDiagonalLength_d
( ElConstDistMatrix_d A, ElInt offset, ElInt* length );
ElError ElDistMatrixDiagonalLength_c
( ElConstDistMatrix_c A, ElInt offset, ElInt* length );
ElError ElDistMatrixDiagonalLength_z
( ElConstDistMatrix_z A, ElInt offset, ElInt* length );

/* bool DistMatrix<T,U,V>::Viewing() const
   --------------------------------------- */
ElError ElDistMatrixViewing_s( ElConstDistMatrix_s A, bool* viewing );
ElError ElDistMatrixViewing_d( ElConstDistMatrix_d A, bool* viewing );
ElError ElDistMatrixViewing_c( ElConstDistMatrix_c A, bool* viewing );
ElError ElDistMatrixViewing_z( ElConstDistMatrix_z A, bool* viewing );

/* bool DistMatrix<T,U,V>::Locked() const
   --------------------------------------- */
ElError ElDistMatrixLocked_s( ElConstDistMatrix_s A, bool* locked );
ElError ElDistMatrixLocked_d( ElConstDistMatrix_d A, bool* locked );
ElError ElDistMatrixLocked_c( ElConstDistMatrix_c A, bool* locked );
ElError ElDistMatrixLocked_z( ElConstDistMatrix_z A, bool* locked );

/* Int DistMatrix<T,U,V>::LocalHeight() const
   ------------------------------------------ */
ElError ElDistMatrixLocalHeight_s( ElConstDistMatrix_s A, ElInt* localHeight );
ElError ElDistMatrixLocalHeight_d( ElConstDistMatrix_d A, ElInt* localHeight );
ElError ElDistMatrixLocalHeight_c( ElConstDistMatrix_c A, ElInt* localHeight );
ElError ElDistMatrixLocalHeight_z( ElConstDistMatrix_z A, ElInt* localHeight );

/* Int DistMatrix<T,U,V>::LocalWidth() const
   ----------------------------------------- */
ElError ElDistMatrixLocalWidth_s( ElConstDistMatrix_s A, ElInt* localWidth );
ElError ElDistMatrixLocalWidth_d( ElConstDistMatrix_d A, ElInt* localWidth );
ElError ElDistMatrixLocalWidth_c( ElConstDistMatrix_c A, ElInt* localWidth );
ElError ElDistMatrixLocalWidth_z( ElConstDistMatrix_z A, ElInt* localWidth );

/* Int DistMatrix<T,U,V>::LDim() const
   ----------------------------------- */
ElError ElDistMatrixLDim_s( ElConstDistMatrix_s A, ElInt* ldim );
ElError ElDistMatrixLDim_d( ElConstDistMatrix_d A, ElInt* ldim );
ElError ElDistMatrixLDim_c( ElConstDistMatrix_c A, ElInt* ldim );
ElError ElDistMatrixLDim_z( ElConstDistMatrix_z A, ElInt* ldim );

/* Matrix<T>& DistMatrix<T,U,V>::Matrix() 
   -------------------------------------- */
ElError ElDistMatrixMatrix_s( ElDistMatrix_s A, ElMatrix_s* ALoc );
ElError ElDistMatrixMatrix_d( ElDistMatrix_d A, ElMatrix_d* ALoc );
ElError ElDistMatrixMatrix_c( ElDistMatrix_c A, ElMatrix_c* ALoc );
ElError ElDistMatrixMatrix_z( ElDistMatrix_z A, ElMatrix_z* ALoc );

/* const Matrix<T>& DistMatrix<T,U,V>::LockedMatrix() const
   -------------------------------------------------------- */
ElError ElDistMatrixLockedMatrix_s
( ElConstDistMatrix_s A, ElConstMatrix_s* ALoc );
ElError ElDistMatrixLockedMatrix_d
( ElConstDistMatrix_d A, ElConstMatrix_d* ALoc );
ElError ElDistMatrixLockedMatrix_c
( ElConstDistMatrix_c A, ElConstMatrix_c* ALoc );
ElError ElDistMatrixLockedMatrix_z
( ElConstDistMatrix_z A, ElConstMatrix_z* ALoc );

/* size_t DistMatrix<T,U,V>::AllocatedMemory() const
   ------------------------------------------------- */
ElError ElDistMatrixAllocatedMemory_s( ElConstDistMatrix_s A, size_t* mem );
ElError ElDistMatrixAllocatedMemory_d( ElConstDistMatrix_d A, size_t* mem );
ElError ElDistMatrixAllocatedMemory_c( ElConstDistMatrix_c A, size_t* mem );
ElError ElDistMatrixAllocatedMemory_z( ElConstDistMatrix_z A, size_t* mem );

/* T* DistMatrix<T,U,V>::Buffer()
   ------------------------------ */
ElError ElDistMatrixBuffer_s( ElDistMatrix_s A, float** buffer );
ElError ElDistMatrixBuffer_d( ElDistMatrix_d A, double** buffer );
ElError ElDistMatrixBuffer_c( ElDistMatrix_c A, complex_float** buffer );
ElError ElDistMatrixBuffer_z( ElDistMatrix_z A, complex_double** buffer );

/* const T* DistMatrix<T,U,V>::LockedBuffer() const
   ------------------------------------------------ */
ElError ElDistMatrixLockedBuffer_s
( ElConstDistMatrix_s A, const float** buffer );
ElError ElDistMatrixLockedBuffer_d
( ElConstDistMatrix_d A, const double** buffer );
ElError ElDistMatrixLockedBuffer_c
( ElConstDistMatrix_c A, const complex_float** buffer );
ElError ElDistMatrixLockedBuffer_z
( ElConstDistMatrix_z A, const complex_double** buffer );

/* const Grid& DistMatrix<T,U,V>::Grid() const
   ------------------------------------------- */
ElError ElDistMatrixGrid_s( ElConstDistMatrix_s A, ElConstGrid* grid );
ElError ElDistMatrixGrid_d( ElConstDistMatrix_d A, ElConstGrid* grid );
ElError ElDistMatrixGrid_c( ElConstDistMatrix_c A, ElConstGrid* grid );
ElError ElDistMatrixGrid_z( ElConstDistMatrix_z A, ElConstGrid* grid );

/* bool DistMatrix<T,U,V>::ColConstrained() const
   ---------------------------------------------- */
ElError ElDistMatrixColConstrained_s( ElConstDistMatrix_s A, bool* colConst );
ElError ElDistMatrixColConstrained_d( ElConstDistMatrix_d A, bool* colConst );
ElError ElDistMatrixColConstrained_c( ElConstDistMatrix_c A, bool* colConst );
ElError ElDistMatrixColConstrained_z( ElConstDistMatrix_z A, bool* colConst );

/* bool DistMatrix<T,U,V>::RowConstrained() const
   ---------------------------------------------- */
ElError ElDistMatrixRowConstrained_s( ElConstDistMatrix_s A, bool* rowConst );
ElError ElDistMatrixRowConstrained_d( ElConstDistMatrix_d A, bool* rowConst );
ElError ElDistMatrixRowConstrained_c( ElConstDistMatrix_c A, bool* rowConst );
ElError ElDistMatrixRowConstrained_z( ElConstDistMatrix_z A, bool* rowConst );

/* bool DistMatrix<T,U,V>::RootConstrained() const
   ----------------------------------------------- */
ElError ElDistMatrixRootConstrained_s( ElConstDistMatrix_s A, bool* rootConst );
ElError ElDistMatrixRootConstrained_d( ElConstDistMatrix_d A, bool* rootConst );
ElError ElDistMatrixRootConstrained_c( ElConstDistMatrix_c A, bool* rootConst );
ElError ElDistMatrixRootConstrained_z( ElConstDistMatrix_z A, bool* rootConst );

/* Int DistMatrix<T,U,V>::ColAlign() const
   --------------------------------------- */
ElError ElDistMatrixColAlign_s( ElConstDistMatrix_s A, ElInt* colAlign );
ElError ElDistMatrixColAlign_d( ElConstDistMatrix_d A, ElInt* colAlign );
ElError ElDistMatrixColAlign_c( ElConstDistMatrix_c A, ElInt* colAlign );
ElError ElDistMatrixColAlign_z( ElConstDistMatrix_z A, ElInt* colAlign );

/* Int DistMatrix<T,U,V>::RowAlign() const
   --------------------------------------- */
ElError ElDistMatrixRowAlign_s( ElConstDistMatrix_s A, ElInt* rowAlign );
ElError ElDistMatrixRowAlign_d( ElConstDistMatrix_d A, ElInt* rowAlign );
ElError ElDistMatrixRowAlign_c( ElConstDistMatrix_c A, ElInt* rowAlign );
ElError ElDistMatrixRowAlign_z( ElConstDistMatrix_z A, ElInt* rowAlign );

/* Int DistMatrix<T,U,V>::ColShift() const
   --------------------------------------- */
ElError ElDistMatrixColShift_s( ElConstDistMatrix_s A, ElInt* colShift );
ElError ElDistMatrixColShift_d( ElConstDistMatrix_d A, ElInt* colShift );
ElError ElDistMatrixColShift_c( ElConstDistMatrix_c A, ElInt* colShift );
ElError ElDistMatrixColShift_z( ElConstDistMatrix_z A, ElInt* colShift );

/* Int DistMatrix<T,U,V>::RowShift() const
   --------------------------------------- */
ElError ElDistMatrixRowShift_s( ElConstDistMatrix_s A, ElInt* rowShift );
ElError ElDistMatrixRowShift_d( ElConstDistMatrix_d A, ElInt* rowShift );
ElError ElDistMatrixRowShift_c( ElConstDistMatrix_c A, ElInt* rowShift );
ElError ElDistMatrixRowShift_z( ElConstDistMatrix_z A, ElInt* rowShift );

/* Int DistMatrix<T,U,V>::ColRank() const
   -------------------------------------- */
ElError ElDistMatrixColRank_s( ElConstDistMatrix_s A, ElInt* rank );
ElError ElDistMatrixColRank_d( ElConstDistMatrix_d A, ElInt* rank );
ElError ElDistMatrixColRank_c( ElConstDistMatrix_c A, ElInt* rank );
ElError ElDistMatrixColRank_z( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::RowRank() const
   -------------------------------------- */
ElError ElDistMatrixRowRank_s( ElConstDistMatrix_s A, ElInt* rank );
ElError ElDistMatrixRowRank_d( ElConstDistMatrix_d A, ElInt* rank );
ElError ElDistMatrixRowRank_c( ElConstDistMatrix_c A, ElInt* rank );
ElError ElDistMatrixRowRank_z( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::PartialColRank() const
   --------------------------------------------- */
ElError ElDistMatrixPartialColRank_s( ElConstDistMatrix_s A, ElInt* rank );
ElError ElDistMatrixPartialColRank_d( ElConstDistMatrix_d A, ElInt* rank );
ElError ElDistMatrixPartialColRank_c( ElConstDistMatrix_c A, ElInt* rank );
ElError ElDistMatrixPartialColRank_z( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::PartialRowRank() const
   --------------------------------------------- */
ElError ElDistMatrixPartialRowRank_s( ElConstDistMatrix_s A, ElInt* rank );
ElError ElDistMatrixPartialRowRank_d( ElConstDistMatrix_d A, ElInt* rank );
ElError ElDistMatrixPartialRowRank_c( ElConstDistMatrix_c A, ElInt* rank );
ElError ElDistMatrixPartialRowRank_z( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::PartialUnionColRank() const
   -------------------------------------------------- */
ElError ElDistMatrixPartialUnionColRank_s( ElConstDistMatrix_s A, ElInt* rank );
ElError ElDistMatrixPartialUnionColRank_d( ElConstDistMatrix_d A, ElInt* rank );
ElError ElDistMatrixPartialUnionColRank_c( ElConstDistMatrix_c A, ElInt* rank );
ElError ElDistMatrixPartialUnionColRank_z( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::PartialUnionRowRank() const
   -------------------------------------------------- */
ElError ElDistMatrixPartialUnionRowRank_s( ElConstDistMatrix_s A, ElInt* rank );
ElError ElDistMatrixPartialUnionRowRank_d( ElConstDistMatrix_d A, ElInt* rank );
ElError ElDistMatrixPartialUnionRowRank_c( ElConstDistMatrix_c A, ElInt* rank );
ElError ElDistMatrixPartialUnionRowRank_z( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::DistRank() const
   --------------------------------------- */
ElError ElDistMatrixDistRank_s( ElConstDistMatrix_s A, ElInt* rank );
ElError ElDistMatrixDistRank_d( ElConstDistMatrix_d A, ElInt* rank );
ElError ElDistMatrixDistRank_c( ElConstDistMatrix_c A, ElInt* rank );
ElError ElDistMatrixDistRank_z( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::CrossRank() const
   ---------------------------------------- */
ElError ElDistMatrixCrossRank_s( ElConstDistMatrix_s A, ElInt* rank );
ElError ElDistMatrixCrossRank_d( ElConstDistMatrix_d A, ElInt* rank );
ElError ElDistMatrixCrossRank_c( ElConstDistMatrix_c A, ElInt* rank );
ElError ElDistMatrixCrossRank_z( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::RedundantRank() const
   -------------------------------------------- */
ElError ElDistMatrixRedundantRank_s( ElConstDistMatrix_s A, ElInt* rank );
ElError ElDistMatrixRedundantRank_d( ElConstDistMatrix_d A, ElInt* rank );
ElError ElDistMatrixRedundantRank_c( ElConstDistMatrix_c A, ElInt* rank );
ElError ElDistMatrixRedundantRank_z( ElConstDistMatrix_z A, ElInt* rank );

/* Int DistMatrix<T,U,V>::Root() const
   ----------------------------------- */
ElError ElDistMatrixRoot_s( ElConstDistMatrix_s A, ElInt* root );
ElError ElDistMatrixRoot_d( ElConstDistMatrix_d A, ElInt* root );
ElError ElDistMatrixRoot_c( ElConstDistMatrix_c A, ElInt* root );
ElError ElDistMatrixRoot_z( ElConstDistMatrix_z A, ElInt* root );

/* bool DistMatrix<T,U,V>::Participating() const
   --------------------------------------------- */
ElError ElDistMatrixParticipating_s
( ElConstDistMatrix_s A, bool* participating );
ElError ElDistMatrixParticipating_d
( ElConstDistMatrix_d A, bool* participating );
ElError ElDistMatrixParticipating_c
( ElConstDistMatrix_c A, bool* participating );
ElError ElDistMatrixParticipating_z
( ElConstDistMatrix_z A, bool* participating );

/* Int DistMatrix<T,U,V>::RowOwner( Int i ) const
   ---------------------------------------------- */
ElError ElDistMatrixRowOwner_s
( ElConstDistMatrix_s A, ElInt i, ElInt* rowOwner );
ElError ElDistMatrixRowOwner_d
( ElConstDistMatrix_d A, ElInt i, ElInt* rowOwner );
ElError ElDistMatrixRowOwner_c
( ElConstDistMatrix_c A, ElInt i, ElInt* rowOwner );
ElError ElDistMatrixRowOwner_z
( ElConstDistMatrix_z A, ElInt i, ElInt* rowOwner );

/* Int DistMatrix<T,U,V>::ColOwner( Int j ) const
   ---------------------------------------------- */
ElError ElDistMatrixColOwner_s
( ElConstDistMatrix_s A, ElInt j, ElInt* colOwner );
ElError ElDistMatrixColOwner_d
( ElConstDistMatrix_d A, ElInt j, ElInt* colOwner );
ElError ElDistMatrixColOwner_c
( ElConstDistMatrix_c A, ElInt j, ElInt* colOwner );
ElError ElDistMatrixColOwner_z
( ElConstDistMatrix_z A, ElInt j, ElInt* colOwner );

/* Int DistMatrix<T,U,V>::Owner( Int i, Int j ) const
   -------------------------------------------------- */
ElError ElDistMatrixOwner_s
( ElConstDistMatrix_s A, ElInt i, ElInt j, ElInt* owner );
ElError ElDistMatrixOwner_d
( ElConstDistMatrix_d A, ElInt i, ElInt j, ElInt* owner );
ElError ElDistMatrixOwner_c
( ElConstDistMatrix_c A, ElInt i, ElInt j, ElInt* owner );
ElError ElDistMatrixOwner_z
( ElConstDistMatrix_z A, ElInt i, ElInt j, ElInt* owner );

/* Int DistMatrix<T,U,V>::LocalRow( Int i ) const
   ---------------------------------------------- */
ElError ElDistMatrixLocalRow_s( ElConstDistMatrix_s A, ElInt i, ElInt* iLoc );
ElError ElDistMatrixLocalRow_d( ElConstDistMatrix_d A, ElInt i, ElInt* iLoc );
ElError ElDistMatrixLocalRow_c( ElConstDistMatrix_c A, ElInt i, ElInt* iLoc );
ElError ElDistMatrixLocalRow_z( ElConstDistMatrix_z A, ElInt i, ElInt* iLoc );

/* Int DistMatrix<T,U,V>::LocalCol( Int j ) const
   ---------------------------------------------- */
ElError ElDistMatrixLocalCol_s( ElConstDistMatrix_s A, ElInt j, ElInt* jLoc );
ElError ElDistMatrixLocalCol_d( ElConstDistMatrix_d A, ElInt j, ElInt* jLoc );
ElError ElDistMatrixLocalCol_c( ElConstDistMatrix_c A, ElInt j, ElInt* jLoc );
ElError ElDistMatrixLocalCol_z( ElConstDistMatrix_z A, ElInt j, ElInt* jLoc );

/* Int DistMatrix<T,U,V>::LocalRowOffset( Int i ) const
   ---------------------------------------------------- */
ElError ElDistMatrixLocalRowOffset_s
( ElConstDistMatrix_s A, ElInt i, ElInt* iLoc );
ElError ElDistMatrixLocalRowOffset_d
( ElConstDistMatrix_d A, ElInt i, ElInt* iLoc );
ElError ElDistMatrixLocalRowOffset_c
( ElConstDistMatrix_c A, ElInt i, ElInt* iLoc );
ElError ElDistMatrixLocalRowOffset_z
( ElConstDistMatrix_z A, ElInt i, ElInt* iLoc );

/* Int DistMatrix<T,U,V>::LocalColOffset( Int j ) const
   ---------------------------------------------------- */
ElError ElDistMatrixLocalColOffset_s
( ElConstDistMatrix_s A, ElInt j, ElInt* jLoc );
ElError ElDistMatrixLocalColOffset_d
( ElConstDistMatrix_d A, ElInt j, ElInt* jLoc );
ElError ElDistMatrixLocalColOffset_c
( ElConstDistMatrix_c A, ElInt j, ElInt* jLoc );
ElError ElDistMatrixLocalColOffset_z
( ElConstDistMatrix_z A, ElInt j, ElInt* jLoc );

/* Int DistMatrix<T,U,V>::GlobalRow( Int iLoc ) const
   -------------------------------------------------- */
ElError ElDistMatrixGlobalRow_s( ElConstDistMatrix_s A, ElInt iLoc, ElInt* i );
ElError ElDistMatrixGlobalRow_d( ElConstDistMatrix_d A, ElInt iLoc, ElInt* i );
ElError ElDistMatrixGlobalRow_c( ElConstDistMatrix_c A, ElInt iLoc, ElInt* i );
ElError ElDistMatrixGlobalRow_z( ElConstDistMatrix_z A, ElInt iLoc, ElInt* i );

/* Int DistMatrix<T,U,V>::GlobalCol( Int jLoc ) const
   -------------------------------------------------- */
ElError ElDistMatrixGlobalCol_s( ElConstDistMatrix_s A, ElInt jLoc, ElInt* j );
ElError ElDistMatrixGlobalCol_d( ElConstDistMatrix_d A, ElInt jLoc, ElInt* j );
ElError ElDistMatrixGlobalCol_c( ElConstDistMatrix_c A, ElInt jLoc, ElInt* j );
ElError ElDistMatrixGlobalCol_z( ElConstDistMatrix_z A, ElInt jLoc, ElInt* j );

/* bool DistMatrix<T,U,V>::IsLocalRow( Int i ) const
   ------------------------------------------------- */
ElError ElDistMatrixIsLocalRow_s
( ElConstDistMatrix_s A, ElInt i, bool* isLocal );
ElError ElDistMatrixIsLocalRow_d
( ElConstDistMatrix_d A, ElInt i, bool* isLocal );
ElError ElDistMatrixIsLocalRow_c
( ElConstDistMatrix_c A, ElInt i, bool* isLocal );
ElError ElDistMatrixIsLocalRow_z
( ElConstDistMatrix_z A, ElInt i, bool* isLocal );

/* bool DistMatrix<T,U,V>::IsLocalCol( Int j ) const
   ------------------------------------------------- */
ElError ElDistMatrixIsLocalCol_s
( ElConstDistMatrix_s A, ElInt j, bool* isLocal );
ElError ElDistMatrixIsLocalCol_d
( ElConstDistMatrix_d A, ElInt j, bool* isLocal );
ElError ElDistMatrixIsLocalCol_c
( ElConstDistMatrix_c A, ElInt j, bool* isLocal );
ElError ElDistMatrixIsLocalCol_z
( ElConstDistMatrix_z A, ElInt j, bool* isLocal );

/* bool DistMatrix<T,U,V>::IsLocal( Int i, Int j ) const
   ----------------------------------------------------- */
ElError ElDistMatrixIsLocal_s
( ElConstDistMatrix_s A, ElInt i, ElInt j, bool* isLocal );
ElError ElDistMatrixIsLocal_d
( ElConstDistMatrix_d A, ElInt i, ElInt j, bool* isLocal );
ElError ElDistMatrixIsLocal_c
( ElConstDistMatrix_c A, ElInt i, ElInt j, bool* isLocal );
ElError ElDistMatrixIsLocal_z
( ElConstDistMatrix_z A, ElInt i, ElInt j, bool* isLocal );

/* DistData DistMatrix<T,U,V>::DistData() const
   -------------------------------------------- */
ElError ElDistMatrixDistData_s( ElConstDistMatrix_s A, ElDistData* distData );
ElError ElDistMatrixDistData_d( ElConstDistMatrix_d A, ElDistData* distData );
ElError ElDistMatrixDistData_c( ElConstDistMatrix_c A, ElDistData* distData );
ElError ElDistMatrixDistData_z( ElConstDistMatrix_z A, ElDistData* distData );

/* TODO: A large number of routines here */

/* T DistMatrix<T,U,V>::Get( Int i, Int j ) const
   ---------------------------------------------- */
ElError ElDistMatrixGet_s
( ElConstDistMatrix_s A, ElInt i, ElInt j, float* val );
ElError ElDistMatrixGet_d
( ElConstDistMatrix_d A, ElInt i, ElInt j, double* val );
ElError ElDistMatrixGet_c
( ElConstDistMatrix_c A, ElInt i, ElInt j, complex_float* val );
ElError ElDistMatrixGet_z
( ElConstDistMatrix_z A, ElInt i, ElInt j, complex_double* val );

/* Base<T> DistMatrix<T,U,V>::GetRealPart( Int i, Int j ) const
   ------------------------------------------------------------ */
ElError ElDistMatrixGetRealPart_c
( ElConstDistMatrix_c A, ElInt i, ElInt j, float* val );
ElError ElDistMatrixGetRealPart_z
( ElConstDistMatrix_z A, ElInt i, ElInt j, double* val );

/* Base<T> DistMatrix<T,U,V>::GetImagPart( Int i, Int j ) const
   ------------------------------------------------------------ */
ElError ElDistMatrixGetImagPart_c
( ElConstDistMatrix_c A, ElInt i, ElInt j, float* val );
ElError ElDistMatrixGetImagPart_z
( ElConstDistMatrix_z A, ElInt i, ElInt j, double* val );

/* void DistMatrix<T,U,V>::Set( Int i, Int j, T alpha )
   ---------------------------------------------------- */
ElError ElDistMatrixSet_s
( ElDistMatrix_s A, ElInt i, ElInt j, float alpha );
ElError ElDistMatrixSet_d
( ElDistMatrix_d A, ElInt i, ElInt j, double alpha );
ElError ElDistMatrixSet_c
( ElDistMatrix_c A, ElInt i, ElInt j, complex_float alpha );
ElError ElDistMatrixSet_z
( ElDistMatrix_z A, ElInt i, ElInt j, complex_double alpha );

/* void DistMatrix<T,U,V>::SetRealPart( Int i, Int j, Base<T> alpha )
   ------------------------------------------------------------------ */
ElError ElDistMatrixSetRealPart_c
( ElDistMatrix_c A, ElInt i, ElInt j, float alpha );
ElError ElDistMatrixSetRealPart_z
( ElDistMatrix_z A, ElInt i, ElInt j, double alpha );

/* void DistMatrix<T,U,V>::SetImagPart( Int i, Int j, Base<T> alpha )
   ------------------------------------------------------------------ */
ElError ElDistMatrixSetImagPart_c
( ElDistMatrix_c A, ElInt i, ElInt j, float alpha );
ElError ElDistMatrixSetImagPart_z
( ElDistMatrix_z A, ElInt i, ElInt j, double alpha );

/* void DistMatrix<T,U,V>::Update( Int i, Int j, T alpha )
   ------------------------------------------------------- */
ElError ElDistMatrixUpdate_s
( ElDistMatrix_s A, ElInt i, ElInt j, float alpha );
ElError ElDistMatrixUpdate_d
( ElDistMatrix_d A, ElInt i, ElInt j, double alpha );
ElError ElDistMatrixUpdate_c
( ElDistMatrix_c A, ElInt i, ElInt j, complex_float alpha );
ElError ElDistMatrixUpdate_z
( ElDistMatrix_z A, ElInt i, ElInt j, complex_double alpha );

/* void DistMatrix<T,U,V>::UpdateRealPart( Int i, Int j, Base<T> alpha )
   --------------------------------------------------------------------- */
ElError ElDistMatrixUpdateRealPart_c
( ElDistMatrix_c A, ElInt i, ElInt j, float alpha );
ElError ElDistMatrixUpdateRealPart_z
( ElDistMatrix_z A, ElInt i, ElInt j, double alpha );

/* void DistMatrix<T,U,V>::UpdateImagPart( Int i, Int j, Base<T> alpha )
   --------------------------------------------------------------------- */
ElError ElDistMatrixUpdateImagPart_c
( ElDistMatrix_c A, ElInt i, ElInt j, float alpha );
ElError ElDistMatrixUpdateImagPart_z
( ElDistMatrix_z A, ElInt i, ElInt j, double alpha );

/* void DistMatrix<T,U,V>::MakeReal( Int i, Int j )
   ------------------------------------------------ */
ElError ElDistMatrixMakeReal_c( ElDistMatrix_c A, ElInt i, ElInt j );
ElError ElDistMatrixMakeReal_z( ElDistMatrix_z A, ElInt i, ElInt j );

/* void DistMatrix<T,U,V>::Conjugate( Int i, Int j )
   ------------------------------------------------ */
ElError ElDistMatrixConjugate_c( ElDistMatrix_c A, ElInt i, ElInt j );
ElError ElDistMatrixConjugate_z( ElDistMatrix_z A, ElInt i, ElInt j );

/* DistMatrix<T,UDiag,VDiag> DistMatrix<T,U,V>::GetDiagonal( Int offset ) const
   ----------------------------------------------------------------------------
*/
ElError ElDistMatrixGetDiagonal_s
( ElConstDistMatrix_s A, ElInt offset, ElDistMatrix_s* d );
ElError ElDistMatrixGetDiagonal_d
( ElConstDistMatrix_d A, ElInt offset, ElDistMatrix_d* d );
ElError ElDistMatrixGetDiagonal_c
( ElConstDistMatrix_c A, ElInt offset, ElDistMatrix_c* d );
ElError ElDistMatrixGetDiagonal_z
( ElConstDistMatrix_z A, ElInt offset, ElDistMatrix_z* d );

/* TODO: More diagonal manipulation
   ================================ */

/* DistMatrix<T,STAR,STAR> DistMatrix<T,U,V>::GetSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
   -------------------------------------------------------------------------- */
ElError ElDistMatrixGetSubmatrix_s
( ElConstDistMatrix_s A,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_s* ASub );
ElError ElDistMatrixGetSubmatrix_d
( ElConstDistMatrix_d A,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_d* ASub );
ElError ElDistMatrixGetSubmatrix_c
( ElConstDistMatrix_c A,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_c* ASub );
ElError ElDistMatrixGetSubmatrix_z
( ElConstDistMatrix_z A,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_z* ASub );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_DISTMATRIX_C_H */
