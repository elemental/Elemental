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

/* mpi::Comm DistMatrix<T,U,V>::DistComm() const
   --------------------------------------------- */
ElError ElDistMatrixDistComm_s( ElConstDistMatrix_s A, MPI_Comm* comm );
ElError ElDistMatrixDistComm_d( ElConstDistMatrix_d A, MPI_Comm* comm );
ElError ElDistMatrixDistComm_c( ElConstDistMatrix_c A, MPI_Comm* comm );
ElError ElDistMatrixDistComm_z( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::CrossComm() const
   ---------------------------------------------- */
ElError ElDistMatrixCrossComm_s( ElConstDistMatrix_s A, MPI_Comm* comm );
ElError ElDistMatrixCrossComm_d( ElConstDistMatrix_d A, MPI_Comm* comm );
ElError ElDistMatrixCrossComm_c( ElConstDistMatrix_c A, MPI_Comm* comm );
ElError ElDistMatrixCrossComm_z( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::RedundantComm() const
   -------------------------------------------------- */
ElError ElDistMatrixRedundantComm_s( ElConstDistMatrix_s A, MPI_Comm* comm );
ElError ElDistMatrixRedundantComm_d( ElConstDistMatrix_d A, MPI_Comm* comm );
ElError ElDistMatrixRedundantComm_c( ElConstDistMatrix_c A, MPI_Comm* comm );
ElError ElDistMatrixRedundantComm_z( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::ColComm() const
   -------------------------------------------- */
ElError ElDistMatrixColComm_s( ElConstDistMatrix_s A, MPI_Comm* comm );
ElError ElDistMatrixColComm_d( ElConstDistMatrix_d A, MPI_Comm* comm );
ElError ElDistMatrixColComm_c( ElConstDistMatrix_c A, MPI_Comm* comm );
ElError ElDistMatrixColComm_z( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::RowComm() const
   -------------------------------------------- */
ElError ElDistMatrixRowComm_s( ElConstDistMatrix_s A, MPI_Comm* comm );
ElError ElDistMatrixRowComm_d( ElConstDistMatrix_d A, MPI_Comm* comm );
ElError ElDistMatrixRowComm_c( ElConstDistMatrix_c A, MPI_Comm* comm );
ElError ElDistMatrixRowComm_z( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::PartialColComm() const
   --------------------------------------------------- */
ElError ElDistMatrixPartialColComm_s( ElConstDistMatrix_s A, MPI_Comm* comm );
ElError ElDistMatrixPartialColComm_d( ElConstDistMatrix_d A, MPI_Comm* comm );
ElError ElDistMatrixPartialColComm_c( ElConstDistMatrix_c A, MPI_Comm* comm );
ElError ElDistMatrixPartialColComm_z( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::PartialRowComm() const
   --------------------------------------------------- */
ElError ElDistMatrixPartialRowComm_s( ElConstDistMatrix_s A, MPI_Comm* comm );
ElError ElDistMatrixPartialRowComm_d( ElConstDistMatrix_d A, MPI_Comm* comm );
ElError ElDistMatrixPartialRowComm_c( ElConstDistMatrix_c A, MPI_Comm* comm );
ElError ElDistMatrixPartialRowComm_z( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::PartialUnionColComm() const
   -------------------------------------------------------- */
ElError ElDistMatrixPartialUnionColComm_s
( ElConstDistMatrix_s A, MPI_Comm* comm );
ElError ElDistMatrixPartialUnionColComm_d
( ElConstDistMatrix_d A, MPI_Comm* comm );
ElError ElDistMatrixPartialUnionColComm_c
( ElConstDistMatrix_c A, MPI_Comm* comm );
ElError ElDistMatrixPartialUnionColComm_z
( ElConstDistMatrix_z A, MPI_Comm* comm );

/* mpi::Comm DistMatrix<T,U,V>::PartialUnionRowComm() const
   -------------------------------------------------------- */
ElError ElDistMatrixPartialUnionRowComm_s
( ElConstDistMatrix_s A, MPI_Comm* comm );
ElError ElDistMatrixPartialUnionRowComm_d
( ElConstDistMatrix_d A, MPI_Comm* comm );
ElError ElDistMatrixPartialUnionRowComm_c
( ElConstDistMatrix_c A, MPI_Comm* comm );
ElError ElDistMatrixPartialUnionRowComm_z
( ElConstDistMatrix_z A, MPI_Comm* comm );

/* Int DistMatrix<T,U,V>::ColStride() const
   ---------------------------------------- */
ElError ElDistMatrixColStride_s( ElConstDistMatrix_s AHandle, ElInt* stride );
ElError ElDistMatrixColStride_d( ElConstDistMatrix_d AHandle, ElInt* stride );
ElError ElDistMatrixColStride_c( ElConstDistMatrix_c AHandle, ElInt* stride );
ElError ElDistMatrixColStride_z( ElConstDistMatrix_z AHandle, ElInt* stride );

/* Int DistMatrix<T,U,V>::RowStride() const
   ---------------------------------------- */
ElError ElDistMatrixRowStride_s( ElConstDistMatrix_s AHandle, ElInt* stride );
ElError ElDistMatrixRowStride_d( ElConstDistMatrix_d AHandle, ElInt* stride );
ElError ElDistMatrixRowStride_c( ElConstDistMatrix_c AHandle, ElInt* stride );
ElError ElDistMatrixRowStride_z( ElConstDistMatrix_z AHandle, ElInt* stride );

/* Int DistMatrix<T,U,V>::PartialColStride() const
   ----------------------------------------------- */
ElError ElDistMatrixPartialColStride_s
( ElConstDistMatrix_s AHandle, ElInt* stride );
ElError ElDistMatrixPartialColStride_d
( ElConstDistMatrix_d AHandle, ElInt* stride );
ElError ElDistMatrixPartialColStride_c
( ElConstDistMatrix_c AHandle, ElInt* stride );
ElError ElDistMatrixPartialColStride_z
( ElConstDistMatrix_z AHandle, ElInt* stride );

/* Int DistMatrix<T,U,V>::PartialRowStride() const
   ----------------------------------------------- */
ElError ElDistMatrixPartialRowStride_s
( ElConstDistMatrix_s AHandle, ElInt* stride );
ElError ElDistMatrixPartialRowStride_d
( ElConstDistMatrix_d AHandle, ElInt* stride );
ElError ElDistMatrixPartialRowStride_c
( ElConstDistMatrix_c AHandle, ElInt* stride );
ElError ElDistMatrixPartialRowStride_z
( ElConstDistMatrix_z AHandle, ElInt* stride );

/* Int DistMatrix<T,U,V>::PartialUnionColStride() const
   ---------------------------------------------------- */
ElError ElDistMatrixPartialUnionColStride_s
( ElConstDistMatrix_s AHandle, ElInt* stride );
ElError ElDistMatrixPartialUnionColStride_d
( ElConstDistMatrix_d AHandle, ElInt* stride );
ElError ElDistMatrixPartialUnionColStride_c
( ElConstDistMatrix_c AHandle, ElInt* stride );
ElError ElDistMatrixPartialUnionColStride_z
( ElConstDistMatrix_z AHandle, ElInt* stride );

/* Int DistMatrix<T,U,V>::PartialUnionRowStride() const
   ---------------------------------------------------- */
ElError ElDistMatrixPartialUnionRowStride_s
( ElConstDistMatrix_s AHandle, ElInt* stride );
ElError ElDistMatrixPartialUnionRowStride_d
( ElConstDistMatrix_d AHandle, ElInt* stride );
ElError ElDistMatrixPartialUnionRowStride_c
( ElConstDistMatrix_c AHandle, ElInt* stride );
ElError ElDistMatrixPartialUnionRowStride_z
( ElConstDistMatrix_z AHandle, ElInt* stride );

/* Int DistMatrix<T,U,V>::DistSize() const
   --------------------------------------- */
ElError ElDistMatrixDistSize_s( ElConstDistMatrix_s AHandle, ElInt* commSize );
ElError ElDistMatrixDistSize_d( ElConstDistMatrix_d AHandle, ElInt* commSize );
ElError ElDistMatrixDistSize_c( ElConstDistMatrix_c AHandle, ElInt* commSize );
ElError ElDistMatrixDistSize_z( ElConstDistMatrix_z AHandle, ElInt* commSize );

/* Int DistMatrix<T,U,V>::CrossSize() const
   ---------------------------------------- */
ElError ElDistMatrixCrossSize_s( ElConstDistMatrix_s AHandle, ElInt* commSize );
ElError ElDistMatrixCrossSize_d( ElConstDistMatrix_d AHandle, ElInt* commSize );
ElError ElDistMatrixCrossSize_c( ElConstDistMatrix_c AHandle, ElInt* commSize );
ElError ElDistMatrixCrossSize_z( ElConstDistMatrix_z AHandle, ElInt* commSize );

/* Int DistMatrix<T,U,V>::RedundantSize() const
   -------------------------------------------- */
ElError ElDistMatrixRedundantSize_s
( ElConstDistMatrix_s AHandle, ElInt* commSize );
ElError ElDistMatrixRedundantSize_d
( ElConstDistMatrix_d AHandle, ElInt* commSize );
ElError ElDistMatrixRedundantSize_c
( ElConstDistMatrix_c AHandle, ElInt* commSize );
ElError ElDistMatrixRedundantSize_z
( ElConstDistMatrix_z AHandle, ElInt* commSize );

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

/* T DistMatrix<T,U,V>::GetLocal( Int iLoc, Int jLoc ) const
   --------------------------------------------------------- */
ElError ElDistMatrixGetLocal_s
( ElConstDistMatrix_s A, ElInt iLoc, ElInt jLoc, float* val );
ElError ElDistMatrixGetLocal_d
( ElConstDistMatrix_d A, ElInt iLoc, ElInt jLoc, double* val );
ElError ElDistMatrixGetLocal_c
( ElConstDistMatrix_c A, ElInt iLoc, ElInt jLoc, complex_float* val );
ElError ElDistMatrixGetLocal_z
( ElConstDistMatrix_z A, ElInt iLoc, ElInt jLoc, complex_double* val );

/* Base<T> DistMatrix<T,U,V>::GetLocalRealPart( Int iLoc, Int jLoc ) const
   ----------------------------------------------------------------------- */
ElError ElDistMatrixGetLocalRealPart_c
( ElConstDistMatrix_c A, ElInt iLoc, ElInt jLoc, float* val );
ElError ElDistMatrixGetLocalRealPart_z
( ElConstDistMatrix_z A, ElInt iLoc, ElInt jLoc, double* val );

/* Base<T> DistMatrix<T,U,V>::GetLocalImagPart( Int iLoc, Int jLoc ) const
   ----------------------------------------------------------------------- */
ElError ElDistMatrixGetLocalImagPart_c
( ElConstDistMatrix_c A, ElInt iLoc, ElInt jLoc, float* val );
ElError ElDistMatrixGetLocalImagPart_z
( ElConstDistMatrix_z A, ElInt iLoc, ElInt jLoc, double* val );

/* void DistMatrix<T,U,V>::SetLocal( Int iLoc, Int jLoc, T val )
   ------------------------------------------------------------- */
ElError ElDistMatrixSetLocal_s
( ElDistMatrix_s A, ElInt iLoc, ElInt jLoc, float val );
ElError ElDistMatrixSetLocal_d
( ElDistMatrix_d A, ElInt iLoc, ElInt jLoc, double val );
ElError ElDistMatrixSetLocal_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, complex_float val );
ElError ElDistMatrixSetLocal_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, complex_double val );

/* void DistMatrix<T,U,V>::SetLocalRealPart( Int iLoc, Int jLoc, Base<T> val )
   --------------------------------------------------------------------------- 
*/
ElError ElDistMatrixSetLocalRealPart_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, float val );
ElError ElDistMatrixSetLocalRealPart_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, double val );

/* void DistMatrix<T,U,V>::SetLocalImagPart( Int iLoc, Int jLoc, Base<T> val )
   --------------------------------------------------------------------------- 
*/
ElError ElDistMatrixSetLocalImagPart_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, float val );
ElError ElDistMatrixSetLocalImagPart_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, double val );

/* void DistMatrix<T,U,V>::UpdateLocal( Int iLoc, Int jLoc, T val )
   ---------------------------------------------------------------- */
ElError ElDistMatrixUpdateLocal_s
( ElDistMatrix_s A, ElInt iLoc, ElInt jLoc, float val );
ElError ElDistMatrixUpdateLocal_d
( ElDistMatrix_d A, ElInt iLoc, ElInt jLoc, double val );
ElError ElDistMatrixUpdateLocal_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, complex_float val );
ElError ElDistMatrixUpdateLocal_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, complex_double val );

/* void DistMatrix<T,U,V>::UpdateLocalRealPart
   ( Int iLoc, Int jLoc, Base<T> val )
   ------------------------------------------- 
*/
ElError ElDistMatrixUpdateLocalRealPart_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, float val );
ElError ElDistMatrixUpdateLocalRealPart_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, double val );

/* void DistMatrix<T,U,V>::UpdateLocalImagPart
   ( Int iLoc, Int jLoc, Base<T> val )
   ------------------------------------------- 
*/
ElError ElDistMatrixUpdateLocalImagPart_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, float val );
ElError ElDistMatrixUpdateLocalImagPart_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, double val );

/* void DistMatrix<T,U,V>::MakeDiagonalReal( Int offset )
   ------------------------------------------------------ */
ElError ElDistMatrixMakeDiagonalReal_c( ElDistMatrix_c A, ElInt offset );
ElError ElDistMatrixMakeDiagonalReal_z( ElDistMatrix_z A, ElInt offset );

/* void DistMatrix<T,U,V>::ConjugateDiagonal( Int offset )
   ------------------------------------------------------- */
ElError ElDistMatrixConjugateDiagonal_c( ElDistMatrix_c A, ElInt offset );
ElError ElDistMatrixConjugateDiagonal_z( ElDistMatrix_z A, ElInt offset );

/* bool DistMatrix<T,U,V>::DiagonalAlignedWith
   ( const DistData& data, Int offset ) const
   ------------------------------------------- */
ElError ElDistMatrixAlignedWith_s
( ElConstDistMatrix_s A, ElDistData distData, ElInt offset, bool* aligned );
ElError ElDistMatrixAlignedWith_d
( ElConstDistMatrix_d A, ElDistData distData, ElInt offset, bool* aligned );
ElError ElDistMatrixAlignedWith_c
( ElConstDistMatrix_c A, ElDistData distData, ElInt offset, bool* aligned );
ElError ElDistMatrixAlignedWith_z
( ElConstDistMatrix_z A, ElDistData distData, ElInt offset, bool* aligned );

/* Int DistMatrix<T,U,V>::DiagonalRoot( Int offset ) const
   ------------------------------------------------------- */
ElError ElDistMatrixDiagonalRoot_s
( ElConstDistMatrix_s A, ElInt offset, ElInt* root );
ElError ElDistMatrixDiagonalRoot_d
( ElConstDistMatrix_d A, ElInt offset, ElInt* root );
ElError ElDistMatrixDiagonalRoot_c
( ElConstDistMatrix_c A, ElInt offset, ElInt* root );
ElError ElDistMatrixDiagonalRoot_z
( ElConstDistMatrix_z A, ElInt offset, ElInt* root );

/* Int DistMatrix<T,U,V>::DiagonalAlign( Int offset ) const
   -------------------------------------------------------- */
ElError ElDistMatrixDiagonalRoot_s
( ElConstDistMatrix_s A, ElInt offset, ElInt* align );
ElError ElDistMatrixDiagonalRoot_d
( ElConstDistMatrix_d A, ElInt offset, ElInt* align );
ElError ElDistMatrixDiagonalRoot_c
( ElConstDistMatrix_c A, ElInt offset, ElInt* align );
ElError ElDistMatrixDiagonalRoot_z
( ElConstDistMatrix_z A, ElInt offset, ElInt* align );

/* DistMatrix<T,UDiag,VDiag> 
   DistMatrix<T,U,V>::GetDiagonal( Int offset ) const
   -------------------------------------------------- */
ElError ElDistMatrixGetDiagonal_s
( ElConstDistMatrix_s A, ElInt offset, ElDistMatrix_s* d );
ElError ElDistMatrixGetDiagonal_d
( ElConstDistMatrix_d A, ElInt offset, ElDistMatrix_d* d );
ElError ElDistMatrixGetDiagonal_c
( ElConstDistMatrix_c A, ElInt offset, ElDistMatrix_c* d );
ElError ElDistMatrixGetDiagonal_z
( ElConstDistMatrix_z A, ElInt offset, ElDistMatrix_z* d );

/* DistMatrix<Base<T>,UDiag,VDiag> 
   DistMatrix<T,U,V>::GetRealPartOfDiagonal( Int offset ) const
   ------------------------------------------------------------ */
ElError ElDistMatrixGetRealPartOfDiagonal_c
( ElConstDistMatrix_c A, ElInt offset, ElDistMatrix_s* d );
ElError ElDistMatrixGetRealPartOfDiagonal_z
( ElConstDistMatrix_z A, ElInt offset, ElDistMatrix_d* d );

/* DistMatrix<Base<T>,UDiag,VDiag> 
   DistMatrix<T,U,V>::GetImagPartOfDiagonal( Int offset ) const
   ------------------------------------------------------------ */
ElError ElDistMatrixGetImagPartOfDiagonal_c
( ElConstDistMatrix_c A, ElInt offset, ElDistMatrix_s* d );
ElError ElDistMatrixGetImagPartOfDiagonal_z
( ElConstDistMatrix_z A, ElInt offset, ElDistMatrix_d* d );

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

/* DistMatrix<Base<T>,STAR,STAR> DistMatrix<T,U,V>::GetRealPartOfSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
   -------------------------------------------------------------------------- */
ElError ElDistMatrixGetRealPartOfSubmatrix_c
( ElConstDistMatrix_c A,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_s* ASub );
ElError ElDistMatrixGetRealPartOfSubmatrix_z
( ElConstDistMatrix_z A,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_d* ASub );

/* DistMatrix<Base<T>,STAR,STAR> DistMatrix<T,U,V>::GetImagPartOfSubmatrix
   ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
   -------------------------------------------------------------------------- */
ElError ElDistMatrixGetImagPartOfSubmatrix_c
( ElConstDistMatrix_c A,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_s* ASub );
ElError ElDistMatrixGetImagPartOfSubmatrix_z
( ElConstDistMatrix_z A,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_d* ASub );

/* void DistMatrix<T,U,V>::SetSubmatrix
   ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
     const DistMatrix<T,STAR,STAR>& ASub )
   ----------------------------------------------------------------- */
ElError ElDistMatrixSetSubmatrix_s
( ElDistMatrix_s A, const ElInt* rowInds, const ElInt* colInds, 
  ElConstDistMatrix_s ASub );
ElError ElDistMatrixSetSubmatrix_d
( ElDistMatrix_d A, const ElInt* rowInds, const ElInt* colInds, 
  ElConstDistMatrix_d ASub );
ElError ElDistMatrixSetSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowInds, const ElInt* colInds, 
  ElConstDistMatrix_c ASub );
ElError ElDistMatrixSetSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowInds, const ElInt* colInds, 
  ElConstDistMatrix_z ASub );

/* void DistMatrix<T,U,V>::SetRealPartOfSubmatrix
   ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
     const DistMatrix<Base<T>,STAR,STAR>& ASub )
   ----------------------------------------------------------------- */
ElError ElDistMatrixSetRealPartOfSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowInds, const ElInt* colInds, 
  ElConstDistMatrix_s ASub );
ElError ElDistMatrixSetRealPartOfSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowInds, const ElInt* colInds, 
  ElConstDistMatrix_d ASub );

/* void DistMatrix<T,U,V>::SetImagPartOfSubmatrix
   ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
     const DistMatrix<Base<T>,STAR,STAR>& ASub )
   ----------------------------------------------------------------- */
ElError ElDistMatrixSetImagPartOfSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowInds, const ElInt* colInds, 
  ElConstDistMatrix_s ASub );
ElError ElDistMatrixSetImagPartOfSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowInds, const ElInt* colInds, 
  ElConstDistMatrix_d ASub );

/* void DistMatrix<T,U,V>::UpdateSubmatrix
   ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
     T alpha, const DistMatrix<T,STAR,STAR>& ASub )
   ----------------------------------------------------------------- */
ElError ElDistMatrixUpdateSubmatrix_s
( ElDistMatrix_s A, const ElInt* rowInds, const ElInt* colInds, 
  float alpha, ElConstDistMatrix_s ASub );
ElError ElDistMatrixUpdateSubmatrix_d
( ElDistMatrix_d A, const ElInt* rowInds, const ElInt* colInds, 
  double alpha, ElConstDistMatrix_d ASub );
ElError ElDistMatrixUpdateSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowInds, const ElInt* colInds, 
  complex_float alpha, ElConstDistMatrix_c ASub );
ElError ElDistMatrixUpdateSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowInds, const ElInt* colInds, 
  complex_double alpha, ElConstDistMatrix_z ASub );

/* void DistMatrix<T,U,V>::UpdateRealPartOfSubmatrix
   ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
     Base<T> alpha, const DistMatrix<Base<T>,STAR,STAR>& ASub )
   ----------------------------------------------------------------- */
ElError ElDistMatrixUpdateRealPartOfSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowInds, const ElInt* colInds, 
  float alpha, ElConstDistMatrix_s ASub );
ElError ElDistMatrixUpdateRealPartOfSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowInds, const ElInt* colInds, 
  double alpha, ElConstDistMatrix_d ASub );

/* void DistMatrix<T,U,V>::UpdateImagPartOfSubmatrix
   ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
     Base<T> alpha, const DistMatrix<Base<T>,STAR,STAR>& ASub )
   ----------------------------------------------------------------- */
ElError ElDistMatrixUpdateImagPartOfSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowInds, const ElInt* colInds, 
  float alpha, ElConstDistMatrix_s ASub );
ElError ElDistMatrixUpdateImagPartOfSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowInds, const ElInt* colInds, 
  double alpha, ElConstDistMatrix_d ASub );

/* void DistMatrix<T,U,V>::MakeSubmatrixReal
   ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
   ------------------------------------------------------------------ */
ElError ElDistMatrixMakeSubmatrixReal_c
( ElDistMatrix_c A, ElInt numRowInds, const ElInt* rowInds, 
                    ElInt numColInds, const ElInt* colInds );
ElError ElDistMatrixMakeSubmatrixReal_z
( ElDistMatrix_z A, ElInt numRowInds, const ElInt* rowInds, 
                    ElInt numColInds, const ElInt* colInds );

/* void DistMatrix<T,U,V>::ConjugateSubmatrix
   ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
   ------------------------------------------------------------------ */
ElError ElDistMatrixConjugateSubmatrix_c
( ElDistMatrix_c A, ElInt numRowInds, const ElInt* rowInds, 
                    ElInt numColInds, const ElInt* colInds );
ElError ElDistMatrixConjugateSubmatrix_z
( ElDistMatrix_z A, ElInt numRowInds, const ElInt* rowInds, 
                    ElInt numColInds, const ElInt* colInds );

/* Matrix<T> DistMatrix<T,U,V>::GetLocalSubmatrix
   ( const std::vector<Int>& rowIndsLoc, 
     const std::vector<Int>& colIndsLoc ) const
   ---------------------------------------------- */
ElError ElDistMatrixGetLocalSubmatrix_s
( ElConstDistMatrix_s A,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_s* ASub );
ElError ElDistMatrixGetLocalSubmatrix_d
( ElConstDistMatrix_d A,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_d* ASub );
ElError ElDistMatrixGetLocalSubmatrix_c
( ElConstDistMatrix_c A,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_c* ASub );
ElError ElDistMatrixGetLocalSubmatrix_z
( ElConstDistMatrix_z A,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_z* ASub );

/* Matrix<Base<T>> DistMatrix<T,U,V>::GetRealPartOfLocalSubmatrix
   ( const std::vector<Int>& rowIndsLoc, 
     const std::vector<Int>& colIndsLoc ) const
   -------------------------------------------------------------- */
ElError ElDistMatrixGetRealPartOfLocalSubmatrix_c
( ElConstDistMatrix_c A,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_s* ASub );
ElError ElDistMatrixGetRealPartOfLocalSubmatrix_z
( ElConstDistMatrix_z A,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_d* ASub );

/* Matrix<Base<T>> DistMatrix<T,U,V>::GetImagPartOfLocalSubmatrix
   ( const std::vector<Int>& rowIndsLoc, 
     const std::vector<Int>& colIndsLoc ) const
   -------------------------------------------------------------- */
ElError ElDistMatrixGetImagPartOfLocalSubmatrix_c
( ElConstDistMatrix_c A,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_s* ASub );
ElError ElDistMatrixGetImagPartOfLocalSubmatrix_z
( ElConstDistMatrix_z A,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_d* ASub );

/* void DistMatrix<T,U,V>::SetLocalSubmatrix
   ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc, 
     const Matrix<T>& ASub )
   ---------------------------------------------------------------------- */
ElError ElDistMatrixSetLocalSubmatrix_s
( ElDistMatrix_s A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  ElConstMatrix_s ASub );
ElError ElDistMatrixSetLocalSubmatrix_d
( ElDistMatrix_d A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  ElConstMatrix_d ASub );
ElError ElDistMatrixSetLocalSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  ElConstMatrix_c ASub );
ElError ElDistMatrixSetLocalSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  ElConstMatrix_z ASub );

/* void DistMatrix<T,U,V>::SetRealPartOfLocalSubmatrix
   ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc, 
     const Matrix<Base<T>>& ASub )
   ---------------------------------------------------------------------- */
ElError ElDistMatrixSetRealPartOfLocalSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  ElConstMatrix_s ASub );
ElError ElDistMatrixSetRealPartOfLocalSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  ElConstMatrix_d ASub );

/* void DistMatrix<T,U,V>::SetImagPartOfLocalSubmatrix
   ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc, 
     const Matrix<Base<T>>& ASub )
   ---------------------------------------------------------------------- */
ElError ElDistMatrixSetImagPartOfLocalSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  ElConstMatrix_s ASub );
ElError ElDistMatrixSetImagPartOfLocalSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  ElConstMatrix_d ASub );

/* void DistMatrix<T,U,V>::UpdateLocalSubmatrix
   ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc, 
     T alpha, const Matrix<T>& ASub )
   ----------------------------------------------------------------------- */
ElError ElDistMatrixUpdateLocalSubmatrix_s
( ElDistMatrix_s A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  float alpha, ElConstMatrix_s ASub );
ElError ElDistMatrixUpdateLocalSubmatrix_d
( ElDistMatrix_d A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  double alpha, ElConstMatrix_d ASub );
ElError ElDistMatrixUpdateLocalSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  complex_float alpha, ElConstMatrix_c ASub );
ElError ElDistMatrixUpdateLocalSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  complex_double alpha, ElConstMatrix_z ASub );

/* void DistMatrix<T,U,V>::UpdateRealPartOfLocalSubmatrix
   ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc, 
     Base<T> alpha, const Matrix<Base<T>>& ASub )
   ----------------------------------------------------------------------- */
ElError ElDistMatrixUpdateRealPartOfLocalSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  float alpha, ElConstMatrix_s ASub );
ElError ElDistMatrixUpdateRealPartOfLocalSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  double alpha, ElConstMatrix_d ASub );

/* void DistMatrix<T,U,V>::UpdateImagPartOfLocalSubmatrix
   ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc, 
     Base<T> alpha, const Matrix<Base<T>>& ASub )
   ----------------------------------------------------------------------- */
ElError ElDistMatrixUpdateImagPartOfLocalSubmatrix_c
( ElDistMatrix_c A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  float alpha, ElConstMatrix_s ASub );
ElError ElDistMatrixUpdateImagPartOfLocalSubmatrix_z
( ElDistMatrix_z A, const ElInt* rowIndsLoc, const ElInt* colIndsLoc, 
  double alpha, ElConstMatrix_d ASub );

/* void DistMatrix<T,U,V>::MakeLocalSubmatrixReal
   ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc )
   ------------------------------------------------------------------------ */
ElError ElDistMatrixMakeLocalSubmatrixReal_c
( ElDistMatrix_c A, ElInt numRowInds, const ElInt* rowIndsLoc, 
                    ElInt numColInds, const ElInt* colIndsLoc );
ElError ElDistMatrixMakeLocalSubmatrixReal_z
( ElDistMatrix_z A, ElInt numRowInds, const ElInt* rowIndsLoc, 
                    ElInt numColInds, const ElInt* colIndsLoc );

/* void DistMatrix<T,U,V>::ConjugateLocalSubmatrix
   ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc )
   ------------------------------------------------------------------------ */
ElError ElDistMatrixConjugateLocalSubmatrix_c
( ElDistMatrix_c A, ElInt numRowInds, const ElInt* rowIndsLoc, 
                    ElInt numColInds, const ElInt* colIndsLoc );
ElError ElDistMatrixConjugateLocalSubmatrix_z
( ElDistMatrix_z A, ElInt numRowInds, const ElInt* rowIndsLoc, 
                    ElInt numColInds, const ElInt* colIndsLoc );

/* void DistMatrix<T,U,V>::SumOver( mpi::Comm comm )
   ------------------------------------------------- */
ElError ElDistMatrixSumOver_s( ElDistMatrix_s A, MPI_Comm comm );
ElError ElDistMatrixSumOver_d( ElDistMatrix_d A, MPI_Comm comm );
ElError ElDistMatrixSumOver_c( ElDistMatrix_c A, MPI_Comm comm );
ElError ElDistMatrixSumOver_z( ElDistMatrix_z A, MPI_Comm comm );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_DISTMATRIX_C_H */
