/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
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
    ElInt blockHeight, blockWidth;
    int colAlign, rowAlign;
    ElInt colCut, rowCut;
    int root;
    ElConstGrid grid;
} ElDistData;

/* An anonymous struct meant as a placeholder for ElementalMatrix<T>
   ----------------------------------------------------------------- */
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

/* ElementalMatrix<T>::~ElementalMatrix()
   -------------------------------------- */
EL_EXPORT ElError ElDistMatrixDestroy_i( ElConstDistMatrix_i A );
EL_EXPORT ElError ElDistMatrixDestroy_s( ElConstDistMatrix_s A );
EL_EXPORT ElError ElDistMatrixDestroy_d( ElConstDistMatrix_d A );
EL_EXPORT ElError ElDistMatrixDestroy_c( ElConstDistMatrix_c A );
EL_EXPORT ElError ElDistMatrixDestroy_z( ElConstDistMatrix_z A );

/* void AbstractDistMatrix<T>::Empty( bool freeMemory=true )
   --------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixEmpty_i( ElDistMatrix_i A, bool freeMemory );
EL_EXPORT ElError ElDistMatrixEmpty_s( ElDistMatrix_s A, bool freeMemory );
EL_EXPORT ElError ElDistMatrixEmpty_d( ElDistMatrix_d A, bool freeMemory );
EL_EXPORT ElError ElDistMatrixEmpty_c( ElDistMatrix_c A, bool freeMemory );
EL_EXPORT ElError ElDistMatrixEmpty_z( ElDistMatrix_z A, bool freeMemory );

/* void AbstractDistMatrix<T>::EmptyData( bool freeMemory=true )
   ------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixEmptyData_i( ElDistMatrix_i A, bool freeMemory );
EL_EXPORT ElError ElDistMatrixEmptyData_s( ElDistMatrix_s A, bool freeMemory );
EL_EXPORT ElError ElDistMatrixEmptyData_d( ElDistMatrix_d A, bool freeMemory );
EL_EXPORT ElError ElDistMatrixEmptyData_c( ElDistMatrix_c A, bool freeMemory );
EL_EXPORT ElError ElDistMatrixEmptyData_z( ElDistMatrix_z A, bool freeMemory );

/* void AbstractDistMatrix<T>::SetGrid( const Grid& g )
   ---------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSetGrid_i( ElDistMatrix_i A, ElConstGrid grid );
EL_EXPORT ElError ElDistMatrixSetGrid_s( ElDistMatrix_s A, ElConstGrid grid );
EL_EXPORT ElError ElDistMatrixSetGrid_d( ElDistMatrix_d A, ElConstGrid grid );
EL_EXPORT ElError ElDistMatrixSetGrid_c( ElDistMatrix_c A, ElConstGrid grid );
EL_EXPORT ElError ElDistMatrixSetGrid_z( ElDistMatrix_z A, ElConstGrid grid );

/* void AbstractDistMatrix<T>::Resize( Int height, Int width )
   ----------------------------------------------------------- */
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

/* void ElementalMatrix<T>::Resize( Int height, Int width, Int ldim )
   ------------------------------------------------------------------ */
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

/* void ElementalMatrix<T>::MakeConsistent( bool includeViewers )
   -------------------------------------------------------------- */
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

/* void ElementalMatrix<T>::MakeSizeConsistent( bool includeViewers )
   ------------------------------------------------------------------ */
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

/* void ElementalMatrix<T>::Align
   ( int colAlign, int rowAlign, bool constrain )
   ---------------------------------------------- */
EL_EXPORT ElError ElDistMatrixAlign_i
( ElDistMatrix_i A, int colAlign, int rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlign_s
( ElDistMatrix_s A, int colAlign, int rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlign_d
( ElDistMatrix_d A, int colAlign, int rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlign_c
( ElDistMatrix_c A, int colAlign, int rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlign_z
( ElDistMatrix_z A, int colAlign, int rowAlign, bool constrain );

/* void ElementalMatrix<T>::AlignCols( int colAlign, bool constrain )
   ------------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixAlignCols_i
( ElDistMatrix_i A, int colAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignCols_s
( ElDistMatrix_s A, int colAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignCols_d
( ElDistMatrix_d A, int colAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignCols_c
( ElDistMatrix_c A, int colAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignCols_z
( ElDistMatrix_z A, int colAlign, bool constrain );

/* void ElementalMatrix<T>::AlignRows( int rowAlign, bool constrain )
   ------------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixAlignRows_i
( ElDistMatrix_i A, int rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRows_s
( ElDistMatrix_s A, int rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRows_d
( ElDistMatrix_d A, int rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRows_c
( ElDistMatrix_c A, int rowAlign, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRows_z
( ElDistMatrix_z A, int rowAlign, bool constrain );

/* void ElementalMatrix<T>::FreeAlignments()
   ----------------------------------------- */
EL_EXPORT ElError ElDistMatrixFreeAlignments_i( ElDistMatrix_i A );
EL_EXPORT ElError ElDistMatrixFreeAlignments_s( ElDistMatrix_s A );
EL_EXPORT ElError ElDistMatrixFreeAlignments_d( ElDistMatrix_d A );
EL_EXPORT ElError ElDistMatrixFreeAlignments_c( ElDistMatrix_c A );
EL_EXPORT ElError ElDistMatrixFreeAlignments_z( ElDistMatrix_z A );

/* void ElementalMatrix<T>::SetRoot( int root, bool constrain )
   ------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixSetRoot_i
( ElDistMatrix_i A, int root, bool constrain );
EL_EXPORT ElError ElDistMatrixSetRoot_s
( ElDistMatrix_s A, int root, bool constrain );
EL_EXPORT ElError ElDistMatrixSetRoot_d
( ElDistMatrix_d A, int root, bool constrain );
EL_EXPORT ElError ElDistMatrixSetRoot_c
( ElDistMatrix_c A, int root, bool constrain );
EL_EXPORT ElError ElDistMatrixSetRoot_z
( ElDistMatrix_z A, int root, bool constrain );

/* void ElementalMatrix<T>::AlignWith
   ( const DistData& data, bool constrain )
   ---------------------------------------- */
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

/* void ElementalMatrix<T>::AlignColsWith
   ( const DistData& data, bool constrain )
   ---------------------------------------- */
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

/* void ElementalMatrix<T>::AlignRowsWith
   ( const DistData& data, bool constrain )
   ---------------------------------------- */
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

/* void ElementalMatrix<T>::AlignAndResize
  ( int colAlign, int rowAlign, Int height, Int width, 
    bool force, bool constrain )
   --------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixAlignAndResize_i
( ElDistMatrix_i A, int colAlign, int rowAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignAndResize_s
( ElDistMatrix_s A, int colAlign, int rowAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignAndResize_d
( ElDistMatrix_d A, int colAlign, int rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignAndResize_c
( ElDistMatrix_c A, int colAlign, int rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignAndResize_z
( ElDistMatrix_z A, int colAlign, int rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );

/* void ElementalMatrix<T>::AlignColsAndResize
   ( int colAlign, Int height, Int width, bool force, bool constrain )
   ------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixAlignColsAndResize_i
( ElDistMatrix_i A, int colAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignColsAndResize_s
( ElDistMatrix_s A, int colAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignColsAndResize_d
( ElDistMatrix_d A, int colAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignColsAndResize_c
( ElDistMatrix_c A, int colAlign, ElInt height, ElInt width, 
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignColsAndResize_z
( ElDistMatrix_z A, int colAlign, ElInt height, ElInt width, 
  bool force, bool constrain );

/* void ElementalMatrix<T>::AlignRowsAndResize
  ( int rowAlign, Int height, Int width, bool force, bool constrain )
   ------------------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixAlignRowsAndResize_i
( ElDistMatrix_i A, int rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRowsAndResize_s
( ElDistMatrix_s A, int rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRowsAndResize_d
( ElDistMatrix_d A, int rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRowsAndResize_c
( ElDistMatrix_c A, int rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );
EL_EXPORT ElError ElDistMatrixAlignRowsAndResize_z
( ElDistMatrix_z A, int rowAlign, ElInt height, ElInt width,
  bool force, bool constrain );

/* void ElementalMatrix<T>::Attach
   ( Int height, Int width, const Grid& grid, int colAlign, int rowAlign, 
     T* buffer, Int ldim, int root )
   ---------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixAttach_i
( ElDistMatrix_i A, ElInt height, ElInt width, ElConstGrid grid,
  int colAlign, int rowAlign, ElInt* buffer, ElInt ldim, 
  int root );
EL_EXPORT ElError ElDistMatrixAttach_s
( ElDistMatrix_s A, ElInt height, ElInt width, ElConstGrid grid,
  int colAlign, int rowAlign, float* buffer, ElInt ldim, 
  int root );
EL_EXPORT ElError ElDistMatrixAttach_d
( ElDistMatrix_d A, ElInt height, ElInt width, ElConstGrid grid,
  int colAlign, int rowAlign, double* buffer, ElInt ldim, 
  int root );
EL_EXPORT ElError ElDistMatrixAttach_c
( ElDistMatrix_c A, ElInt height, ElInt width, ElConstGrid grid,
  int colAlign, int rowAlign, complex_float* buffer, ElInt ldim, 
  int root );
EL_EXPORT ElError ElDistMatrixAttach_z
( ElDistMatrix_z A, ElInt height, ElInt width, ElConstGrid grid,
  int colAlign, int rowAlign, complex_double* buffer, ElInt ldim, 
  int root );

/* void ElementalMatrix<T>::LockedAttach
   ( Int height, Int width, const Grid& grid, int colAlign, int rowAlign, 
     const T* buffer, Int ldim, int root )
   ---------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixLockedAttach_i
( ElDistMatrix_i A, ElInt height, ElInt width, ElConstGrid grid,
  int colAlign, int rowAlign, const ElInt* buffer, 
  ElInt ldim, int root );
EL_EXPORT ElError ElDistMatrixLockedAttach_s
( ElDistMatrix_s A, ElInt height, ElInt width, ElConstGrid grid,
  int colAlign, int rowAlign, const float* buffer, 
  ElInt ldim, int root );
EL_EXPORT ElError ElDistMatrixLockedAttach_d
( ElDistMatrix_d A, ElInt height, ElInt width, ElConstGrid grid,
  int colAlign, int rowAlign, const double* buffer, 
  ElInt ldim, int root );
EL_EXPORT ElError ElDistMatrixLockedAttach_c
( ElDistMatrix_c A, ElInt height, ElInt width, ElConstGrid grid,
  int colAlign, int rowAlign, const complex_float* buffer, 
  ElInt ldim, int root );
EL_EXPORT ElError ElDistMatrixLockedAttach_z
( ElDistMatrix_z A, ElInt height, ElInt width, ElConstGrid grid,
  int colAlign, int rowAlign, const complex_double* buffer, 
  ElInt ldim, int root );

/* Int AbstractDistMatrix<T>::Height() const
   ----------------------------------------- */
EL_EXPORT ElError ElDistMatrixHeight_i( ElConstDistMatrix_i A, ElInt* height );
EL_EXPORT ElError ElDistMatrixHeight_s( ElConstDistMatrix_s A, ElInt* height );
EL_EXPORT ElError ElDistMatrixHeight_d( ElConstDistMatrix_d A, ElInt* height );
EL_EXPORT ElError ElDistMatrixHeight_c( ElConstDistMatrix_c A, ElInt* height );
EL_EXPORT ElError ElDistMatrixHeight_z( ElConstDistMatrix_z A, ElInt* height );

/* Int AbstractDistMatrix<T>::Width() const
   ---------------------------------------- */
EL_EXPORT ElError ElDistMatrixWidth_i( ElConstDistMatrix_i A, ElInt* width );
EL_EXPORT ElError ElDistMatrixWidth_s( ElConstDistMatrix_s A, ElInt* width );
EL_EXPORT ElError ElDistMatrixWidth_d( ElConstDistMatrix_d A, ElInt* width );
EL_EXPORT ElError ElDistMatrixWidth_c( ElConstDistMatrix_c A, ElInt* width );
EL_EXPORT ElError ElDistMatrixWidth_z( ElConstDistMatrix_z A, ElInt* width );

/* Int AbstractDistMatrix<T>::DiagonalLength( Int offset ) const
   ------------------------------------------------------------- */
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

/* bool AbstractDistMatrix<T>::Viewing() const
   -------------------------------------------- */
EL_EXPORT ElError ElDistMatrixViewing_i( ElConstDistMatrix_i A, bool* viewing );
EL_EXPORT ElError ElDistMatrixViewing_s( ElConstDistMatrix_s A, bool* viewing );
EL_EXPORT ElError ElDistMatrixViewing_d( ElConstDistMatrix_d A, bool* viewing );
EL_EXPORT ElError ElDistMatrixViewing_c( ElConstDistMatrix_c A, bool* viewing );
EL_EXPORT ElError ElDistMatrixViewing_z( ElConstDistMatrix_z A, bool* viewing );

/* bool AbstractDistMatrix<T>::Locked() const
   ------------------------------------------ */
EL_EXPORT ElError ElDistMatrixLocked_i( ElConstDistMatrix_i A, bool* locked );
EL_EXPORT ElError ElDistMatrixLocked_s( ElConstDistMatrix_s A, bool* locked );
EL_EXPORT ElError ElDistMatrixLocked_d( ElConstDistMatrix_d A, bool* locked );
EL_EXPORT ElError ElDistMatrixLocked_c( ElConstDistMatrix_c A, bool* locked );
EL_EXPORT ElError ElDistMatrixLocked_z( ElConstDistMatrix_z A, bool* locked );

/* Int AbstractDistMatrix<T>::LocalHeight() const
   ---------------------------------------------- */
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

/* Int AbstractDistMatrix<T>::LocalWidth() const
   --------------------------------------------- */
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

/* Int AbstractDistMatrix<T>::LDim() const
   --------------------------------------- */
EL_EXPORT ElError ElDistMatrixLDim_i( ElConstDistMatrix_i A, ElInt* ldim );
EL_EXPORT ElError ElDistMatrixLDim_s( ElConstDistMatrix_s A, ElInt* ldim );
EL_EXPORT ElError ElDistMatrixLDim_d( ElConstDistMatrix_d A, ElInt* ldim );
EL_EXPORT ElError ElDistMatrixLDim_c( ElConstDistMatrix_c A, ElInt* ldim );
EL_EXPORT ElError ElDistMatrixLDim_z( ElConstDistMatrix_z A, ElInt* ldim );

/* Matrix<T>& AbstractDistMatrix<T>::Matrix() 
   ------------------------------------------ */
EL_EXPORT ElError ElDistMatrixMatrix_i( ElDistMatrix_i A, ElMatrix_i* ALoc );
EL_EXPORT ElError ElDistMatrixMatrix_s( ElDistMatrix_s A, ElMatrix_s* ALoc );
EL_EXPORT ElError ElDistMatrixMatrix_d( ElDistMatrix_d A, ElMatrix_d* ALoc );
EL_EXPORT ElError ElDistMatrixMatrix_c( ElDistMatrix_c A, ElMatrix_c* ALoc );
EL_EXPORT ElError ElDistMatrixMatrix_z( ElDistMatrix_z A, ElMatrix_z* ALoc );

/* const Matrix<T>& AbstractDistMatrix<T>::LockedMatrix() const
   ------------------------------------------------------------ */
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

/* size_t AbstractDistMatrix<T>::AllocatedMemory() const
   ----------------------------------------------------- */
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

/* T* AbstractDistMatrix<T>::Buffer()
   ---------------------------------- */
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

/* const T* AbstractDistMatrix<T>::LockedBuffer() const
   ---------------------------------------------------- */
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

/* const Grid& AbstractDistMatrix<T>::Grid() const
   ----------------------------------------------- */
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

/* bool AbstractDistMatrix<T>::ColConstrained() const
   -------------------------------------------------- */
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

/* bool AbstractDistMatrix<T>::RowConstrained() const
   -------------------------------------------------- */
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

/* bool AbstractDistMatrix<T>::RootConstrained() const
   --------------------------------------------------- */
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

/* int AbstractDistMatrix<T>::ColAlign() const
   ---------------------------------------- */
EL_EXPORT ElError ElDistMatrixColAlign_i
( ElConstDistMatrix_i A, int* colAlign );
EL_EXPORT ElError ElDistMatrixColAlign_s
( ElConstDistMatrix_s A, int* colAlign );
EL_EXPORT ElError ElDistMatrixColAlign_d
( ElConstDistMatrix_d A, int* colAlign );
EL_EXPORT ElError ElDistMatrixColAlign_c
( ElConstDistMatrix_c A, int* colAlign );
EL_EXPORT ElError ElDistMatrixColAlign_z
( ElConstDistMatrix_z A, int* colAlign );

/* int AbstractDistMatrix<T>::RowAlign() const
   ------------------------------------------- */
EL_EXPORT ElError ElDistMatrixRowAlign_i
( ElConstDistMatrix_i A, int* rowAlign );
EL_EXPORT ElError ElDistMatrixRowAlign_s
( ElConstDistMatrix_s A, int* rowAlign );
EL_EXPORT ElError ElDistMatrixRowAlign_d
( ElConstDistMatrix_d A, int* rowAlign );
EL_EXPORT ElError ElDistMatrixRowAlign_c
( ElConstDistMatrix_c A, int* rowAlign );
EL_EXPORT ElError ElDistMatrixRowAlign_z
( ElConstDistMatrix_z A, int* rowAlign );

/* int AbstractDistMatrix<T>::ColShift() const
   ------------------------------------------- */
EL_EXPORT ElError ElDistMatrixColShift_i
( ElConstDistMatrix_i A, int* colShift );
EL_EXPORT ElError ElDistMatrixColShift_s
( ElConstDistMatrix_s A, int* colShift );
EL_EXPORT ElError ElDistMatrixColShift_d
( ElConstDistMatrix_d A, int* colShift );
EL_EXPORT ElError ElDistMatrixColShift_c
( ElConstDistMatrix_c A, int* colShift );
EL_EXPORT ElError ElDistMatrixColShift_z
( ElConstDistMatrix_z A, int* colShift );

/* int AbstractDistMatrix<T>::RowShift() const
   ------------------------------------------- */
EL_EXPORT ElError ElDistMatrixRowShift_i
( ElConstDistMatrix_i A, int* rowShift );
EL_EXPORT ElError ElDistMatrixRowShift_s
( ElConstDistMatrix_s A, int* rowShift );
EL_EXPORT ElError ElDistMatrixRowShift_d
( ElConstDistMatrix_d A, int* rowShift );
EL_EXPORT ElError ElDistMatrixRowShift_c
( ElConstDistMatrix_c A, int* rowShift );
EL_EXPORT ElError ElDistMatrixRowShift_z
( ElConstDistMatrix_z A, int* rowShift );

/* int AbstractDistMatrix<T>::ColRank() const
   ------------------------------------------ */
EL_EXPORT ElError ElDistMatrixColRank_i( ElConstDistMatrix_i A, int* rank );
EL_EXPORT ElError ElDistMatrixColRank_s( ElConstDistMatrix_s A, int* rank );
EL_EXPORT ElError ElDistMatrixColRank_d( ElConstDistMatrix_d A, int* rank );
EL_EXPORT ElError ElDistMatrixColRank_c( ElConstDistMatrix_c A, int* rank );
EL_EXPORT ElError ElDistMatrixColRank_z( ElConstDistMatrix_z A, int* rank );

/* int AbstractDistMatrix<T>::RowRank() const
   ------------------------------------------ */
EL_EXPORT ElError ElDistMatrixRowRank_i( ElConstDistMatrix_i A, int* rank );
EL_EXPORT ElError ElDistMatrixRowRank_s( ElConstDistMatrix_s A, int* rank );
EL_EXPORT ElError ElDistMatrixRowRank_d( ElConstDistMatrix_d A, int* rank );
EL_EXPORT ElError ElDistMatrixRowRank_c( ElConstDistMatrix_c A, int* rank );
EL_EXPORT ElError ElDistMatrixRowRank_z( ElConstDistMatrix_z A, int* rank );

/* int AbstractDistMatrix<T>::PartialColRank() const
   ------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialColRank_i
( ElConstDistMatrix_i A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialColRank_s
( ElConstDistMatrix_s A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialColRank_d
( ElConstDistMatrix_d A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialColRank_c
( ElConstDistMatrix_c A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialColRank_z
( ElConstDistMatrix_z A, int* rank );

/* int AbstractDistMatrix<T>::PartialRowRank() const
   ------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialRowRank_i
( ElConstDistMatrix_i A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialRowRank_s
( ElConstDistMatrix_s A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialRowRank_d
( ElConstDistMatrix_d A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialRowRank_c
( ElConstDistMatrix_c A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialRowRank_z
( ElConstDistMatrix_z A, int* rank );

/* int AbstractDistMatrix<T>::PartialUnionColRank() const
   ------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixPartialUnionColRank_i
( ElConstDistMatrix_i A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionColRank_s
( ElConstDistMatrix_s A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionColRank_d
( ElConstDistMatrix_d A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionColRank_c
( ElConstDistMatrix_c A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionColRank_z
( ElConstDistMatrix_z A, int* rank );

/* int AbstractDistMatrix<T>::PartialUnionRowRank() const
   ------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixPartialUnionRowRank_i
( ElConstDistMatrix_i A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionRowRank_s
( ElConstDistMatrix_s A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionRowRank_d
( ElConstDistMatrix_d A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionRowRank_c
( ElConstDistMatrix_c A, int* rank );
EL_EXPORT ElError ElDistMatrixPartialUnionRowRank_z
( ElConstDistMatrix_z A, int* rank );

/* int AbstractDistMatrix<T>::DistRank() const
   ------------------------------------------- */
EL_EXPORT ElError ElDistMatrixDistRank_i
( ElConstDistMatrix_i A, int* rank );
EL_EXPORT ElError ElDistMatrixDistRank_s
( ElConstDistMatrix_s A, int* rank );
EL_EXPORT ElError ElDistMatrixDistRank_d
( ElConstDistMatrix_d A, int* rank );
EL_EXPORT ElError ElDistMatrixDistRank_c
( ElConstDistMatrix_c A, int* rank );
EL_EXPORT ElError ElDistMatrixDistRank_z
( ElConstDistMatrix_z A, int* rank );

/* int AbstractDistMatrix<T>::CrossRank() const
   -------------------------------------------- */
EL_EXPORT ElError ElDistMatrixCrossRank_i( ElConstDistMatrix_i A, int* rank );
EL_EXPORT ElError ElDistMatrixCrossRank_s( ElConstDistMatrix_s A, int* rank );
EL_EXPORT ElError ElDistMatrixCrossRank_d( ElConstDistMatrix_d A, int* rank );
EL_EXPORT ElError ElDistMatrixCrossRank_c( ElConstDistMatrix_c A, int* rank );
EL_EXPORT ElError ElDistMatrixCrossRank_z( ElConstDistMatrix_z A, int* rank );

/* int AbstractDistMatrix<T>::RedundantRank() const
   ------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixRedundantRank_i
( ElConstDistMatrix_i A, int* rank );
EL_EXPORT ElError ElDistMatrixRedundantRank_s
( ElConstDistMatrix_s A, int* rank );
EL_EXPORT ElError ElDistMatrixRedundantRank_d
( ElConstDistMatrix_d A, int* rank );
EL_EXPORT ElError ElDistMatrixRedundantRank_c
( ElConstDistMatrix_c A, int* rank );
EL_EXPORT ElError ElDistMatrixRedundantRank_z
( ElConstDistMatrix_z A, int* rank );

/* int AbstractDistMatrix<T>::Root() const
   --------------------------------------- */
EL_EXPORT ElError ElDistMatrixRoot_i( ElConstDistMatrix_i A, int* root );
EL_EXPORT ElError ElDistMatrixRoot_s( ElConstDistMatrix_s A, int* root );
EL_EXPORT ElError ElDistMatrixRoot_d( ElConstDistMatrix_d A, int* root );
EL_EXPORT ElError ElDistMatrixRoot_c( ElConstDistMatrix_c A, int* root );
EL_EXPORT ElError ElDistMatrixRoot_z( ElConstDistMatrix_z A, int* root );

/* bool AbstractDistMatrix<T>::Participating() const
   ------------------------------------------------- */
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

/* int AbstractDistMatrix<T>::RowOwner( Int i ) const
   -------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixRowOwner_i
( ElConstDistMatrix_i A, ElInt i, int* rowOwner );
EL_EXPORT ElError ElDistMatrixRowOwner_s
( ElConstDistMatrix_s A, ElInt i, int* rowOwner );
EL_EXPORT ElError ElDistMatrixRowOwner_d
( ElConstDistMatrix_d A, ElInt i, int* rowOwner );
EL_EXPORT ElError ElDistMatrixRowOwner_c
( ElConstDistMatrix_c A, ElInt i, int* rowOwner );
EL_EXPORT ElError ElDistMatrixRowOwner_z
( ElConstDistMatrix_z A, ElInt i, int* rowOwner );

/* int AbstractDistMatrix<T>::ColOwner( Int j ) const
   -------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixColOwner_i
( ElConstDistMatrix_i A, ElInt j, int* colOwner );
EL_EXPORT ElError ElDistMatrixColOwner_s
( ElConstDistMatrix_s A, ElInt j, int* colOwner );
EL_EXPORT ElError ElDistMatrixColOwner_d
( ElConstDistMatrix_d A, ElInt j, int* colOwner );
EL_EXPORT ElError ElDistMatrixColOwner_c
( ElConstDistMatrix_c A, ElInt j, int* colOwner );
EL_EXPORT ElError ElDistMatrixColOwner_z
( ElConstDistMatrix_z A, ElInt j, int* colOwner );

/* int AbstractDistMatrix<T>::Owner( Int i, Int j ) const
   ------------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixOwner_i
( ElConstDistMatrix_i A, ElInt i, ElInt j, int* owner );
EL_EXPORT ElError ElDistMatrixOwner_s
( ElConstDistMatrix_s A, ElInt i, ElInt j, int* owner );
EL_EXPORT ElError ElDistMatrixOwner_d
( ElConstDistMatrix_d A, ElInt i, ElInt j, int* owner );
EL_EXPORT ElError ElDistMatrixOwner_c
( ElConstDistMatrix_c A, ElInt i, ElInt j, int* owner );
EL_EXPORT ElError ElDistMatrixOwner_z
( ElConstDistMatrix_z A, ElInt i, ElInt j, int* owner );

/* Int AbstractDistMatrix<T>::LocalRow( Int i ) const
   -------------------------------------------------- */
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

/* Int AbstractDistMatrix<T>::LocalCol( Int j ) const
   -------------------------------------------------- */
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

/* Int AbstractDistMatrix<T>::LocalRowOffset( Int i ) const
   -------------------------------------------------------- */
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

/* Int AbstractDistMatrix<T>::LocalColOffset( Int j ) const
   -------------------------------------------------------- */
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

/* Int AbstractDistMatrix<T>::GlobalRow( Int iLoc ) const
   ------------------------------------------------------ */
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

/* Int AbstractDistMatrix<T>::GlobalCol( Int jLoc ) const
   ------------------------------------------------------ */
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

/* bool AbstractDistMatrix<T>::IsLocalRow( Int i ) const
   ------------------------------------------------------ */
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

/* bool AbstractDistMatrix<T>::IsLocalCol( Int j ) const
   ----------------------------------------------------- */
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

/* bool AbstractDistMatrix<T>::IsLocal( Int i, Int j ) const
   --------------------------------------------------------- */
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

/* DistData ElementalMatrix<T>::DistData() const
   --------------------------------------------- */
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

/* mpi::Comm AbstractDistMatrix<T>::DistComm() const
   ------------------------------------------------- */
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

/* mpi::Comm AbstractDistMatrix<T>::CrossComm() const
   -------------------------------------------------- */
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

/* mpi::Comm AbstractDistMatrix<T>::RedundantComm() const
   ------------------------------------------------------ */
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

/* mpi::Comm AbstractDistMatrix<T>::ColComm() const
   ------------------------------------------------ */
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

/* mpi::Comm AbstractDistMatrix<T>::RowComm() const
   ------------------------------------------------ */
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

/* mpi::Comm AbstractDistMatrix<T>::PartialColComm() const
   ------------------------------------------------------- */
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

/* mpi::Comm AbstractDistMatrix<T>::PartialRowComm() const
   ------------------------------------------------------- */
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

/* mpi::Comm AbstractDistMatrix<T>::PartialUnionColComm() const
   ------------------------------------------------------------ */
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

/* mpi::Comm AbstractDistMatrix<T>::PartialUnionRowComm() const
   ------------------------------------------------------------ */
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

/* int AbstractDistMatrix<T>::ColStride() const
   -------------------------------------------- */
EL_EXPORT ElError ElDistMatrixColStride_i
( ElConstDistMatrix_i A, int* stride );
EL_EXPORT ElError ElDistMatrixColStride_s
( ElConstDistMatrix_s A, int* stride );
EL_EXPORT ElError ElDistMatrixColStride_d
( ElConstDistMatrix_d A, int* stride );
EL_EXPORT ElError ElDistMatrixColStride_c
( ElConstDistMatrix_c A, int* stride );
EL_EXPORT ElError ElDistMatrixColStride_z
( ElConstDistMatrix_z A, int* stride );

/* int AbstractDistMatrix<T>::RowStride() const
   -------------------------------------------- */
EL_EXPORT ElError ElDistMatrixRowStride_i
( ElConstDistMatrix_i A, int* stride );
EL_EXPORT ElError ElDistMatrixRowStride_s
( ElConstDistMatrix_s A, int* stride );
EL_EXPORT ElError ElDistMatrixRowStride_d
( ElConstDistMatrix_d A, int* stride );
EL_EXPORT ElError ElDistMatrixRowStride_c
( ElConstDistMatrix_c A, int* stride );
EL_EXPORT ElError ElDistMatrixRowStride_z
( ElConstDistMatrix_z A, int* stride );

/* int AbstractDistMatrix<T>::PartialColStride() const
   --------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialColStride_i
( ElConstDistMatrix_i A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialColStride_s
( ElConstDistMatrix_s A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialColStride_d
( ElConstDistMatrix_d A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialColStride_c
( ElConstDistMatrix_c A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialColStride_z
( ElConstDistMatrix_z A, int* stride );

/* int AbstractDistMatrix<T>::PartialRowStride() const
   --------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialRowStride_i
( ElConstDistMatrix_i A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialRowStride_s
( ElConstDistMatrix_s A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialRowStride_d
( ElConstDistMatrix_d A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialRowStride_c
( ElConstDistMatrix_c A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialRowStride_z
( ElConstDistMatrix_z A, int* stride );

/* int AbstractDistMatrix<T>::PartialUnionColStride() const
   -------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialUnionColStride_i
( ElConstDistMatrix_i A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionColStride_s
( ElConstDistMatrix_s A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionColStride_d
( ElConstDistMatrix_d A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionColStride_c
( ElConstDistMatrix_c A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionColStride_z
( ElConstDistMatrix_z A, int* stride );

/* int AbstractDistMatrix<T>::PartialUnionRowStride() const
   -------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixPartialUnionRowStride_i
( ElConstDistMatrix_i A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionRowStride_s
( ElConstDistMatrix_s A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionRowStride_d
( ElConstDistMatrix_d A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionRowStride_c
( ElConstDistMatrix_c A, int* stride );
EL_EXPORT ElError ElDistMatrixPartialUnionRowStride_z
( ElConstDistMatrix_z A, int* stride );

/* int AbstractDistMatrix<T>::DistSize() const
   ------------------------------------------- */
EL_EXPORT ElError ElDistMatrixDistSize_i
( ElConstDistMatrix_i A, int* commSize );
EL_EXPORT ElError ElDistMatrixDistSize_s
( ElConstDistMatrix_s A, int* commSize );
EL_EXPORT ElError ElDistMatrixDistSize_d
( ElConstDistMatrix_d A, int* commSize );
EL_EXPORT ElError ElDistMatrixDistSize_c
( ElConstDistMatrix_c A, int* commSize );
EL_EXPORT ElError ElDistMatrixDistSize_z
( ElConstDistMatrix_z A, int* commSize );

/* int AbstractDistMatrix<T>::CrossSize() const
   -------------------------------------------- */
EL_EXPORT ElError ElDistMatrixCrossSize_i
( ElConstDistMatrix_i A, int* commSize );
EL_EXPORT ElError ElDistMatrixCrossSize_s
( ElConstDistMatrix_s A, int* commSize );
EL_EXPORT ElError ElDistMatrixCrossSize_d
( ElConstDistMatrix_d A, int* commSize );
EL_EXPORT ElError ElDistMatrixCrossSize_c
( ElConstDistMatrix_c A, int* commSize );
EL_EXPORT ElError ElDistMatrixCrossSize_z
( ElConstDistMatrix_z A, int* commSize );

/* int AbstractDistMatrix<T>::RedundantSize() const
   ------------------------------------------------ */
EL_EXPORT ElError ElDistMatrixRedundantSize_i
( ElConstDistMatrix_i A, int* commSize );
EL_EXPORT ElError ElDistMatrixRedundantSize_s
( ElConstDistMatrix_s A, int* commSize );
EL_EXPORT ElError ElDistMatrixRedundantSize_d
( ElConstDistMatrix_d A, int* commSize );
EL_EXPORT ElError ElDistMatrixRedundantSize_c
( ElConstDistMatrix_c A, int* commSize );
EL_EXPORT ElError ElDistMatrixRedundantSize_z
( ElConstDistMatrix_z A, int* commSize );

/* T AbstractDistMatrix<T>::Get( Int i, Int j ) const
   -------------------------------------------------- */
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

/* Base<T> AbstractDistMatrix<T>::GetRealPart( Int i, Int j ) const
   ---------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGetRealPart_c
( ElConstDistMatrix_c A, ElInt i, ElInt j, float* val );
EL_EXPORT ElError ElDistMatrixGetRealPart_z
( ElConstDistMatrix_z A, ElInt i, ElInt j, double* val );

/* Base<T> AbstractDistMatrix<T>::GetImagPart( Int i, Int j ) const
   ---------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixGetImagPart_c
( ElConstDistMatrix_c A, ElInt i, ElInt j, float* val );
EL_EXPORT ElError ElDistMatrixGetImagPart_z
( ElConstDistMatrix_z A, ElInt i, ElInt j, double* val );

/* void AbstractDistMatrix<T>::Set( Int i, Int j, T alpha )
   -------------------------------------------------------- */
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

/* void AbstractDistMatrix<T>::SetRealPart( Int i, Int j, Base<T> alpha )
   ---------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSetRealPart_c
( ElDistMatrix_c A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElDistMatrixSetRealPart_z
( ElDistMatrix_z A, ElInt i, ElInt j, double alpha );

/* void AbstractDistMatrix<T>::SetImagPart( Int i, Int j, Base<T> alpha )
   ---------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSetImagPart_c
( ElDistMatrix_c A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElDistMatrixSetImagPart_z
( ElDistMatrix_z A, ElInt i, ElInt j, double alpha );

/* void AbstractDistMatrix<T>::Update( Int i, Int j, T alpha )
   ----------------------------------------------------------- */
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

/* void AbstractDistMatrix<T>::UpdateRealPart( Int i, Int j, Base<T> alpha )
   ------------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixUpdateRealPart_c
( ElDistMatrix_c A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElDistMatrixUpdateRealPart_z
( ElDistMatrix_z A, ElInt i, ElInt j, double alpha );

/* void AbstractDistMatrix<T>::UpdateImagPart( Int i, Int j, Base<T> alpha )
   ------------------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixUpdateImagPart_c
( ElDistMatrix_c A, ElInt i, ElInt j, float alpha );
EL_EXPORT ElError ElDistMatrixUpdateImagPart_z
( ElDistMatrix_z A, ElInt i, ElInt j, double alpha );

/* void AbstractDistMatrix<T>::MakeReal( Int i, Int j )
   ---------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixMakeReal_c( ElDistMatrix_c A, ElInt i, ElInt j );
EL_EXPORT ElError ElDistMatrixMakeReal_z( ElDistMatrix_z A, ElInt i, ElInt j );

/* void AbstractDistMatrix<T>::Conjugate( Int i, Int j )
   ----------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixConjugate_c( ElDistMatrix_c A, ElInt i, ElInt j );
EL_EXPORT ElError ElDistMatrixConjugate_z( ElDistMatrix_z A, ElInt i, ElInt j );

/* void AbstractDistMatrix<T>::Reserve( Int numRemoteEntries )
   ----------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixReserve_i
( ElDistMatrix_i A, ElInt numRemoteEntries );
EL_EXPORT ElError ElDistMatrixReserve_s
( ElDistMatrix_s A, ElInt numRemoteEntries );
EL_EXPORT ElError ElDistMatrixReserve_d
( ElDistMatrix_d A, ElInt numRemoteEntries );
EL_EXPORT ElError ElDistMatrixReserve_c
( ElDistMatrix_c A, ElInt numRemoteEntries );
EL_EXPORT ElError ElDistMatrixReserve_z
( ElDistMatrix_z A, ElInt numRemoteEntries );

/* void AbstractDistMatrix<T>::QueueUpdate( Int i, Int j, T value )
   ---------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixQueueUpdate_i
( ElDistMatrix_i A, ElInt i, ElInt j, ElInt value );
EL_EXPORT ElError ElDistMatrixQueueUpdate_s
( ElDistMatrix_s A, ElInt i, ElInt j, float value );
EL_EXPORT ElError ElDistMatrixQueueUpdate_d
( ElDistMatrix_d A, ElInt i, ElInt j, double value );
EL_EXPORT ElError ElDistMatrixQueueUpdate_c
( ElDistMatrix_c A, ElInt i, ElInt j, complex_float value );
EL_EXPORT ElError ElDistMatrixQueueUpdate_z
( ElDistMatrix_z A, ElInt i, ElInt j, complex_double value );

/* void AbstractDistMatrix<T>::ProcessQueues()
   ------------------------------------------- */
EL_EXPORT ElError ElDistMatrixProcessQueues_i( ElDistMatrix_i A );
EL_EXPORT ElError ElDistMatrixProcessQueues_s( ElDistMatrix_s A );
EL_EXPORT ElError ElDistMatrixProcessQueues_d( ElDistMatrix_d A );
EL_EXPORT ElError ElDistMatrixProcessQueues_c( ElDistMatrix_c A );
EL_EXPORT ElError ElDistMatrixProcessQueues_z( ElDistMatrix_z A );

/* void AbstractDistMatrix<T>::ReservePulls( Int numPulls ) const
   -------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixReservePulls_i
( ElConstDistMatrix_i A, ElInt numPulls );
EL_EXPORT ElError ElDistMatrixReservePulls_s
( ElConstDistMatrix_s A, ElInt numPulls );
EL_EXPORT ElError ElDistMatrixReservePulls_d
( ElConstDistMatrix_d A, ElInt numPulls );
EL_EXPORT ElError ElDistMatrixReservePulls_c
( ElConstDistMatrix_c A, ElInt numPulls );
EL_EXPORT ElError ElDistMatrixReservePulls_z
( ElConstDistMatrix_z A, ElInt numPulls );

/* void AbstractDistMatrix<T>::QueuePull( Int i, Int j ) const
   ----------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixQueuePull_i
( ElConstDistMatrix_i A, ElInt i, ElInt j );
EL_EXPORT ElError ElDistMatrixQueuePull_s
( ElConstDistMatrix_s A, ElInt i, ElInt j );
EL_EXPORT ElError ElDistMatrixQueuePull_d
( ElConstDistMatrix_d A, ElInt i, ElInt j );
EL_EXPORT ElError ElDistMatrixQueuePull_c
( ElConstDistMatrix_c A, ElInt i, ElInt j );
EL_EXPORT ElError ElDistMatrixQueuePull_z
( ElConstDistMatrix_z A, ElInt i, ElInt j );

/* void AbstractDistMatrix<T>::ProcessPullQueue( T* pullBuf ) const
   ---------------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixProcessPullQueue_i
( ElConstDistMatrix_i A, ElInt* pullBuf );
EL_EXPORT ElError ElDistMatrixProcessPullQueue_s
( ElConstDistMatrix_s A, float* pullBuf );
EL_EXPORT ElError ElDistMatrixProcessPullQueue_d
( ElConstDistMatrix_d A, double* pullBuf );
EL_EXPORT ElError ElDistMatrixProcessPullQueue_c
( ElConstDistMatrix_c A, complex_float* pullBuf );
EL_EXPORT ElError ElDistMatrixProcessPullQueue_z
( ElConstDistMatrix_z A, complex_double* pullBuf );

/* T AbstractDistMatrix<T>::GetLocal( Int iLoc, Int jLoc ) const
   ------------------------------------------------------------- */
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

/* Base<T> AbstractDistMatrix<T>::GetLocalRealPart( Int iLoc, Int jLoc ) const
   ---------------------------------------------------------------------------*/
EL_EXPORT ElError ElDistMatrixGetLocalRealPart_c
( ElConstDistMatrix_c A, ElInt iLoc, ElInt jLoc, float* val );
EL_EXPORT ElError ElDistMatrixGetLocalRealPart_z
( ElConstDistMatrix_z A, ElInt iLoc, ElInt jLoc, double* val );

/* Base<T> AbstractDistMatrix<T>::GetLocalImagPart( Int iLoc, Int jLoc ) const
   ---------------------------------------------------------------------------*/
EL_EXPORT ElError ElDistMatrixGetLocalImagPart_c
( ElConstDistMatrix_c A, ElInt iLoc, ElInt jLoc, float* val );
EL_EXPORT ElError ElDistMatrixGetLocalImagPart_z
( ElConstDistMatrix_z A, ElInt iLoc, ElInt jLoc, double* val );

/* void AbstractDistMatrix<T>::SetLocal( Int iLoc, Int jLoc, T val )
   ----------------------------------------------------------------- */
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

/* void AbstractDistMatrix<T>::SetLocalRealPart
   ( Int iLoc, Int jLoc, Base<T> val )
   -------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSetLocalRealPart_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, float val );
EL_EXPORT ElError ElDistMatrixSetLocalRealPart_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, double val );

/* void AbstractDistMatrix<T>::SetLocalRealPart
   ( Int iLoc, Int jLoc, Base<T> val )
   -------------------------------------------- */
EL_EXPORT ElError ElDistMatrixSetLocalImagPart_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, float val );
EL_EXPORT ElError ElDistMatrixSetLocalImagPart_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, double val );

/* void AbstractDistMatrix<T>::UpdateLocal( Int iLoc, Int jLoc, T val )
   -------------------------------------------------------------------- */
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

/* void AbstractDistMatrix<T>::UpdateLocalRealPart
   ( Int iLoc, Int jLoc, Base<T> val )
   ----------------------------------------------- */
EL_EXPORT ElError ElDistMatrixUpdateLocalRealPart_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, float val );
EL_EXPORT ElError ElDistMatrixUpdateLocalRealPart_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, double val );

/* void AbstractDistMatrix<T>::UpdateLocalImagPart
   ( Int iLoc, Int jLoc, Base<T> val )
   ----------------------------------------------- */
EL_EXPORT ElError ElDistMatrixUpdateLocalImagPart_c
( ElDistMatrix_c A, ElInt iLoc, ElInt jLoc, float val );
EL_EXPORT ElError ElDistMatrixUpdateLocalImagPart_z
( ElDistMatrix_z A, ElInt iLoc, ElInt jLoc, double val );

/* bool ElementalMatrix<T>::DiagonalAlignedWith
   ( const DistData& data, Int offset ) const
   -------------------------------------------- */
EL_EXPORT ElError ElDistMatrixDiagonalAlignedWith_i
( ElConstDistMatrix_i A, ElDistData distData, ElInt offset, bool* aligned );
EL_EXPORT ElError ElDistMatrixDiagonalAlignedWith_s
( ElConstDistMatrix_s A, ElDistData distData, ElInt offset, bool* aligned );
EL_EXPORT ElError ElDistMatrixDiagonalAlignedWith_d
( ElConstDistMatrix_d A, ElDistData distData, ElInt offset, bool* aligned );
EL_EXPORT ElError ElDistMatrixDiagonalAlignedWith_c
( ElConstDistMatrix_c A, ElDistData distData, ElInt offset, bool* aligned );
EL_EXPORT ElError ElDistMatrixDiagonalAlignedWith_z
( ElConstDistMatrix_z A, ElDistData distData, ElInt offset, bool* aligned );

/* int ElementalMatrix<T>::DiagonalRoot( Int offset ) const
   -------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixDiagonalRoot_i
( ElConstDistMatrix_i A, ElInt offset, int* root );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_s
( ElConstDistMatrix_s A, ElInt offset, int* root );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_d
( ElConstDistMatrix_d A, ElInt offset, int* root );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_c
( ElConstDistMatrix_c A, ElInt offset, int* root );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_z
( ElConstDistMatrix_z A, ElInt offset, int* root );

/* int ElementalMatrix<T>::DiagonalAlign( Int offset ) const
   --------------------------------------------------------- */
EL_EXPORT ElError ElDistMatrixDiagonalRoot_i
( ElConstDistMatrix_i A, ElInt offset, int* align );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_s
( ElConstDistMatrix_s A, ElInt offset, int* align );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_d
( ElConstDistMatrix_d A, ElInt offset, int* align );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_c
( ElConstDistMatrix_c A, ElInt offset, int* align );
EL_EXPORT ElError ElDistMatrixDiagonalRoot_z
( ElConstDistMatrix_z A, ElInt offset, int* align );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_DISTMATRIX_C_H */
