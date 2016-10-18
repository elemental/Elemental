/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SPARSEMATRIX_C_H
#define EL_SPARSEMATRIX_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* An anonymous struct meant as a placeholder for SparseMatrix<T>
   -------------------------------------------------------------- */
typedef struct ElSparseMatrix_iDummy* ElSparseMatrix_i;
typedef struct ElSparseMatrix_sDummy* ElSparseMatrix_s;
typedef struct ElSparseMatrix_dDummy* ElSparseMatrix_d;
typedef struct ElSparseMatrix_cDummy* ElSparseMatrix_c;
typedef struct ElSparseMatrix_zDummy* ElSparseMatrix_z;

typedef const struct ElSparseMatrix_iDummy* ElConstSparseMatrix_i;
typedef const struct ElSparseMatrix_sDummy* ElConstSparseMatrix_s;
typedef const struct ElSparseMatrix_dDummy* ElConstSparseMatrix_d;
typedef const struct ElSparseMatrix_cDummy* ElConstSparseMatrix_c;
typedef const struct ElSparseMatrix_zDummy* ElConstSparseMatrix_z;

/* Constructors and destructors
   ============================ */

/* SparseMatrix<T>::SparseMatrix()
   ------------------------------- */
EL_EXPORT ElError ElSparseMatrixCreate_i( ElSparseMatrix_i* A );
EL_EXPORT ElError ElSparseMatrixCreate_s( ElSparseMatrix_s* A );
EL_EXPORT ElError ElSparseMatrixCreate_d( ElSparseMatrix_d* A );
EL_EXPORT ElError ElSparseMatrixCreate_c( ElSparseMatrix_c* A );
EL_EXPORT ElError ElSparseMatrixCreate_z( ElSparseMatrix_z* A );

/* SparseMatrix<T>::~SparseMatrix()
   -------------------------------- */
EL_EXPORT ElError ElSparseMatrixDestroy_i( ElConstSparseMatrix_i A );
EL_EXPORT ElError ElSparseMatrixDestroy_s( ElConstSparseMatrix_s A );
EL_EXPORT ElError ElSparseMatrixDestroy_d( ElConstSparseMatrix_d A );
EL_EXPORT ElError ElSparseMatrixDestroy_c( ElConstSparseMatrix_c A );
EL_EXPORT ElError ElSparseMatrixDestroy_z( ElConstSparseMatrix_z A );

/* Assignment and reconfiguration
   ============================== */

/* void SparseMatrix<T>::Empty()
   ----------------------------- */
EL_EXPORT ElError ElSparseMatrixEmpty_i( ElSparseMatrix_i A );
EL_EXPORT ElError ElSparseMatrixEmpty_s( ElSparseMatrix_s A );
EL_EXPORT ElError ElSparseMatrixEmpty_d( ElSparseMatrix_d A );
EL_EXPORT ElError ElSparseMatrixEmpty_c( ElSparseMatrix_c A );
EL_EXPORT ElError ElSparseMatrixEmpty_z( ElSparseMatrix_z A );

/* void SparseMatrix<T>::Resize( Int height, Int width )
   ----------------------------------------------------- */
EL_EXPORT ElError ElSparseMatrixResize_i
( ElSparseMatrix_i A, ElInt height, ElInt width );
EL_EXPORT ElError ElSparseMatrixResize_s
( ElSparseMatrix_s A, ElInt height, ElInt width );
EL_EXPORT ElError ElSparseMatrixResize_d
( ElSparseMatrix_d A, ElInt height, ElInt width );
EL_EXPORT ElError ElSparseMatrixResize_c
( ElSparseMatrix_c A, ElInt height, ElInt width );
EL_EXPORT ElError ElSparseMatrixResize_z
( ElSparseMatrix_z A, ElInt height, ElInt width );

/* void SparseMatrix<T>::Reserve( Int numEntries )
   ----------------------------------------------- */
EL_EXPORT ElError ElSparseMatrixReserve_i
( ElSparseMatrix_i A, ElInt numEntries );
EL_EXPORT ElError ElSparseMatrixReserve_s
( ElSparseMatrix_s A, ElInt numEntries );
EL_EXPORT ElError ElSparseMatrixReserve_d
( ElSparseMatrix_d A, ElInt numEntries );
EL_EXPORT ElError ElSparseMatrixReserve_c
( ElSparseMatrix_c A, ElInt numEntries );
EL_EXPORT ElError ElSparseMatrixReserve_z
( ElSparseMatrix_z A, ElInt numEntries );

/* void SparseMatrix<T>::Update( Int row, Int col, T value )
   --------------------------------------------------------- */
EL_EXPORT ElError ElSparseMatrixUpdate_i
( ElSparseMatrix_i A, ElInt row, ElInt col, ElInt value );
EL_EXPORT ElError ElSparseMatrixUpdate_s
( ElSparseMatrix_s A, ElInt row, ElInt col, float value );
EL_EXPORT ElError ElSparseMatrixUpdate_d
( ElSparseMatrix_d A, ElInt row, ElInt col, double value );
EL_EXPORT ElError ElSparseMatrixUpdate_c
( ElSparseMatrix_c A, ElInt row, ElInt col, complex_float value );
EL_EXPORT ElError ElSparseMatrixUpdate_z
( ElSparseMatrix_z A, ElInt row, ElInt col, complex_double value );

/* void SparseMatrix<T>::Zero( Int row, Int col )
   ---------------------------------------------- */
EL_EXPORT ElError ElSparseMatrixZero_i
( ElSparseMatrix_i A, ElInt row, ElInt col );
EL_EXPORT ElError ElSparseMatrixZero_s
( ElSparseMatrix_s A, ElInt row, ElInt col );
EL_EXPORT ElError ElSparseMatrixZero_d
( ElSparseMatrix_d A, ElInt row, ElInt col );
EL_EXPORT ElError ElSparseMatrixZero_c
( ElSparseMatrix_c A, ElInt row, ElInt col );
EL_EXPORT ElError ElSparseMatrixZero_z
( ElSparseMatrix_z A, ElInt row, ElInt col );

/* void SparseMatrix<T>::QueueUpdate( Int row, Int col, T value )
   -------------------------------------------------------------- */
EL_EXPORT ElError ElSparseMatrixQueueUpdate_i
( ElSparseMatrix_i A, ElInt row, ElInt col, ElInt value );
EL_EXPORT ElError ElSparseMatrixQueueUpdate_s
( ElSparseMatrix_s A, ElInt row, ElInt col, float value );
EL_EXPORT ElError ElSparseMatrixQueueUpdate_d
( ElSparseMatrix_d A, ElInt row, ElInt col, double value );
EL_EXPORT ElError ElSparseMatrixQueueUpdate_c
( ElSparseMatrix_c A, ElInt row, ElInt col, complex_float value );
EL_EXPORT ElError ElSparseMatrixQueueUpdate_z
( ElSparseMatrix_z A, ElInt row, ElInt col, complex_double value );

/* void SparseMatrix<T>::QueueZero( Int row, Int col )
   --------------------------------------------------- */
EL_EXPORT ElError ElSparseMatrixQueueZero_i
( ElSparseMatrix_i A, ElInt row, ElInt col );
EL_EXPORT ElError ElSparseMatrixQueueZero_s
( ElSparseMatrix_s A, ElInt row, ElInt col );
EL_EXPORT ElError ElSparseMatrixQueueZero_d
( ElSparseMatrix_d A, ElInt row, ElInt col );
EL_EXPORT ElError ElSparseMatrixQueueZero_c
( ElSparseMatrix_c A, ElInt row, ElInt col );
EL_EXPORT ElError ElSparseMatrixQueueZero_z
( ElSparseMatrix_z A, ElInt row, ElInt col );

/* void SparseMatrix<T>::ProcessQueues()
   -------------------------------------- */ 
EL_EXPORT ElError ElSparseMatrixProcessQueues_i( ElSparseMatrix_i A );
EL_EXPORT ElError ElSparseMatrixProcessQueues_s( ElSparseMatrix_s A );
EL_EXPORT ElError ElSparseMatrixProcessQueues_d( ElSparseMatrix_d A );
EL_EXPORT ElError ElSparseMatrixProcessQueues_c( ElSparseMatrix_c A );
EL_EXPORT ElError ElSparseMatrixProcessQueues_z( ElSparseMatrix_z A );

/* Queries
   ======= */

/* Int SparseMatrix<T>::Height() const
   ----------------------------------- */
EL_EXPORT ElError ElSparseMatrixHeight_i
( ElConstSparseMatrix_i A, ElInt* height );
EL_EXPORT ElError ElSparseMatrixHeight_s
( ElConstSparseMatrix_s A, ElInt* height );
EL_EXPORT ElError ElSparseMatrixHeight_d
( ElConstSparseMatrix_d A, ElInt* height );
EL_EXPORT ElError ElSparseMatrixHeight_c
( ElConstSparseMatrix_c A, ElInt* height );
EL_EXPORT ElError ElSparseMatrixHeight_z
( ElConstSparseMatrix_z A, ElInt* height );

/* Int SparseMatrix<T>::Width() const
   ---------------------------------- */
EL_EXPORT ElError ElSparseMatrixWidth_i
( ElConstSparseMatrix_i A, ElInt* width );
EL_EXPORT ElError ElSparseMatrixWidth_s
( ElConstSparseMatrix_s A, ElInt* width );
EL_EXPORT ElError ElSparseMatrixWidth_d
( ElConstSparseMatrix_d A, ElInt* width );
EL_EXPORT ElError ElSparseMatrixWidth_c
( ElConstSparseMatrix_c A, ElInt* width );
EL_EXPORT ElError ElSparseMatrixWidth_z
( ElConstSparseMatrix_z A, ElInt* width );

/* Int SparseMatrix<T>::NumEntries() const
   --------------------------------------- */
EL_EXPORT ElError ElSparseMatrixNumEntries_i
( ElConstSparseMatrix_i A, ElInt* numEntries );
EL_EXPORT ElError ElSparseMatrixNumEntries_s
( ElConstSparseMatrix_s A, ElInt* numEntries );
EL_EXPORT ElError ElSparseMatrixNumEntries_d
( ElConstSparseMatrix_d A, ElInt* numEntries );
EL_EXPORT ElError ElSparseMatrixNumEntries_c
( ElConstSparseMatrix_c A, ElInt* numEntries );
EL_EXPORT ElError ElSparseMatrixNumEntries_z
( ElConstSparseMatrix_z A, ElInt* numEntries );

/* Int SparseMatrix<T>::Capacity() const
   ------------------------------------- */
EL_EXPORT ElError ElSparseMatrixCapacity_i
( ElConstSparseMatrix_i A, ElInt* capacity );
EL_EXPORT ElError ElSparseMatrixCapacity_s
( ElConstSparseMatrix_s A, ElInt* capacity );
EL_EXPORT ElError ElSparseMatrixCapacity_d
( ElConstSparseMatrix_d A, ElInt* capacity );
EL_EXPORT ElError ElSparseMatrixCapacity_c
( ElConstSparseMatrix_c A, ElInt* capacity );
EL_EXPORT ElError ElSparseMatrixCapacity_z
( ElConstSparseMatrix_z A, ElInt* capacity );

/* bool SparseMatrix<T>::Consistent() const
   ---------------------------------------- */
EL_EXPORT ElError ElSparseMatrixConsistent_i
( ElConstSparseMatrix_i A, bool* consistent );
EL_EXPORT ElError ElSparseMatrixConsistent_s
( ElConstSparseMatrix_s A, bool* consistent );
EL_EXPORT ElError ElSparseMatrixConsistent_d
( ElConstSparseMatrix_d A, bool* consistent );
EL_EXPORT ElError ElSparseMatrixConsistent_c
( ElConstSparseMatrix_c A, bool* consistent );
EL_EXPORT ElError ElSparseMatrixConsistent_z
( ElConstSparseMatrix_z A, bool* consistent );

/* Graph& SparseMatrix<T>::Graph()
   ------------------------------- */
EL_EXPORT ElError ElSparseMatrixGraph_i( ElSparseMatrix_i A, ElGraph* graph );
EL_EXPORT ElError ElSparseMatrixGraph_s( ElSparseMatrix_s A, ElGraph* graph );
EL_EXPORT ElError ElSparseMatrixGraph_d( ElSparseMatrix_d A, ElGraph* graph );
EL_EXPORT ElError ElSparseMatrixGraph_c( ElSparseMatrix_c A, ElGraph* graph );
EL_EXPORT ElError ElSparseMatrixGraph_z( ElSparseMatrix_z A, ElGraph* graph );

/* const Graph& SparseMatrix<T>::LockedGraph() const
   ------------------------------------------------- */
EL_EXPORT ElError ElSparseMatrixLockedGraph_i
( ElConstSparseMatrix_i A, ElConstGraph* graph );
EL_EXPORT ElError ElSparseMatrixLockedGraph_s
( ElConstSparseMatrix_s A, ElConstGraph* graph );
EL_EXPORT ElError ElSparseMatrixLockedGraph_d
( ElConstSparseMatrix_d A, ElConstGraph* graph );
EL_EXPORT ElError ElSparseMatrixLockedGraph_c
( ElConstSparseMatrix_c A, ElConstGraph* graph );
EL_EXPORT ElError ElSparseMatrixLockedGraph_z
( ElConstSparseMatrix_z A, ElConstGraph* graph );

/* Int SparseMatrix<T>::Row( Int index ) const
   ------------------------------------------- */
EL_EXPORT ElError ElSparseMatrixRow_i
( ElConstSparseMatrix_i A, ElInt index, ElInt* row );
EL_EXPORT ElError ElSparseMatrixRow_s
( ElConstSparseMatrix_s A, ElInt index, ElInt* row );
EL_EXPORT ElError ElSparseMatrixRow_d
( ElConstSparseMatrix_d A, ElInt index, ElInt* row );
EL_EXPORT ElError ElSparseMatrixRow_c
( ElConstSparseMatrix_c A, ElInt index, ElInt* row );
EL_EXPORT ElError ElSparseMatrixRow_z
( ElConstSparseMatrix_z A, ElInt index, ElInt* row );

/* Int SparseMatrix<T>::Col( Int index ) const
   ------------------------------------------- */
EL_EXPORT ElError ElSparseMatrixCol_i
( ElConstSparseMatrix_i A, ElInt index, ElInt* col );
EL_EXPORT ElError ElSparseMatrixCol_s
( ElConstSparseMatrix_s A, ElInt index, ElInt* col );
EL_EXPORT ElError ElSparseMatrixCol_d
( ElConstSparseMatrix_d A, ElInt index, ElInt* col );
EL_EXPORT ElError ElSparseMatrixCol_c
( ElConstSparseMatrix_c A, ElInt index, ElInt* col );
EL_EXPORT ElError ElSparseMatrixCol_z
( ElConstSparseMatrix_z A, ElInt index, ElInt* col );

/* T SparseMatrix<T>::Value( Int index ) const
   ------------------------------------------- */
EL_EXPORT ElError ElSparseMatrixValue_i
( ElConstSparseMatrix_i A, ElInt index, ElInt* value );
EL_EXPORT ElError ElSparseMatrixValue_s
( ElConstSparseMatrix_s A, ElInt index, float* value );
EL_EXPORT ElError ElSparseMatrixValue_d
( ElConstSparseMatrix_d A, ElInt index, double* value );
EL_EXPORT ElError ElSparseMatrixValue_c
( ElConstSparseMatrix_c A, ElInt index, complex_float* value );
EL_EXPORT ElError ElSparseMatrixValue_z
( ElConstSparseMatrix_z A, ElInt index, complex_double* value );

/* Int SparseMatrix<T>::RowOffset( Int row ) const
   ----------------------------------------------- */
EL_EXPORT ElError ElSparseMatrixRowOffset_i
( ElConstSparseMatrix_i A, ElInt row, ElInt* rowOffset );
EL_EXPORT ElError ElSparseMatrixRowOffset_s
( ElConstSparseMatrix_s A, ElInt row, ElInt* rowOffset );
EL_EXPORT ElError ElSparseMatrixRowOffset_d
( ElConstSparseMatrix_d A, ElInt row, ElInt* rowOffset );
EL_EXPORT ElError ElSparseMatrixRowOffset_c
( ElConstSparseMatrix_c A, ElInt row, ElInt* rowOffset );
EL_EXPORT ElError ElSparseMatrixRowOffset_z
( ElConstSparseMatrix_z A, ElInt row, ElInt* rowOffset );

/* Int SparseMatrix<T>::Offset( Int row, Int col ) const
   ----------------------------------------------------- */
EL_EXPORT ElError ElSparseMatrixOffset_i
( ElConstSparseMatrix_i A, ElInt row, ElInt col, ElInt* offset );
EL_EXPORT ElError ElSparseMatrixOffset_s
( ElConstSparseMatrix_s A, ElInt row, ElInt col, ElInt* offset );
EL_EXPORT ElError ElSparseMatrixOffset_d
( ElConstSparseMatrix_d A, ElInt row, ElInt col, ElInt* offset );
EL_EXPORT ElError ElSparseMatrixOffset_c
( ElConstSparseMatrix_c A, ElInt row, ElInt col, ElInt* offset );
EL_EXPORT ElError ElSparseMatrixOffset_z
( ElConstSparseMatrix_z A, ElInt row, ElInt col, ElInt* offset );

/* Int SparseMatrix<T>::NumConnections( Int row ) const
   ---------------------------------------------------- */
EL_EXPORT ElError ElSparseMatrixNumConnections_i
( ElConstSparseMatrix_i A, ElInt row, ElInt* numConnections );
EL_EXPORT ElError ElSparseMatrixNumConnections_s
( ElConstSparseMatrix_s A, ElInt row, ElInt* numConnections );
EL_EXPORT ElError ElSparseMatrixNumConnections_d
( ElConstSparseMatrix_d A, ElInt row, ElInt* numConnections );
EL_EXPORT ElError ElSparseMatrixNumConnections_c
( ElConstSparseMatrix_c A, ElInt row, ElInt* numConnections );
EL_EXPORT ElError ElSparseMatrixNumConnections_z
( ElConstSparseMatrix_z A, ElInt row, ElInt* numConnections );

/* Int* SparseMatrix<T>::SourceBuffer()
   ------------------------------------ */
EL_EXPORT ElError ElSparseMatrixSourceBuffer_i
( ElSparseMatrix_i A, ElInt** sourceBuffer );
EL_EXPORT ElError ElSparseMatrixSourceBuffer_s
( ElSparseMatrix_s A, ElInt** sourceBuffer );
EL_EXPORT ElError ElSparseMatrixSourceBuffer_d
( ElSparseMatrix_d A, ElInt** sourceBuffer );
EL_EXPORT ElError ElSparseMatrixSourceBuffer_c
( ElSparseMatrix_c A, ElInt** sourceBuffer );
EL_EXPORT ElError ElSparseMatrixSourceBuffer_z
( ElSparseMatrix_z A, ElInt** sourceBuffer );

/* const Int* SparseMatrix<T>::LockedSourceBuffer() const
   ------------------------------------------------------ */
EL_EXPORT ElError ElSparseMatrixLockedSourceBuffer_i
( ElConstSparseMatrix_i A, const ElInt** sourceBuffer );
EL_EXPORT ElError ElSparseMatrixLockedSourceBuffer_s
( ElConstSparseMatrix_s A, const ElInt** sourceBuffer );
EL_EXPORT ElError ElSparseMatrixLockedSourceBuffer_d
( ElConstSparseMatrix_d A, const ElInt** sourceBuffer );
EL_EXPORT ElError ElSparseMatrixLockedSourceBuffer_c
( ElConstSparseMatrix_c A, const ElInt** sourceBuffer );
EL_EXPORT ElError ElSparseMatrixLockedSourceBuffer_z
( ElConstSparseMatrix_z A, const ElInt** sourceBuffer );

/* Int* SparseMatrix<T>::TargetBuffer()
   ------------------------------------ */
EL_EXPORT ElError ElSparseMatrixTargetBuffer_i
( ElSparseMatrix_i A, ElInt** targetBuffer );
EL_EXPORT ElError ElSparseMatrixTargetBuffer_s
( ElSparseMatrix_s A, ElInt** targetBuffer );
EL_EXPORT ElError ElSparseMatrixTargetBuffer_d
( ElSparseMatrix_d A, ElInt** targetBuffer );
EL_EXPORT ElError ElSparseMatrixTargetBuffer_c
( ElSparseMatrix_c A, ElInt** targetBuffer );
EL_EXPORT ElError ElSparseMatrixTargetBuffer_z
( ElSparseMatrix_z A, ElInt** targetBuffer );

/* const Int* SparseMatrix<T>::LockedTargetBuffer() const
   ------------------------------------------------------ */
EL_EXPORT ElError ElSparseMatrixLockedTargetBuffer_i
( ElConstSparseMatrix_i A, const ElInt** targetBuffer );
EL_EXPORT ElError ElSparseMatrixLockedTargetBuffer_s
( ElConstSparseMatrix_s A, const ElInt** targetBuffer );
EL_EXPORT ElError ElSparseMatrixLockedTargetBuffer_d
( ElConstSparseMatrix_d A, const ElInt** targetBuffer );
EL_EXPORT ElError ElSparseMatrixLockedTargetBuffer_c
( ElConstSparseMatrix_c A, const ElInt** targetBuffer );
EL_EXPORT ElError ElSparseMatrixLockedTargetBuffer_z
( ElConstSparseMatrix_z A, const ElInt** targetBuffer );

/* T* SparseMatrix<T>::ValueBuffer()
   --------------------------------- */
EL_EXPORT ElError ElSparseMatrixValueBuffer_i
( ElSparseMatrix_i A, ElInt** valueBuffer );
EL_EXPORT ElError ElSparseMatrixValueBuffer_s
( ElSparseMatrix_s A, float** valueBuffer );
EL_EXPORT ElError ElSparseMatrixValueBuffer_d
( ElSparseMatrix_d A, double** valueBuffer );
EL_EXPORT ElError ElSparseMatrixValueBuffer_c
( ElSparseMatrix_c A, complex_float** valueBuffer );
EL_EXPORT ElError ElSparseMatrixValueBuffer_z
( ElSparseMatrix_z A, complex_double** valueBuffer );

/* const T* SparseMatrix<T>::LockedValueBuffer() const
   --------------------------------------------------- */
EL_EXPORT ElError ElSparseMatrixLockedValueBuffer_i
( ElConstSparseMatrix_i A, const ElInt** valueBuffer );
EL_EXPORT ElError ElSparseMatrixLockedValueBuffer_s
( ElConstSparseMatrix_s A, const float** valueBuffer );
EL_EXPORT ElError ElSparseMatrixLockedValueBuffer_d
( ElConstSparseMatrix_d A, const double** valueBuffer );
EL_EXPORT ElError ElSparseMatrixLockedValueBuffer_c
( ElConstSparseMatrix_c A, const complex_float** valueBuffer );
EL_EXPORT ElError ElSparseMatrixLockedValueBuffer_z
( ElConstSparseMatrix_z A, const complex_double** valueBuffer );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_SPARSEMATRIX_C_H */
