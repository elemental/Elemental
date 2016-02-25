/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_DISTSPARSEMATRIX_C_H
#define EL_DISTSPARSEMATRIX_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* An anonymous struct meant as a placeholder for DistSparseMatrix<T>
   ------------------------------------------------------------------ */
typedef struct ElDistSparseMatrix_iDummy* ElDistSparseMatrix_i;
typedef struct ElDistSparseMatrix_sDummy* ElDistSparseMatrix_s;
typedef struct ElDistSparseMatrix_dDummy* ElDistSparseMatrix_d;
typedef struct ElDistSparseMatrix_cDummy* ElDistSparseMatrix_c;
typedef struct ElDistSparseMatrix_zDummy* ElDistSparseMatrix_z;

typedef const struct ElDistSparseMatrix_iDummy* ElConstDistSparseMatrix_i;
typedef const struct ElDistSparseMatrix_sDummy* ElConstDistSparseMatrix_s;
typedef const struct ElDistSparseMatrix_dDummy* ElConstDistSparseMatrix_d;
typedef const struct ElDistSparseMatrix_cDummy* ElConstDistSparseMatrix_c;
typedef const struct ElDistSparseMatrix_zDummy* ElConstDistSparseMatrix_z;

/* Constructors and destructors
   ============================ */

/* DistSparseMatrix<T>::DistSparseMatrix()
   --------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixCreate_i
( ElDistSparseMatrix_i* A, MPI_Comm comm );
EL_EXPORT ElError ElDistSparseMatrixCreate_s
( ElDistSparseMatrix_s* A, MPI_Comm comm );
EL_EXPORT ElError ElDistSparseMatrixCreate_d
( ElDistSparseMatrix_d* A, MPI_Comm comm );
EL_EXPORT ElError ElDistSparseMatrixCreate_c
( ElDistSparseMatrix_c* A, MPI_Comm comm );
EL_EXPORT ElError ElDistSparseMatrixCreate_z
( ElDistSparseMatrix_z* A, MPI_Comm comm );

/* DistSparseMatrix<T>::~DistSparseMatrix()
   ---------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixDestroy_i( ElConstDistSparseMatrix_i A );
EL_EXPORT ElError ElDistSparseMatrixDestroy_s( ElConstDistSparseMatrix_s A );
EL_EXPORT ElError ElDistSparseMatrixDestroy_d( ElConstDistSparseMatrix_d A );
EL_EXPORT ElError ElDistSparseMatrixDestroy_c( ElConstDistSparseMatrix_c A );
EL_EXPORT ElError ElDistSparseMatrixDestroy_z( ElConstDistSparseMatrix_z A );

/* Assignment and reconfiguration
   ============================== */

/* void DistSparseMatrix<T>::Empty()
   --------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixEmpty_i( ElDistSparseMatrix_i A );
EL_EXPORT ElError ElDistSparseMatrixEmpty_s( ElDistSparseMatrix_s A );
EL_EXPORT ElError ElDistSparseMatrixEmpty_d( ElDistSparseMatrix_d A );
EL_EXPORT ElError ElDistSparseMatrixEmpty_c( ElDistSparseMatrix_c A );
EL_EXPORT ElError ElDistSparseMatrixEmpty_z( ElDistSparseMatrix_z A );

/* void DistSparseMatrix<T>::Resize( Int height, Int width )
   --------------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixResize_i
( ElDistSparseMatrix_i A, ElInt height, ElInt width );
EL_EXPORT ElError ElDistSparseMatrixResize_s
( ElDistSparseMatrix_s A, ElInt height, ElInt width );
EL_EXPORT ElError ElDistSparseMatrixResize_d
( ElDistSparseMatrix_d A, ElInt height, ElInt width );
EL_EXPORT ElError ElDistSparseMatrixResize_c
( ElDistSparseMatrix_c A, ElInt height, ElInt width );
EL_EXPORT ElError ElDistSparseMatrixResize_z
( ElDistSparseMatrix_z A, ElInt height, ElInt width );

/* void DistSparseMatrix<T>::SetComm( mpi::Comm comm )
   --------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixSetComm_i
( ElDistSparseMatrix_i A, MPI_Comm comm );
EL_EXPORT ElError ElDistSparseMatrixSetComm_s
( ElDistSparseMatrix_s A, MPI_Comm comm );
EL_EXPORT ElError ElDistSparseMatrixSetComm_d
( ElDistSparseMatrix_d A, MPI_Comm comm );
EL_EXPORT ElError ElDistSparseMatrixSetComm_c
( ElDistSparseMatrix_c A, MPI_Comm comm );
EL_EXPORT ElError ElDistSparseMatrixSetComm_z
( ElDistSparseMatrix_z A, MPI_Comm comm );

/* void DistSparseMatrix<T>::Reserve
   ( Int numLocalEntries, Int numRemoteEntries )
  ---------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixReserve_i
( ElDistSparseMatrix_i A, ElInt numLocalEntries, ElInt numRemoteEntries );
EL_EXPORT ElError ElDistSparseMatrixReserve_s
( ElDistSparseMatrix_s A, ElInt numLocalEntries, ElInt numRemoteEntries );
EL_EXPORT ElError ElDistSparseMatrixReserve_d
( ElDistSparseMatrix_d A, ElInt numLocalEntries, ElInt numRemoteEntries );
EL_EXPORT ElError ElDistSparseMatrixReserve_c
( ElDistSparseMatrix_c A, ElInt numLocalEntries, ElInt numRemoteEntries );
EL_EXPORT ElError ElDistSparseMatrixReserve_z
( ElDistSparseMatrix_z A, ElInt numLocalEntries, ElInt numRemoteEntries );

/* void DistSparseMatrix<T>::Update( Int row, Int col, T value )
   -------------------------------------------------------------*/
EL_EXPORT ElError ElDistSparseMatrixUpdate_i
( ElDistSparseMatrix_i A, ElInt row, ElInt col, ElInt value );
EL_EXPORT ElError ElDistSparseMatrixUpdate_s
( ElDistSparseMatrix_s A, ElInt row, ElInt col, float value );
EL_EXPORT ElError ElDistSparseMatrixUpdate_d
( ElDistSparseMatrix_d A, ElInt row, ElInt col, double value );
EL_EXPORT ElError ElDistSparseMatrixUpdate_c
( ElDistSparseMatrix_c A, ElInt row, ElInt col, complex_float value );
EL_EXPORT ElError ElDistSparseMatrixUpdate_z
( ElDistSparseMatrix_z A, ElInt row, ElInt col, complex_double value );

/* void DistSparseMatrix<T>::UpdateLocal( Int localRow, Int col, T value )
   ----------------------------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixUpdateLocal_i
( ElDistSparseMatrix_i A, ElInt localRow, ElInt col, ElInt value );
EL_EXPORT ElError ElDistSparseMatrixUpdateLocal_s
( ElDistSparseMatrix_s A, ElInt localRow, ElInt col, float value );
EL_EXPORT ElError ElDistSparseMatrixUpdateLocal_d
( ElDistSparseMatrix_d A, ElInt localRow, ElInt col, double value );
EL_EXPORT ElError ElDistSparseMatrixUpdateLocal_c
( ElDistSparseMatrix_c A, ElInt localRow, ElInt col, complex_float value );
EL_EXPORT ElError ElDistSparseMatrixUpdateLocal_z
( ElDistSparseMatrix_z A, ElInt localRow, ElInt col, complex_double value );

/* void DistSparseMatrix<T>::Zero( Int row, Int col )
   -------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixZero_i
( ElDistSparseMatrix_i A, ElInt row, ElInt col );
EL_EXPORT ElError ElDistSparseMatrixZero_s
( ElDistSparseMatrix_s A, ElInt row, ElInt col );
EL_EXPORT ElError ElDistSparseMatrixZero_d
( ElDistSparseMatrix_d A, ElInt row, ElInt col );
EL_EXPORT ElError ElDistSparseMatrixZero_c
( ElDistSparseMatrix_c A, ElInt row, ElInt col );
EL_EXPORT ElError ElDistSparseMatrixZero_z
( ElDistSparseMatrix_z A, ElInt row, ElInt col );

/* void DistSparseMatrix<T>::ZeroLocal( Int localRow, Int col )
   ------------------------------------------------------------ */
EL_EXPORT ElError ElDistSparseMatrixZeroLocal_i
( ElDistSparseMatrix_i A, ElInt localRow, ElInt col );
EL_EXPORT ElError ElDistSparseMatrixZeroLocal_s
( ElDistSparseMatrix_s A, ElInt localRow, ElInt col );
EL_EXPORT ElError ElDistSparseMatrixZeroLocal_d
( ElDistSparseMatrix_d A, ElInt localRow, ElInt col );
EL_EXPORT ElError ElDistSparseMatrixZeroLocal_c
( ElDistSparseMatrix_c A, ElInt localRow, ElInt col );
EL_EXPORT ElError ElDistSparseMatrixZeroLocal_z
( ElDistSparseMatrix_z A, ElInt localRow, ElInt col );

/* void DistSparseMatrix<T>::QueueUpdate
   ( Int row, Int col, T value, bool passive )
   ------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixQueueUpdate_i
( ElDistSparseMatrix_i A, 
  ElInt row, ElInt col, ElInt value, bool passive );
EL_EXPORT ElError ElDistSparseMatrixQueueUpdate_s
( ElDistSparseMatrix_s A, 
  ElInt row, ElInt col, float value, bool passive );
EL_EXPORT ElError ElDistSparseMatrixQueueUpdate_d
( ElDistSparseMatrix_d A, 
  ElInt row, ElInt col, double value, bool passive );
EL_EXPORT ElError ElDistSparseMatrixQueueUpdate_c
( ElDistSparseMatrix_c A, 
  ElInt row, ElInt col, complex_float value, bool passive );
EL_EXPORT ElError ElDistSparseMatrixQueueUpdate_z
( ElDistSparseMatrix_z A, 
  ElInt row, ElInt col, complex_double value, bool passive );

/* void DistSparseMatrix<T>::QueueLocalUpdate
   ( Int localRow, Int col, T value )
   ------------------------------------------ */
EL_EXPORT ElError ElDistSparseMatrixQueueLocalUpdate_i
( ElDistSparseMatrix_i A, ElInt localRow, ElInt col, ElInt value );
EL_EXPORT ElError ElDistSparseMatrixQueueLocalUpdate_s
( ElDistSparseMatrix_s A, ElInt localRow, ElInt col, float value );
EL_EXPORT ElError ElDistSparseMatrixQueueLocalUpdate_d
( ElDistSparseMatrix_d A, ElInt localRow, ElInt col, double value );
EL_EXPORT ElError ElDistSparseMatrixQueueLocalUpdate_c
( ElDistSparseMatrix_c A, ElInt localRow, ElInt col, complex_float value );
EL_EXPORT ElError ElDistSparseMatrixQueueLocalUpdate_z
( ElDistSparseMatrix_z A, ElInt localRow, ElInt col, complex_double value );

/* void DistSparseMatrix<T>::QueueZero( Int row, Int col, bool passive )
   --------------------------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixQueueZero_i
( ElDistSparseMatrix_i A, ElInt row, ElInt col, bool passive );
EL_EXPORT ElError ElDistSparseMatrixQueueZero_s
( ElDistSparseMatrix_s A, ElInt row, ElInt col, bool passive );
EL_EXPORT ElError ElDistSparseMatrixQueueZero_d
( ElDistSparseMatrix_d A, ElInt row, ElInt col, bool passive );
EL_EXPORT ElError ElDistSparseMatrixQueueZero_c
( ElDistSparseMatrix_c A, ElInt row, ElInt col, bool passive );
EL_EXPORT ElError ElDistSparseMatrixQueueZero_z
( ElDistSparseMatrix_z A, ElInt row, ElInt col, bool passive );

/* void DistSparseMatrix<T>::QueueLocalZero( Int localRow, Int col )
   ----------------------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixQueueLocalZero_i
( ElDistSparseMatrix_i A, ElInt localRow, ElInt col );
EL_EXPORT ElError ElDistSparseMatrixQueueLocalZero_s
( ElDistSparseMatrix_s A, ElInt localRow, ElInt col );
EL_EXPORT ElError ElDistSparseMatrixQueueLocalZero_d
( ElDistSparseMatrix_d A, ElInt localRow, ElInt col );
EL_EXPORT ElError ElDistSparseMatrixQueueLocalZero_c
( ElDistSparseMatrix_c A, ElInt localRow, ElInt col );
EL_EXPORT ElError ElDistSparseMatrixQueueLocalZero_z
( ElDistSparseMatrix_z A, ElInt localRow, ElInt col );

/* void DistSparseMatrix<T>::ProcessQueues()
   ----------------------------------------- */ 
EL_EXPORT ElError ElDistSparseMatrixProcessQueues_i( ElDistSparseMatrix_i A );
EL_EXPORT ElError ElDistSparseMatrixProcessQueues_s( ElDistSparseMatrix_s A );
EL_EXPORT ElError ElDistSparseMatrixProcessQueues_d( ElDistSparseMatrix_d A );
EL_EXPORT ElError ElDistSparseMatrixProcessQueues_c( ElDistSparseMatrix_c A );
EL_EXPORT ElError ElDistSparseMatrixProcessQueues_z( ElDistSparseMatrix_z A );

/* void DistSparseMatrix<T>::ProcessLocalQueues()
   ---------------------------------------------- */
EL_EXPORT ElError 
ElDistSparseMatrixProcessLocalQueues_i( ElDistSparseMatrix_i A );
EL_EXPORT ElError 
ElDistSparseMatrixProcessLocalQueues_s( ElDistSparseMatrix_s A );
EL_EXPORT ElError 
ElDistSparseMatrixProcessLocalQueues_d( ElDistSparseMatrix_d A );
EL_EXPORT ElError 
ElDistSparseMatrixProcessLocalQueues_c( ElDistSparseMatrix_c A );
EL_EXPORT ElError 
ElDistSparseMatrixProcessLocalQueues_z( ElDistSparseMatrix_z A );

/* Queries
   ======= */

/* Int DistSparseMatrix<T>::Height() const
   --------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixHeight_i
( ElConstDistSparseMatrix_i A, ElInt* height );
EL_EXPORT ElError ElDistSparseMatrixHeight_s
( ElConstDistSparseMatrix_s A, ElInt* height );
EL_EXPORT ElError ElDistSparseMatrixHeight_d
( ElConstDistSparseMatrix_d A, ElInt* height );
EL_EXPORT ElError ElDistSparseMatrixHeight_c
( ElConstDistSparseMatrix_c A, ElInt* height );
EL_EXPORT ElError ElDistSparseMatrixHeight_z
( ElConstDistSparseMatrix_z A, ElInt* height );

/* Int DistSparseMatrix<T>::Width() const
   -------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixWidth_i
( ElConstDistSparseMatrix_i A, ElInt* width );
EL_EXPORT ElError ElDistSparseMatrixWidth_s
( ElConstDistSparseMatrix_s A, ElInt* width );
EL_EXPORT ElError ElDistSparseMatrixWidth_d
( ElConstDistSparseMatrix_d A, ElInt* width );
EL_EXPORT ElError ElDistSparseMatrixWidth_c
( ElConstDistSparseMatrix_c A, ElInt* width );
EL_EXPORT ElError ElDistSparseMatrixWidth_z
( ElConstDistSparseMatrix_z A, ElInt* width );

/* DistGraph& DistSparseMatrix<T>::DistGraph()
   ------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixDistGraph_i
( ElDistSparseMatrix_i A, ElDistGraph* graph );
EL_EXPORT ElError ElDistSparseMatrixDistGraph_s
( ElDistSparseMatrix_s A, ElDistGraph* graph );
EL_EXPORT ElError ElDistSparseMatrixDistGraph_d
( ElDistSparseMatrix_d A, ElDistGraph* graph );
EL_EXPORT ElError ElDistSparseMatrixDistGraph_c
( ElDistSparseMatrix_c A, ElDistGraph* graph );
EL_EXPORT ElError ElDistSparseMatrixDistGraph_z
( ElDistSparseMatrix_z A, ElDistGraph* graph );

/* const DistGraph& DistSparseMatrix<T>::LockedDistGraph() const
   ------------------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixLockedDistGraph_i
( ElConstDistSparseMatrix_i A, ElConstDistGraph* graph );
EL_EXPORT ElError ElDistSparseMatrixLockedDistGraph_s
( ElConstDistSparseMatrix_s A, ElConstDistGraph* graph );
EL_EXPORT ElError ElDistSparseMatrixLockedDistGraph_d
( ElConstDistSparseMatrix_d A, ElConstDistGraph* graph );
EL_EXPORT ElError ElDistSparseMatrixLockedDistGraph_c
( ElConstDistSparseMatrix_c A, ElConstDistGraph* graph );
EL_EXPORT ElError ElDistSparseMatrixLockedDistGraph_z
( ElConstDistSparseMatrix_z A, ElConstDistGraph* graph );

/* Int DistSparseMatrix<T>::FirstLocalRow() const
   ---------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixFirstLocalRow_i
( ElConstDistSparseMatrix_i A, ElInt* firstLocalRow );
EL_EXPORT ElError ElDistSparseMatrixFirstLocalRow_s
( ElConstDistSparseMatrix_s A, ElInt* firstLocalRow );
EL_EXPORT ElError ElDistSparseMatrixFirstLocalRow_d
( ElConstDistSparseMatrix_d A, ElInt* firstLocalRow );
EL_EXPORT ElError ElDistSparseMatrixFirstLocalRow_c
( ElConstDistSparseMatrix_c A, ElInt* firstLocalRow );
EL_EXPORT ElError ElDistSparseMatrixFirstLocalRow_z
( ElConstDistSparseMatrix_z A, ElInt* firstLocalRow );

/* Int DistSparseMatrix<T>::LocalHeight() const
   -------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixLocalHeight_i
( ElConstDistSparseMatrix_i A, ElInt* localHeight );
EL_EXPORT ElError ElDistSparseMatrixLocalHeight_s
( ElConstDistSparseMatrix_s A, ElInt* localHeight );
EL_EXPORT ElError ElDistSparseMatrixLocalHeight_d
( ElConstDistSparseMatrix_d A, ElInt* localHeight );
EL_EXPORT ElError ElDistSparseMatrixLocalHeight_c
( ElConstDistSparseMatrix_c A, ElInt* localHeight );
EL_EXPORT ElError ElDistSparseMatrixLocalHeight_z
( ElConstDistSparseMatrix_z A, ElInt* localHeight );

/* Int DistSparseMatrix<T>::NumLocalEntries() const
   ------------------------------------------------ */
EL_EXPORT ElError ElDistSparseMatrixNumLocalEntries_i
( ElConstDistSparseMatrix_i A, ElInt* numLocalEntries );
EL_EXPORT ElError ElDistSparseMatrixNumLocalEntries_s
( ElConstDistSparseMatrix_s A, ElInt* numLocalEntries );
EL_EXPORT ElError ElDistSparseMatrixNumLocalEntries_d
( ElConstDistSparseMatrix_d A, ElInt* numLocalEntries );
EL_EXPORT ElError ElDistSparseMatrixNumLocalEntries_c
( ElConstDistSparseMatrix_c A, ElInt* numLocalEntries );
EL_EXPORT ElError ElDistSparseMatrixNumLocalEntries_z
( ElConstDistSparseMatrix_z A, ElInt* numLocalEntries );

/* Int DistSparseMatrix<T>::Capacity() const
   ----------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixCapacity_i
( ElConstDistSparseMatrix_i A, ElInt* capacity );
EL_EXPORT ElError ElDistSparseMatrixCapacity_s
( ElConstDistSparseMatrix_s A, ElInt* capacity );
EL_EXPORT ElError ElDistSparseMatrixCapacity_d
( ElConstDistSparseMatrix_d A, ElInt* capacity );
EL_EXPORT ElError ElDistSparseMatrixCapacity_c
( ElConstDistSparseMatrix_c A, ElInt* capacity );
EL_EXPORT ElError ElDistSparseMatrixCapacity_z
( ElConstDistSparseMatrix_z A, ElInt* capacity );

/* bool DistSparseMatrix<T>::Consistent() const
   -------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixConsistent_i
( ElConstDistSparseMatrix_i A, bool* consistent );
EL_EXPORT ElError ElDistSparseMatrixConsistent_s
( ElConstDistSparseMatrix_s A, bool* consistent );
EL_EXPORT ElError ElDistSparseMatrixConsistent_d
( ElConstDistSparseMatrix_d A, bool* consistent );
EL_EXPORT ElError ElDistSparseMatrixConsistent_c
( ElConstDistSparseMatrix_c A, bool* consistent );
EL_EXPORT ElError ElDistSparseMatrixConsistent_z
( ElConstDistSparseMatrix_z A, bool* consistent );

/* mpi::Comm DistSparseMatrix<T>::Comm() const
   ------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixComm_i
( ElConstDistSparseMatrix_i A, MPI_Comm* comm );
EL_EXPORT ElError ElDistSparseMatrixComm_s
( ElConstDistSparseMatrix_s A, MPI_Comm* comm );
EL_EXPORT ElError ElDistSparseMatrixComm_d
( ElConstDistSparseMatrix_d A, MPI_Comm* comm );
EL_EXPORT ElError ElDistSparseMatrixComm_c
( ElConstDistSparseMatrix_c A, MPI_Comm* comm );
EL_EXPORT ElError ElDistSparseMatrixComm_z
( ElConstDistSparseMatrix_z A, MPI_Comm* comm );

/* Int DistSparseMatrix<T>::Blocksize() const
   ------------------------------------------ */
EL_EXPORT ElError ElDistSparseMatrixBlocksize_i
( ElConstDistSparseMatrix_i A, ElInt* blocksize );
EL_EXPORT ElError ElDistSparseMatrixBlocksize_s
( ElConstDistSparseMatrix_s A, ElInt* blocksize );
EL_EXPORT ElError ElDistSparseMatrixBlocksize_d
( ElConstDistSparseMatrix_d A, ElInt* blocksize );
EL_EXPORT ElError ElDistSparseMatrixBlocksize_c
( ElConstDistSparseMatrix_c A, ElInt* blocksize );
EL_EXPORT ElError ElDistSparseMatrixBlocksize_z
( ElConstDistSparseMatrix_z A, ElInt* blocksize );

/* int DistSparseMatrix<T>::RowOwner( Int i ) const
   ------------------------------------------------ */
EL_EXPORT ElError ElDistSparseMatrixRowOwner_i
( ElConstDistSparseMatrix_i A, ElInt i, int* owner );
EL_EXPORT ElError ElDistSparseMatrixRowOwner_s
( ElConstDistSparseMatrix_s A, ElInt i, int* owner );
EL_EXPORT ElError ElDistSparseMatrixRowOwner_d
( ElConstDistSparseMatrix_d A, ElInt i, int* owner );
EL_EXPORT ElError ElDistSparseMatrixRowOwner_c
( ElConstDistSparseMatrix_c A, ElInt i, int* owner );
EL_EXPORT ElError ElDistSparseMatrixRowOwner_z
( ElConstDistSparseMatrix_z A, ElInt i, int* owner );

/* Int DistSparseMatrix<T>::GlobalRow( Int iLoc ) const
   ---------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixGlobalRow_i
( ElConstDistSparseMatrix_i A, ElInt iLoc, ElInt* i );
EL_EXPORT ElError ElDistSparseMatrixGlobalRow_s
( ElConstDistSparseMatrix_s A, ElInt iLoc, ElInt* i );
EL_EXPORT ElError ElDistSparseMatrixGlobalRow_d
( ElConstDistSparseMatrix_d A, ElInt iLoc, ElInt* i );
EL_EXPORT ElError ElDistSparseMatrixGlobalRow_c
( ElConstDistSparseMatrix_c A, ElInt iLoc, ElInt* i );
EL_EXPORT ElError ElDistSparseMatrixGlobalRow_z
( ElConstDistSparseMatrix_z A, ElInt iLoc, ElInt* i );

/* Int DistSparseMatrix<T>::Row( Int localInd ) const
   -------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixRow_i
( ElConstDistSparseMatrix_i A, ElInt localInd, ElInt* row );
EL_EXPORT ElError ElDistSparseMatrixRow_s
( ElConstDistSparseMatrix_s A, ElInt localInd, ElInt* row );
EL_EXPORT ElError ElDistSparseMatrixRow_d
( ElConstDistSparseMatrix_d A, ElInt localInd, ElInt* row );
EL_EXPORT ElError ElDistSparseMatrixRow_c
( ElConstDistSparseMatrix_c A, ElInt localInd, ElInt* row );
EL_EXPORT ElError ElDistSparseMatrixRow_z
( ElConstDistSparseMatrix_z A, ElInt localInd, ElInt* row );

/* Int DistSparseMatrix<T>::Col( Int localInd ) const
   -------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixCol_i
( ElConstDistSparseMatrix_i A, ElInt localInd, ElInt* col );
EL_EXPORT ElError ElDistSparseMatrixCol_s
( ElConstDistSparseMatrix_s A, ElInt localInd, ElInt* col );
EL_EXPORT ElError ElDistSparseMatrixCol_d
( ElConstDistSparseMatrix_d A, ElInt localInd, ElInt* col );
EL_EXPORT ElError ElDistSparseMatrixCol_c
( ElConstDistSparseMatrix_c A, ElInt localInd, ElInt* col );
EL_EXPORT ElError ElDistSparseMatrixCol_z
( ElConstDistSparseMatrix_z A, ElInt localInd, ElInt* col );

/* T DistSparseMatrix<T>::Value( Int localInd ) const
   -------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixValue_i
( ElConstDistSparseMatrix_i A, ElInt localInd, ElInt* value );
EL_EXPORT ElError ElDistSparseMatrixValue_s
( ElConstDistSparseMatrix_s A, ElInt localInd, float* value );
EL_EXPORT ElError ElDistSparseMatrixValue_d
( ElConstDistSparseMatrix_d A, ElInt localInd, double* value );
EL_EXPORT ElError ElDistSparseMatrixValue_c
( ElConstDistSparseMatrix_c A, ElInt localInd, complex_float* value );
EL_EXPORT ElError ElDistSparseMatrixValue_z
( ElConstDistSparseMatrix_z A, ElInt localInd, complex_double* value );

/* Int DistSparseMatrix<T>::RowOffset( Int localRow ) const
   -------------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixRowOffset_i
( ElConstDistSparseMatrix_i A, ElInt localRow, ElInt* rowOffset );
EL_EXPORT ElError ElDistSparseMatrixRowOffset_s
( ElConstDistSparseMatrix_s A, ElInt localRow, ElInt* rowOffset );
EL_EXPORT ElError ElDistSparseMatrixRowOffset_d
( ElConstDistSparseMatrix_d A, ElInt localRow, ElInt* rowOffset );
EL_EXPORT ElError ElDistSparseMatrixRowOffset_c
( ElConstDistSparseMatrix_c A, ElInt localRow, ElInt* rowOffset );
EL_EXPORT ElError ElDistSparseMatrixRowOffset_z
( ElConstDistSparseMatrix_z A, ElInt localRow, ElInt* rowOffset );

/* Int DistSparseMatrix<T>::Offset( Int localRow, Int col ) const
   -------------------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixOffset_i
( ElConstDistSparseMatrix_i A, ElInt localRow, ElInt col, ElInt* offset );
EL_EXPORT ElError ElDistSparseMatrixOffset_s
( ElConstDistSparseMatrix_s A, ElInt localRow, ElInt col, ElInt* offset );
EL_EXPORT ElError ElDistSparseMatrixOffset_d
( ElConstDistSparseMatrix_d A, ElInt localRow, ElInt col, ElInt* offset );
EL_EXPORT ElError ElDistSparseMatrixOffset_c
( ElConstDistSparseMatrix_c A, ElInt localRow, ElInt col, ElInt* offset );
EL_EXPORT ElError ElDistSparseMatrixOffset_z
( ElConstDistSparseMatrix_z A, ElInt localRow, ElInt col, ElInt* offset );

/* Int DistSparseMatrix<T>::NumConnections( Int localRow ) const
   ------------------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixNumConnections_i
( ElConstDistSparseMatrix_i A, ElInt localRow, ElInt* numConnections );
EL_EXPORT ElError ElDistSparseMatrixNumConnections_s
( ElConstDistSparseMatrix_s A, ElInt localRow, ElInt* numConnections );
EL_EXPORT ElError ElDistSparseMatrixNumConnections_d
( ElConstDistSparseMatrix_d A, ElInt localRow, ElInt* numConnections );
EL_EXPORT ElError ElDistSparseMatrixNumConnections_c
( ElConstDistSparseMatrix_c A, ElInt localRow, ElInt* numConnections );
EL_EXPORT ElError ElDistSparseMatrixNumConnections_z
( ElConstDistSparseMatrix_z A, ElInt localRow, ElInt* numConnections );

/* double DistSparseMatrix<T>::Imbalance() const
   --------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixImbalance_i
( ElConstDistSparseMatrix_i A, double* imbalance );
EL_EXPORT ElError ElDistSparseMatrixImbalance_s
( ElConstDistSparseMatrix_s A, double* imbalance );
EL_EXPORT ElError ElDistSparseMatrixImbalance_d
( ElConstDistSparseMatrix_d A, double* imbalance );
EL_EXPORT ElError ElDistSparseMatrixImbalance_c
( ElConstDistSparseMatrix_c A, double* imbalance );
EL_EXPORT ElError ElDistSparseMatrixImbalance_z
( ElConstDistSparseMatrix_z A, double* imbalance );

/* Int* DistSparseMatrix<T>::SourceBuffer()
   ---------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixSourceBuffer_i
( ElDistSparseMatrix_i A, ElInt** sourceBuffer );
EL_EXPORT ElError ElDistSparseMatrixSourceBuffer_s
( ElDistSparseMatrix_s A, ElInt** sourceBuffer );
EL_EXPORT ElError ElDistSparseMatrixSourceBuffer_d
( ElDistSparseMatrix_d A, ElInt** sourceBuffer );
EL_EXPORT ElError ElDistSparseMatrixSourceBuffer_c
( ElDistSparseMatrix_c A, ElInt** sourceBuffer );
EL_EXPORT ElError ElDistSparseMatrixSourceBuffer_z
( ElDistSparseMatrix_z A, ElInt** sourceBuffer );

/* const Int* DistSparseMatrix<T>::LockedSourceBuffer() const
   ---------------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixLockedSourceBuffer_i
( ElConstDistSparseMatrix_i A, const ElInt** sourceBuffer );
EL_EXPORT ElError ElDistSparseMatrixLockedSourceBuffer_s
( ElConstDistSparseMatrix_s A, const ElInt** sourceBuffer );
EL_EXPORT ElError ElDistSparseMatrixLockedSourceBuffer_d
( ElConstDistSparseMatrix_d A, const ElInt** sourceBuffer );
EL_EXPORT ElError ElDistSparseMatrixLockedSourceBuffer_c
( ElConstDistSparseMatrix_c A, const ElInt** sourceBuffer );
EL_EXPORT ElError ElDistSparseMatrixLockedSourceBuffer_z
( ElConstDistSparseMatrix_z A, const ElInt** sourceBuffer );

/* Int* DistSparseMatrix<T>::TargetBuffer()
   ---------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixTargetBuffer_i
( ElDistSparseMatrix_i A, ElInt** targetBuffer );
EL_EXPORT ElError ElDistSparseMatrixTargetBuffer_s
( ElDistSparseMatrix_s A, ElInt** targetBuffer );
EL_EXPORT ElError ElDistSparseMatrixTargetBuffer_d
( ElDistSparseMatrix_d A, ElInt** targetBuffer );
EL_EXPORT ElError ElDistSparseMatrixTargetBuffer_c
( ElDistSparseMatrix_c A, ElInt** targetBuffer );
EL_EXPORT ElError ElDistSparseMatrixTargetBuffer_z
( ElDistSparseMatrix_z A, ElInt** targetBuffer );

/* const Int* DistSparseMatrix<T>::LockedTargetBuffer() const
   ---------------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixLockedTargetBuffer_i
( ElConstDistSparseMatrix_i A, const ElInt** targetBuffer );
EL_EXPORT ElError ElDistSparseMatrixLockedTargetBuffer_s
( ElConstDistSparseMatrix_s A, const ElInt** targetBuffer );
EL_EXPORT ElError ElDistSparseMatrixLockedTargetBuffer_d
( ElConstDistSparseMatrix_d A, const ElInt** targetBuffer );
EL_EXPORT ElError ElDistSparseMatrixLockedTargetBuffer_c
( ElConstDistSparseMatrix_c A, const ElInt** targetBuffer );
EL_EXPORT ElError ElDistSparseMatrixLockedTargetBuffer_z
( ElConstDistSparseMatrix_z A, const ElInt** targetBuffer );

/* T* DistSparseMatrix<T>::ValueBuffer()
   ------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixValueBuffer_i
( ElDistSparseMatrix_i A, ElInt** valueBuffer );
EL_EXPORT ElError ElDistSparseMatrixValueBuffer_s
( ElDistSparseMatrix_s A, float** valueBuffer );
EL_EXPORT ElError ElDistSparseMatrixValueBuffer_d
( ElDistSparseMatrix_d A, double** valueBuffer );
EL_EXPORT ElError ElDistSparseMatrixValueBuffer_c
( ElDistSparseMatrix_c A, complex_float** valueBuffer );
EL_EXPORT ElError ElDistSparseMatrixValueBuffer_z
( ElDistSparseMatrix_z A, complex_double** valueBuffer );

/* const T* DistSparseMatrix<T>::LockedValueBuffer() const
   ------------------------------------------------------- */
EL_EXPORT ElError ElDistSparseMatrixLockedValueBuffer_i
( ElConstDistSparseMatrix_i A, const ElInt** valueBuffer );
EL_EXPORT ElError ElDistSparseMatrixLockedValueBuffer_s
( ElConstDistSparseMatrix_s A, const float** valueBuffer );
EL_EXPORT ElError ElDistSparseMatrixLockedValueBuffer_d
( ElConstDistSparseMatrix_d A, const double** valueBuffer );
EL_EXPORT ElError ElDistSparseMatrixLockedValueBuffer_c
( ElConstDistSparseMatrix_c A, const complex_float** valueBuffer );
EL_EXPORT ElError ElDistSparseMatrixLockedValueBuffer_z
( ElConstDistSparseMatrix_z A, const complex_double** valueBuffer );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_DISTSPARSEMATRIX_C_H */
