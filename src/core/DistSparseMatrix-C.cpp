/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El.h"
using namespace El;

extern "C" {

#define C_PROTO(SIG,SIGBASE,T) \
  ElError ElDistSparseMatrixCreate_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG * A, MPI_Comm comm ) \
  { EL_TRY( *A = CReflect( new DistSparseMatrix<T>(mpi::Comm(comm)) ) ) } \
  ElError ElDistSparseMatrixDestroy_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A ) \
  { EL_TRY( delete CReflect(A) ) } \
  ElError ElDistSparseMatrixEmpty_ ## SIG ( ElDistSparseMatrix_ ## SIG A ) \
  { EL_TRY( CReflect(A)->Empty() ) } \
  ElError ElDistSparseMatrixResize_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, ElInt height, ElInt width ) \
  { EL_TRY( CReflect(A)->Resize(height,width) ) } \
  ElError ElDistSparseMatrixSetComm_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, MPI_Comm comm ) \
  { EL_TRY( CReflect(A)->SetComm(mpi::Comm(comm)) ) } \
  ElError ElDistSparseMatrixReserve_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, \
    ElInt numLocalEntries, ElInt numRemoteEntries ) \
  { EL_TRY( CReflect(A)->Reserve(numLocalEntries,numRemoteEntries) ) } \
  ElError ElDistSparseMatrixUpdate_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, \
    ElInt row, ElInt col, CREFLECT(T) value ) \
  { EL_TRY( CReflect(A)->Update(row,col,CReflect(value)) ) } \
  ElError ElDistSparseMatrixUpdateLocal_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, \
    ElInt localRow, ElInt col, CREFLECT(T) value ) \
  { EL_TRY( CReflect(A)->UpdateLocal(localRow,col,CReflect(value)) ) } \
  ElError ElDistSparseMatrixZero_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, ElInt row, ElInt col ) \
  { EL_TRY( CReflect(A)->Zero(row,col) ) } \
  ElError ElDistSparseMatrixZeroLocal_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, ElInt localRow, ElInt col ) \
  { EL_TRY( CReflect(A)->ZeroLocal(localRow,col) ) } \
  ElError ElDistSparseMatrixQueueUpdate_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, \
    ElInt row, ElInt col, CREFLECT(T) value, bool passive ) \
  { EL_TRY( CReflect(A)->QueueUpdate(row,col,CReflect(value),passive) ) } \
  ElError ElDistSparseMatrixQueueLocalUpdate_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, \
    ElInt localRow, ElInt col, CREFLECT(T) value ) \
  { EL_TRY( CReflect(A)->QueueLocalUpdate(localRow,col,CReflect(value)) ) } \
  ElError ElDistSparseMatrixQueueZero_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, ElInt row, ElInt col, bool passive ) \
  { EL_TRY( CReflect(A)->QueueZero(row,col,passive) ) } \
  ElError ElDistSparseMatrixQueueLocalZero_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, ElInt localRow, ElInt col ) \
  { EL_TRY( CReflect(A)->QueueLocalZero(localRow,col) ) } \
  ElError ElDistSparseMatrixProcessQueues_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A ) \
  { EL_TRY( CReflect(A)->ProcessQueues() ) } \
  ElError ElDistSparseMatrixProcessLocalQueues_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A ) \
  { EL_TRY( CReflect(A)->ProcessLocalQueues() ) } \
  ElError ElDistSparseMatrixHeight_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt* height ) \
  { EL_TRY( *height = CReflect(A)->Height() ) } \
  ElError ElDistSparseMatrixWidth_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt* width ) \
  { EL_TRY( *width = CReflect(A)->Width() ) } \
  ElError ElDistSparseMatrixDistGraph_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, ElDistGraph* graph) \
  { EL_TRY( *graph = CReflect(&CReflect(A)->DistGraph()) ) } \
  ElError ElDistSparseMatrixLockedDistGraph_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistGraph* graph) \
  { EL_TRY( *graph = CReflect(&CReflect(A)->LockedDistGraph()) ) } \
  ElError ElDistSparseMatrixFirstLocalRow_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt* firstLocalRow ) \
  { EL_TRY( *firstLocalRow = CReflect(A)->FirstLocalRow() ) } \
  ElError ElDistSparseMatrixLocalHeight_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt* localHeight ) \
  { EL_TRY( *localHeight = CReflect(A)->LocalHeight() ) } \
  ElError ElDistSparseMatrixNumLocalEntries_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt* numLocalEntries ) \
  { EL_TRY( *numLocalEntries = CReflect(A)->NumLocalEntries() ) } \
  ElError ElDistSparseMatrixCapacity_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt* capacity ) \
  { EL_TRY( *capacity = CReflect(A)->Capacity() ) } \
  ElError ElDistSparseMatrixLocallyConsistent_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, bool* consistent ) \
  { EL_TRY( *consistent = CReflect(A)->LocallyConsistent() ) } \
  ElError ElDistSparseMatrixComm_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, MPI_Comm* comm ) \
  { EL_TRY( *comm = CReflect(A)->Comm().comm ) } \
  ElError ElDistSparseMatrixBlocksize_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt* blocksize ) \
  { EL_TRY( *blocksize = CReflect(A)->Blocksize() ) } \
  ElError ElDistSparseMatrixRowOwner_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt i, int* owner ) \
  { EL_TRY( *owner = CReflect(A)->RowOwner(i) ) } \
  ElError ElDistSparseMatrixGlobalRow_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt iLoc, ElInt* i ) \
  { EL_TRY( *i = CReflect(A)->GlobalRow(iLoc) ) } \
  ElError ElDistSparseMatrixRow_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt localInd, ElInt* row ) \
  { EL_TRY( *row = CReflect(A)->Row(localInd) ) } \
  ElError ElDistSparseMatrixCol_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt localInd, ElInt* col ) \
  { EL_TRY( *col = CReflect(A)->Col(localInd) ) } \
  ElError ElDistSparseMatrixValue_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt localInd, CREFLECT(T)* value ) \
  { EL_TRY( *value = CReflect(CReflect(A)->Value(localInd)) ) } \
  ElError ElDistSparseMatrixRowOffset_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt localRow, \
    ElInt* localRowOffset ) \
  { EL_TRY( *localRowOffset = CReflect(A)->RowOffset(localRow) ) } \
  ElError ElDistSparseMatrixOffset_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt localRow, ElInt col, \
    ElInt* localOffset ) \
  { EL_TRY( *localOffset = CReflect(A)->Offset(localRow,col) ) } \
  ElError ElDistSparseMatrixNumConnections_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt localRow, ElInt* numConnections ) \
  { EL_TRY( *numConnections = CReflect(A)->NumConnections(localRow) ) } \
  ElError ElDistSparseMatrixImbalance_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, double* imbalance ) \
  { EL_TRY( *imbalance = CReflect(A)->Imbalance() ) } \
  ElError ElDistSparseMatrixSourceBuffer_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, ElInt** sourceBuffer ) \
  { EL_TRY( *sourceBuffer = CReflect(A)->SourceBuffer() ) } \
  ElError ElDistSparseMatrixLockedSourceBuffer_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, const ElInt** sourceBuffer ) \
  { EL_TRY( *sourceBuffer = CReflect(A)->LockedSourceBuffer() ) } \
  ElError ElDistSparseMatrixTargetBuffer_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, ElInt** targetBuffer ) \
  { EL_TRY( *targetBuffer = CReflect(A)->TargetBuffer() ) } \
  ElError ElDistSparseMatrixLockedTargetBuffer_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, const ElInt** targetBuffer ) \
  { EL_TRY( *targetBuffer = CReflect(A)->LockedTargetBuffer() ) } \
  ElError ElDistSparseMatrixValueBuffer_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, CREFLECT(T)** valueBuffer ) \
  { EL_TRY( *valueBuffer = CReflect(CReflect(A)->ValueBuffer()) ) } \
  ElError ElDistSparseMatrixLockedValueBuffer_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, const CREFLECT(T)** valueBuffer ) \
  { EL_TRY( *valueBuffer = CReflect(CReflect(A)->LockedValueBuffer()) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
