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
  ElError ElSparseMatrixCreate_ ## SIG ( ElSparseMatrix_ ## SIG * A ) \
  { EL_TRY( *A = CReflect( new SparseMatrix<T>() ) ) } \
  ElError ElSparseMatrixDestroy_ ## SIG ( ElConstSparseMatrix_ ## SIG A ) \
  { EL_TRY( delete CReflect(A) ) } \
  ElError ElSparseMatrixEmpty_ ## SIG ( ElSparseMatrix_ ## SIG A ) \
  { EL_TRY( CReflect(A)->Empty() ) } \
  ElError ElSparseMatrixResize_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, ElInt height, ElInt width ) \
  { EL_TRY( CReflect(A)->Resize(height,width) ) } \
  ElError ElSparseMatrixReserve_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, ElInt numEntries ) \
  { EL_TRY( CReflect(A)->Reserve(numEntries) ) } \
  ElError ElSparseMatrixUpdate_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, ElInt row, ElInt col, CREFLECT(T) value ) \
  { EL_TRY( CReflect(A)->Update(row,col,CReflect(value)) ) } \
  ElError ElSparseMatrixZero_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, ElInt row, ElInt col ) \
  { EL_TRY( CReflect(A)->Zero(row,col) ) } \
  ElError ElSparseMatrixQueueUpdate_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, ElInt row, ElInt col, CREFLECT(T) value ) \
  { EL_TRY( CReflect(A)->QueueUpdate(row,col,CReflect(value)) ) } \
  ElError ElSparseMatrixQueueZero_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, ElInt row, ElInt col ) \
  { EL_TRY( CReflect(A)->QueueZero(row,col) ) } \
  ElError ElSparseMatrixProcessQueues_ ## SIG ( ElSparseMatrix_ ## SIG A ) \
  { EL_TRY( CReflect(A)->ProcessQueues() ) } \
  ElError ElSparseMatrixHeight_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElInt* height ) \
  { EL_TRY( *height = CReflect(A)->Height() ) } \
  ElError ElSparseMatrixWidth_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElInt* width ) \
  { EL_TRY( *width = CReflect(A)->Width() ) } \
  ElError ElSparseMatrixNumEntries_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElInt* numEntries ) \
  { EL_TRY( *numEntries = CReflect(A)->NumEntries() ) } \
  ElError ElSparseMatrixCapacity_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElInt* capacity ) \
  { EL_TRY( *capacity = CReflect(A)->Capacity() ) } \
  ElError ElSparseMatrixConsistent_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, bool* consistent ) \
  { EL_TRY( *consistent = CReflect(A)->Consistent() ) } \
  ElError ElSparseMatrixGraph_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, ElGraph* graph) \
  { EL_TRY( *graph = CReflect(&CReflect(A)->Graph()) ) } \
  ElError ElSparseMatrixLockedGraph_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstGraph* graph) \
  { EL_TRY( *graph = CReflect(&CReflect(A)->LockedGraph()) ) } \
  ElError ElSparseMatrixRow_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElInt index, ElInt* row ) \
  { EL_TRY( *row = CReflect(A)->Row(index) ) } \
  ElError ElSparseMatrixCol_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElInt index, ElInt* col ) \
  { EL_TRY( *col = CReflect(A)->Col(index) ) } \
  ElError ElSparseMatrixValue_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElInt index, CREFLECT(T)* value ) \
  { EL_TRY( *value = CReflect(CReflect(A)->Value(index)) ) } \
  ElError ElSparseMatrixRowOffset_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElInt row, ElInt* rowOffset ) \
  { EL_TRY( *rowOffset = CReflect(A)->RowOffset(row) ) } \
  ElError ElSparseMatrixOffset_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElInt row, ElInt col, ElInt* offset ) \
  { EL_TRY( *offset = CReflect(A)->Offset(row,col) ) } \
  ElError ElSparseMatrixNumConnections_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElInt row, ElInt* numConnections ) \
  { EL_TRY( *numConnections = CReflect(A)->NumConnections(row) ) } \
  ElError ElSparseMatrixSourceBuffer_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, ElInt** sourceBuffer ) \
  { EL_TRY( *sourceBuffer = CReflect(A)->SourceBuffer() ) } \
  ElError ElSparseMatrixLockedSourceBuffer_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, const ElInt** sourceBuffer ) \
  { EL_TRY( *sourceBuffer = CReflect(A)->LockedSourceBuffer() ) } \
  ElError ElSparseMatrixTargetBuffer_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, ElInt** targetBuffer ) \
  { EL_TRY( *targetBuffer = CReflect(A)->TargetBuffer() ) } \
  ElError ElSparseMatrixLockedTargetBuffer_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, const ElInt** targetBuffer ) \
  { EL_TRY( *targetBuffer = CReflect(A)->LockedTargetBuffer() ) } \
  ElError ElSparseMatrixValueBuffer_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, CREFLECT(T)** valueBuffer ) \
  { EL_TRY( *valueBuffer = CReflect(CReflect(A)->ValueBuffer()) ) } \
  ElError ElSparseMatrixLockedValueBuffer_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, const CREFLECT(T)** valueBuffer ) \
  { EL_TRY( *valueBuffer = CReflect(CReflect(A)->LockedValueBuffer()) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
