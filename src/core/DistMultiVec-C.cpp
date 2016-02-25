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
  /* Constructors and destructors */ \
  ElError ElDistMultiVecCreate_ ## SIG \
  ( ElDistMultiVec_ ## SIG * A, MPI_Comm comm ) \
  { EL_TRY( *A = CReflect( new DistMultiVec<T>(mpi::Comm(comm)) ) ) } \
  ElError ElDistMultiVecDestroy_ ## SIG ( ElConstDistMultiVec_ ## SIG A ) \
  { EL_TRY( delete CReflect(A) ) } \
  /* Assignment and reconfiguration */ \
  ElError ElDistMultiVecEmpty_ ## SIG ( ElDistMultiVec_ ## SIG A ) \
  { EL_TRY( CReflect(A)->Empty() ) } \
  ElError ElDistMultiVecResize_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, ElInt height, ElInt width ) \
  { EL_TRY( CReflect(A)->Resize(height,width) ) } \
  ElError ElDistMultiVecSetComm_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, MPI_Comm comm ) \
  { EL_TRY( CReflect(A)->SetComm(mpi::Comm(comm)) ) } \
  /* Queries */ \
  ElError ElDistMultiVecHeight_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElInt* height ) \
  { EL_TRY( *height = CReflect(A)->Height() ) } \
  ElError ElDistMultiVecWidth_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElInt* width ) \
  { EL_TRY( *width = CReflect(A)->Width() ) } \
  ElError ElDistMultiVecFirstLocalRow_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElInt* firstLocalRow ) \
  { EL_TRY( *firstLocalRow = CReflect(A)->FirstLocalRow() ) } \
  ElError ElDistMultiVecLocalHeight_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElInt* localHeight ) \
  { EL_TRY( *localHeight = CReflect(A)->LocalHeight() ) } \
  ElError ElDistMultiVecMatrix_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, ElMatrix_ ## SIG * ALoc ) \
  { EL_TRY( *ALoc = CReflect(&CReflect(A)->Matrix()) ) } \
  ElError ElDistMultiVecLockedMatrix_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElConstMatrix_ ## SIG * ALoc ) \
  { EL_TRY( *ALoc = CReflect(&CReflect(A)->LockedMatrix()) ) } \
  ElError ElDistMultiVecComm_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, MPI_Comm* comm ) \
  { EL_TRY( *comm = CReflect(A)->Comm().comm ) } \
  ElError ElDistMultiVecBlocksize_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElInt* blocksize ) \
  { EL_TRY( *blocksize = CReflect(A)->Blocksize() ) } \
  ElError ElDistMultiVecRowOwner_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElInt i, int* owner ) \
  { EL_TRY( *owner = CReflect(A)->RowOwner(i) ) } \
  ElError ElDistMultiVecGlobalRow_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElInt iLoc, ElInt* i ) \
  { EL_TRY( *i = CReflect(A)->GlobalRow(iLoc) ) } \
  /* Entrywise manipulation */ \
  ElError ElDistMultiVecGet_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElInt i, ElInt j, CREFLECT(T)* value ) \
  { EL_TRY( *value = CReflect(CReflect(A)->Get(i,j)) ) } \
  ElError ElDistMultiVecSet_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, ElInt i, ElInt j, CREFLECT(T) value ) \
  { EL_TRY( CReflect(A)->Set(i,j,CReflect(value)) ) } \
  ElError ElDistMultiVecGetLocal_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElInt iLocal, ElInt j, CREFLECT(T)* value ) \
  { EL_TRY( *value = CReflect(CReflect(A)->GetLocal(iLocal,j)) ) } \
  ElError ElDistMultiVecSetLocal_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, ElInt iLocal, ElInt j, CREFLECT(T) value ) \
  { EL_TRY( CReflect(A)->SetLocal(iLocal,j,CReflect(value)) ) } \
  ElError ElDistMultiVecUpdateLocal_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, ElInt iLocal, ElInt j, CREFLECT(T) value ) \
  { EL_TRY( CReflect(A)->UpdateLocal(iLocal,j,CReflect(value)) ) } \
  /* void Reserve( Int numRemoteEntries ) */ \
  ElError ElDistMultiVecReserve_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, ElInt numRemoteEntries ) \
  { EL_TRY( CReflect(A)->Reserve(numRemoteEntries) ) } \
  /* void QueueUpdate( Int i, Int j, T value ) */ \
  ElError ElDistMultiVecQueueUpdate_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, ElInt i, ElInt j, CREFLECT(T) value ) \
  { EL_TRY( CReflect(A)->QueueUpdate(i,j,CReflect(value)) ) } \
  /* void ProcessQueues() */ \
  ElError ElDistMultiVecProcessQueues_ ## SIG( ElDistMultiVec_ ## SIG A ) \
  { EL_TRY( CReflect(A)->ProcessQueues() ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
