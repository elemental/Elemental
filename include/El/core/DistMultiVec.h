/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_DISTMULTIVEC_C_H
#define EL_DISTMULTIVEC_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* An anonymous struct meant as a placeholder for DistMultiVec<T>
   -------------------------------------------------------------- */
typedef struct ElDistMultiVec_iDummy* ElDistMultiVec_i;
typedef struct ElDistMultiVec_sDummy* ElDistMultiVec_s;
typedef struct ElDistMultiVec_dDummy* ElDistMultiVec_d;
typedef struct ElDistMultiVec_cDummy* ElDistMultiVec_c;
typedef struct ElDistMultiVec_zDummy* ElDistMultiVec_z;

typedef const struct ElDistMultiVec_iDummy* ElConstDistMultiVec_i;
typedef const struct ElDistMultiVec_sDummy* ElConstDistMultiVec_s;
typedef const struct ElDistMultiVec_dDummy* ElConstDistMultiVec_d;
typedef const struct ElDistMultiVec_cDummy* ElConstDistMultiVec_c;
typedef const struct ElDistMultiVec_zDummy* ElConstDistMultiVec_z;

/* Constructors and destructors
   ============================ */

/* DistMultiVec<T>::DistMultiVec( const Grid& grid )
   ------------------------------------------------- */
EL_EXPORT ElError ElDistMultiVecCreate_i
( ElDistMultiVec_i* A, ElConstGrid grid );
EL_EXPORT ElError ElDistMultiVecCreate_s
( ElDistMultiVec_s* A, ElConstGrid grid );
EL_EXPORT ElError ElDistMultiVecCreate_d
( ElDistMultiVec_d* A, ElConstGrid grid );
EL_EXPORT ElError ElDistMultiVecCreate_c
( ElDistMultiVec_c* A, ElConstGrid grid );
EL_EXPORT ElError ElDistMultiVecCreate_z
( ElDistMultiVec_z* A, ElConstGrid grid );

/* DistMultiVec<T>::~DistMultiVec()
   -------------------------------- */
EL_EXPORT ElError ElDistMultiVecDestroy_i( ElConstDistMultiVec_i A );
EL_EXPORT ElError ElDistMultiVecDestroy_s( ElConstDistMultiVec_s A );
EL_EXPORT ElError ElDistMultiVecDestroy_d( ElConstDistMultiVec_d A );
EL_EXPORT ElError ElDistMultiVecDestroy_c( ElConstDistMultiVec_c A );
EL_EXPORT ElError ElDistMultiVecDestroy_z( ElConstDistMultiVec_z A );

/* Assignment and reconfiguration
   ============================== */

/* void DistMultiVec<T>::Empty()
   ----------------------------- */
EL_EXPORT ElError ElDistMultiVecEmpty_i( ElDistMultiVec_i A );
EL_EXPORT ElError ElDistMultiVecEmpty_s( ElDistMultiVec_s A );
EL_EXPORT ElError ElDistMultiVecEmpty_d( ElDistMultiVec_d A );
EL_EXPORT ElError ElDistMultiVecEmpty_c( ElDistMultiVec_c A );
EL_EXPORT ElError ElDistMultiVecEmpty_z( ElDistMultiVec_z A );

/* void DistMultiVec<T>::Resize( Int height, Int width )
   ----------------------------------------------------- */
EL_EXPORT ElError ElDistMultiVecResize_i
( ElDistMultiVec_i A, ElInt height, ElInt width );
EL_EXPORT ElError ElDistMultiVecResize_s
( ElDistMultiVec_s A, ElInt height, ElInt width );
EL_EXPORT ElError ElDistMultiVecResize_d
( ElDistMultiVec_d A, ElInt height, ElInt width );
EL_EXPORT ElError ElDistMultiVecResize_c
( ElDistMultiVec_c A, ElInt height, ElInt width );
EL_EXPORT ElError ElDistMultiVecResize_z
( ElDistMultiVec_z A, ElInt height, ElInt width );

/* void DistMultiVec<T>::SetGrid( const Grid& grid )
   ------------------------------------------------- */
EL_EXPORT ElError ElDistMultiVecSetGrid_i
( ElDistMultiVec_i A, ElConstGrid grid );
EL_EXPORT ElError ElDistMultiVecSetGrid_s
( ElDistMultiVec_s A, ElConstGrid grid );
EL_EXPORT ElError ElDistMultiVecSetGrid_d
( ElDistMultiVec_d A, ElConstGrid grid );
EL_EXPORT ElError ElDistMultiVecSetGrid_c
( ElDistMultiVec_c A, ElConstGrid grid );
EL_EXPORT ElError ElDistMultiVecSetGrid_z
( ElDistMultiVec_z A, ElConstGrid grid );

/* Queries
   ======= */

/* Int DistMultiVec<T>::Height() const
   ----------------------------------- */
EL_EXPORT ElError ElDistMultiVecHeight_i
( ElConstDistMultiVec_i A, ElInt* height );
EL_EXPORT ElError ElDistMultiVecHeight_s
( ElConstDistMultiVec_s A, ElInt* height );
EL_EXPORT ElError ElDistMultiVecHeight_d
( ElConstDistMultiVec_d A, ElInt* height );
EL_EXPORT ElError ElDistMultiVecHeight_c
( ElConstDistMultiVec_c A, ElInt* height );
EL_EXPORT ElError ElDistMultiVecHeight_z
( ElConstDistMultiVec_z A, ElInt* height );

/* Int DistMultiVec<T>::Width() const
   ---------------------------------- */
EL_EXPORT ElError ElDistMultiVecWidth_i
( ElConstDistMultiVec_i A, ElInt* width );
EL_EXPORT ElError ElDistMultiVecWidth_s
( ElConstDistMultiVec_s A, ElInt* width );
EL_EXPORT ElError ElDistMultiVecWidth_d
( ElConstDistMultiVec_d A, ElInt* width );
EL_EXPORT ElError ElDistMultiVecWidth_c
( ElConstDistMultiVec_c A, ElInt* width );
EL_EXPORT ElError ElDistMultiVecWidth_z
( ElConstDistMultiVec_z A, ElInt* width );

/* Int DistMultiVec<T>::FirstLocalRow() const
   ------------------------------------------ */
EL_EXPORT ElError ElDistMultiVecFirstLocalRow_i
( ElConstDistMultiVec_i A, ElInt* firstLocalRow );
EL_EXPORT ElError ElDistMultiVecFirstLocalRow_s
( ElConstDistMultiVec_s A, ElInt* firstLocalRow );
EL_EXPORT ElError ElDistMultiVecFirstLocalRow_d
( ElConstDistMultiVec_d A, ElInt* firstLocalRow );
EL_EXPORT ElError ElDistMultiVecFirstLocalRow_c
( ElConstDistMultiVec_c A, ElInt* firstLocalRow );
EL_EXPORT ElError ElDistMultiVecFirstLocalRow_z
( ElConstDistMultiVec_z A, ElInt* firstLocalRow );

/* Int DistMultiVec<T>::LocalHeight() const
   ---------------------------------------- */
EL_EXPORT ElError ElDistMultiVecLocalHeight_i
( ElConstDistMultiVec_i A, ElInt* localHeight );
EL_EXPORT ElError ElDistMultiVecLocalHeight_s
( ElConstDistMultiVec_s A, ElInt* localHeight );
EL_EXPORT ElError ElDistMultiVecLocalHeight_d
( ElConstDistMultiVec_d A, ElInt* localHeight );
EL_EXPORT ElError ElDistMultiVecLocalHeight_c
( ElConstDistMultiVec_c A, ElInt* localHeight );
EL_EXPORT ElError ElDistMultiVecLocalHeight_z
( ElConstDistMultiVec_z A, ElInt* localHeight );

/* Matrix<T>& DistMultiVec<T>::Matrix()
   ------------------------------------ */
EL_EXPORT ElError ElDistMultiVecMatrix_i
( ElDistMultiVec_i A, ElMatrix_i* ALoc );
EL_EXPORT ElError ElDistMultiVecMatrix_s
( ElDistMultiVec_s A, ElMatrix_s* ALoc );
EL_EXPORT ElError ElDistMultiVecMatrix_d
( ElDistMultiVec_d A, ElMatrix_d* ALoc );
EL_EXPORT ElError ElDistMultiVecMatrix_c
( ElDistMultiVec_c A, ElMatrix_c* ALoc );
EL_EXPORT ElError ElDistMultiVecMatrix_z
( ElDistMultiVec_z A, ElMatrix_z* ALoc );

/* const Matrix<T>& DistMultiVec<T>::LockedMatrix() const
   ------------------------------------------------------ */
EL_EXPORT ElError ElDistMultiVecLockedMatrix_i
( ElConstDistMultiVec_i A, ElConstMatrix_i* ALoc );
EL_EXPORT ElError ElDistMultiVecLockedMatrix_s
( ElConstDistMultiVec_s A, ElConstMatrix_s* ALoc );
EL_EXPORT ElError ElDistMultiVecLockedMatrix_d
( ElConstDistMultiVec_d A, ElConstMatrix_d* ALoc );
EL_EXPORT ElError ElDistMultiVecLockedMatrix_c
( ElConstDistMultiVec_c A, ElConstMatrix_c* ALoc );
EL_EXPORT ElError ElDistMultiVecLockedMatrix_z
( ElConstDistMultiVec_z A, ElConstMatrix_z* ALoc );

/* const Grid& DistMultiVec<T>::Grid() const
   ----------------------------------------- */
EL_EXPORT ElError ElDistMultiVecGrid_i
( ElConstDistMultiVec_i A, ElConstGrid* grid );
EL_EXPORT ElError ElDistMultiVecGrid_s
( ElConstDistMultiVec_s A, ElConstGrid* grid );
EL_EXPORT ElError ElDistMultiVecGrid_d
( ElConstDistMultiVec_d A, ElConstGrid* grid );
EL_EXPORT ElError ElDistMultiVecGrid_c
( ElConstDistMultiVec_c A, ElConstGrid* grid );
EL_EXPORT ElError ElDistMultiVecGrid_z
( ElConstDistMultiVec_z A, ElConstGrid* grid );

/* Int DistMultiVec<T>::Blocksize() const
   -------------------------------------- */
EL_EXPORT ElError ElDistMultiVecBlocksize_i
( ElConstDistMultiVec_i A, ElInt* blocksize );
EL_EXPORT ElError ElDistMultiVecBlocksize_s
( ElConstDistMultiVec_s A, ElInt* blocksize );
EL_EXPORT ElError ElDistMultiVecBlocksize_d
( ElConstDistMultiVec_d A, ElInt* blocksize );
EL_EXPORT ElError ElDistMultiVecBlocksize_c
( ElConstDistMultiVec_c A, ElInt* blocksize );
EL_EXPORT ElError ElDistMultiVecBlocksize_z
( ElConstDistMultiVec_z A, ElInt* blocksize );

/* int DistMultiVec<T>::RowOwner( Int i ) const
   -------------------------------------------- */
EL_EXPORT ElError ElDistMultiVecRowOwner_i
( ElConstDistMultiVec_i A, ElInt i, int* owner );
EL_EXPORT ElError ElDistMultiVecRowOwner_s
( ElConstDistMultiVec_s A, ElInt i, int* owner );
EL_EXPORT ElError ElDistMultiVecRowOwner_d
( ElConstDistMultiVec_d A, ElInt i, int* owner );
EL_EXPORT ElError ElDistMultiVecRowOwner_c
( ElConstDistMultiVec_c A, ElInt i, int* owner );
EL_EXPORT ElError ElDistMultiVecRowOwner_z
( ElConstDistMultiVec_z A, ElInt i, int* owner );

/* Int DistMultiVec<T>::GlobalRow( Int iLoc ) const
   ------------------------------------------------ */
EL_EXPORT ElError ElDistMultiVecGlobalRow_i
( ElConstDistMultiVec_i A, ElInt iLoc, ElInt* i );
EL_EXPORT ElError ElDistMultiVecGlobalRow_s
( ElConstDistMultiVec_s A, ElInt iLoc, ElInt* i );
EL_EXPORT ElError ElDistMultiVecGlobalRow_d
( ElConstDistMultiVec_d A, ElInt iLoc, ElInt* i );
EL_EXPORT ElError ElDistMultiVecGlobalRow_c
( ElConstDistMultiVec_c A, ElInt iLoc, ElInt* i );
EL_EXPORT ElError ElDistMultiVecGlobalRow_z
( ElConstDistMultiVec_z A, ElInt iLoc, ElInt* i );

/* Entrywise manipulation
   ====================== */
/* T DistMultiVec<T>::Get( Int i, Int j ) const
   -------------------------------------------- */
EL_EXPORT ElError ElDistMultiVecGet_i
( ElConstDistMultiVec_i A, ElInt i, ElInt j, ElInt* value );
EL_EXPORT ElError ElDistMultiVecGet_s
( ElConstDistMultiVec_s A, ElInt i, ElInt j, float* value );
EL_EXPORT ElError ElDistMultiVecGet_d
( ElConstDistMultiVec_d A, ElInt i, ElInt j, double* value );
EL_EXPORT ElError ElDistMultiVecGet_c
( ElConstDistMultiVec_c A, ElInt i, ElInt j, complex_float* value );
EL_EXPORT ElError ElDistMultiVecGet_z
( ElConstDistMultiVec_z A, ElInt i, ElInt j, complex_double* value );

/* void DistMultiVec<T>::Set( Int i, Int j, T value )
   -------------------------------------------------- */
EL_EXPORT ElError ElDistMultiVecSet_i
( ElDistMultiVec_i A, ElInt i, ElInt j, ElInt value );
EL_EXPORT ElError ElDistMultiVecSet_s
( ElDistMultiVec_s A, ElInt i, ElInt j, float value );
EL_EXPORT ElError ElDistMultiVecSet_d
( ElDistMultiVec_d A, ElInt i, ElInt j, double value );
EL_EXPORT ElError ElDistMultiVecSet_c
( ElDistMultiVec_c A, ElInt i, ElInt j, complex_float value );
EL_EXPORT ElError ElDistMultiVecSet_z
( ElDistMultiVec_z A, ElInt i, ElInt j, complex_double value );

/* void DistMultiVec<T>::Reserve( Int numRemoteEntries )
   ----------------------------------------------------- */
EL_EXPORT ElError ElDistMultiVecReserve_i
( ElDistMultiVec_i A, ElInt numRemoteEntries );
EL_EXPORT ElError ElDistMultiVecReserve_s
( ElDistMultiVec_s A, ElInt numRemoteEntries );
EL_EXPORT ElError ElDistMultiVecReserve_d
( ElDistMultiVec_d A, ElInt numRemoteEntries );
EL_EXPORT ElError ElDistMultiVecReserve_c
( ElDistMultiVec_c A, ElInt numRemoteEntries );
EL_EXPORT ElError ElDistMultiVecReserve_z
( ElDistMultiVec_z A, ElInt numRemoteEntries );

/* void DistMultiVec<T>::QueueUpdate( Int i, Int j, T value )
   ---------------------------------------------------------- */
EL_EXPORT ElError ElDistMultiVecQueueUpdate_i
( ElDistMultiVec_i A, ElInt i, ElInt j, ElInt value );
EL_EXPORT ElError ElDistMultiVecQueueUpdate_s
( ElDistMultiVec_s A, ElInt i, ElInt j, float value );
EL_EXPORT ElError ElDistMultiVecQueueUpdate_d
( ElDistMultiVec_d A, ElInt i, ElInt j, double value );
EL_EXPORT ElError ElDistMultiVecQueueUpdate_c
( ElDistMultiVec_c A, ElInt i, ElInt j, complex_float value );
EL_EXPORT ElError ElDistMultiVecQueueUpdate_z
( ElDistMultiVec_z A, ElInt i, ElInt j, complex_double value );

/* void DistMultiVec<T>::ProcessQueues()
   ------------------------------------- */
EL_EXPORT ElError ElDistMultiVecProcessQueues_i( ElDistMultiVec_i A );
EL_EXPORT ElError ElDistMultiVecProcessQueues_s( ElDistMultiVec_s A );
EL_EXPORT ElError ElDistMultiVecProcessQueues_d( ElDistMultiVec_d A );
EL_EXPORT ElError ElDistMultiVecProcessQueues_c( ElDistMultiVec_c A );
EL_EXPORT ElError ElDistMultiVecProcessQueues_z( ElDistMultiVec_z A );

/* T DistMultiVec<T>::GetLocal( Int iLocal, Int j ) const
   ------------------------------------------------------ */
EL_EXPORT ElError ElDistMultiVecGetLocal_i
( ElConstDistMultiVec_i A, ElInt iLocal, ElInt j, ElInt* value );
EL_EXPORT ElError ElDistMultiVecGetLocal_s
( ElConstDistMultiVec_s A, ElInt iLocal, ElInt j, float* value );
EL_EXPORT ElError ElDistMultiVecGetLocal_d
( ElConstDistMultiVec_d A, ElInt iLocal, ElInt j, double* value );
EL_EXPORT ElError ElDistMultiVecGetLocal_c
( ElConstDistMultiVec_c A, ElInt iLocal, ElInt j, complex_float* value );
EL_EXPORT ElError ElDistMultiVecGetLocal_z
( ElConstDistMultiVec_z A, ElInt iLocal, ElInt j, complex_double* value );

/* void DistMultiVec<T>::SetLocal( Int iLocal, Int j, T value )
   ------------------------------------------------------------ */
EL_EXPORT ElError ElDistMultiVecSetLocal_i
( ElDistMultiVec_i A, ElInt iLocal, ElInt j, ElInt value );
EL_EXPORT ElError ElDistMultiVecSetLocal_s
( ElDistMultiVec_s A, ElInt iLocal, ElInt j, float value );
EL_EXPORT ElError ElDistMultiVecSetLocal_d
( ElDistMultiVec_d A, ElInt iLocal, ElInt j, double value );
EL_EXPORT ElError ElDistMultiVecSetLocal_c
( ElDistMultiVec_c A, ElInt iLocal, ElInt j, complex_float value );
EL_EXPORT ElError ElDistMultiVecSetLocal_z
( ElDistMultiVec_z A, ElInt iLocal, ElInt j, complex_double value );

/* void DistMultiVec<T>::UpdateLocal( Int iLocal, Int j, T value )
   --------------------------------------------------------------- */
EL_EXPORT ElError ElDistMultiVecUpdateLocal_i
( ElDistMultiVec_i A, ElInt iLocal, ElInt j, ElInt value );
EL_EXPORT ElError ElDistMultiVecUpdateLocal_s
( ElDistMultiVec_s A, ElInt iLocal, ElInt j, float value );
EL_EXPORT ElError ElDistMultiVecUpdateLocal_d
( ElDistMultiVec_d A, ElInt iLocal, ElInt j, double value );
EL_EXPORT ElError ElDistMultiVecUpdateLocal_c
( ElDistMultiVec_c A, ElInt iLocal, ElInt j, complex_float value );
EL_EXPORT ElError ElDistMultiVecUpdateLocal_z
( ElDistMultiVec_z A, ElInt iLocal, ElInt j, complex_double value );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_DISTMULTIVEC_C_H */
