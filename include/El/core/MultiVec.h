/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_MULTIVEC_C_H
#define EL_MULTIVEC_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* An anonymous struct meant as a placeholder for MultiVec<T>
   -------------------------------------------------------------- */
typedef struct ElMultiVec_iDummy* ElMultiVec_i;
typedef struct ElMultiVec_sDummy* ElMultiVec_s;
typedef struct ElMultiVec_dDummy* ElMultiVec_d;
typedef struct ElMultiVec_cDummy* ElMultiVec_c;
typedef struct ElMultiVec_zDummy* ElMultiVec_z;

typedef const struct ElMultiVec_iDummy* ElConstMultiVec_i;
typedef const struct ElMultiVec_sDummy* ElConstMultiVec_s;
typedef const struct ElMultiVec_dDummy* ElConstMultiVec_d;
typedef const struct ElMultiVec_cDummy* ElConstMultiVec_c;
typedef const struct ElMultiVec_zDummy* ElConstMultiVec_z;

/* Constructors and destructors
   ============================ */

/* MultiVec<T>::MultiVec()
   ----------------------- */
EL_EXPORT ElError ElMultiVecCreate_i( ElMultiVec_i* A );
EL_EXPORT ElError ElMultiVecCreate_s( ElMultiVec_s* A );
EL_EXPORT ElError ElMultiVecCreate_d( ElMultiVec_d* A );
EL_EXPORT ElError ElMultiVecCreate_c( ElMultiVec_c* A );
EL_EXPORT ElError ElMultiVecCreate_z( ElMultiVec_z* A );

/* MultiVec<T>::~MultiVec()
   ------------------------ */
EL_EXPORT ElError ElMultiVecDestroy_i( ElConstMultiVec_i A );
EL_EXPORT ElError ElMultiVecDestroy_s( ElConstMultiVec_s A );
EL_EXPORT ElError ElMultiVecDestroy_d( ElConstMultiVec_d A );
EL_EXPORT ElError ElMultiVecDestroy_c( ElConstMultiVec_c A );
EL_EXPORT ElError ElMultiVecDestroy_z( ElConstMultiVec_z A );

/* Assignment and reconfiguration
   ============================== */

/* void MultiVec<T>::Empty()
   ------------------------- */
EL_EXPORT ElError ElMultiVecEmpty_i( ElMultiVec_i A );
EL_EXPORT ElError ElMultiVecEmpty_s( ElMultiVec_s A );
EL_EXPORT ElError ElMultiVecEmpty_d( ElMultiVec_d A );
EL_EXPORT ElError ElMultiVecEmpty_c( ElMultiVec_c A );
EL_EXPORT ElError ElMultiVecEmpty_z( ElMultiVec_z A );

/* void MultiVec<T>::Resize( Int height, Int width )
   ------------------------------------------------- */
EL_EXPORT ElError ElMultiVecResize_i
( ElMultiVec_i A, ElInt height, ElInt width );
EL_EXPORT ElError ElMultiVecResize_s
( ElMultiVec_s A, ElInt height, ElInt width );
EL_EXPORT ElError ElMultiVecResize_d
( ElMultiVec_d A, ElInt height, ElInt width );
EL_EXPORT ElError ElMultiVecResize_c
( ElMultiVec_c A, ElInt height, ElInt width );
EL_EXPORT ElError ElMultiVecResize_z
( ElMultiVec_z A, ElInt height, ElInt width );

/* Queries
   ======= */

/* Int MultiVec<T>::Height() const
   ------------------------------- */
EL_EXPORT ElError ElMultiVecHeight_i( ElConstMultiVec_i A, ElInt* height );
EL_EXPORT ElError ElMultiVecHeight_s( ElConstMultiVec_s A, ElInt* height );
EL_EXPORT ElError ElMultiVecHeight_d( ElConstMultiVec_d A, ElInt* height );
EL_EXPORT ElError ElMultiVecHeight_c( ElConstMultiVec_c A, ElInt* height );
EL_EXPORT ElError ElMultiVecHeight_z( ElConstMultiVec_z A, ElInt* height );

/* Int MultiVec<T>::Width() const
   ------------------------------ */
EL_EXPORT ElError ElMultiVecWidth_i( ElConstMultiVec_i A, ElInt* width );
EL_EXPORT ElError ElMultiVecWidth_s( ElConstMultiVec_s A, ElInt* width );
EL_EXPORT ElError ElMultiVecWidth_d( ElConstMultiVec_d A, ElInt* width );
EL_EXPORT ElError ElMultiVecWidth_c( ElConstMultiVec_c A, ElInt* width );
EL_EXPORT ElError ElMultiVecWidth_z( ElConstMultiVec_z A, ElInt* width );

/* Entrywise manipulation
   ====================== */

/* T MultiVec<T>::Get( Int i, Int j ) const
   ---------------------------------------- */
EL_EXPORT ElError ElMultiVecGet_i
( ElConstMultiVec_i A, ElInt i, ElInt j, ElInt* value );
EL_EXPORT ElError ElMultiVecGet_s
( ElConstMultiVec_s A, ElInt i, ElInt j, float* value );
EL_EXPORT ElError ElMultiVecGet_d
( ElConstMultiVec_d A, ElInt i, ElInt j, double* value );
EL_EXPORT ElError ElMultiVecGet_c
( ElConstMultiVec_c A, ElInt i, ElInt j, complex_float* value );
EL_EXPORT ElError ElMultiVecGet_z
( ElConstMultiVec_z A, ElInt i, ElInt j, complex_double* value );

/* void MultiVec<T>::Set( Int i, Int j, T value )
   ---------------------------------------------- */
EL_EXPORT ElError ElMultiVecSet_i
( ElMultiVec_i A, ElInt i, ElInt j, ElInt value );
EL_EXPORT ElError ElMultiVecSet_s
( ElMultiVec_s A, ElInt i, ElInt j, float value );
EL_EXPORT ElError ElMultiVecSet_d
( ElMultiVec_d A, ElInt i, ElInt j, double value );
EL_EXPORT ElError ElMultiVecSet_c
( ElMultiVec_c A, ElInt i, ElInt j, complex_float value );
EL_EXPORT ElError ElMultiVecSet_z
( ElMultiVec_z A, ElInt i, ElInt j, complex_double value );

/* void MultiVec<T>::Update( Int i, Int j, T value )
   ------------------------------------------------- */
EL_EXPORT ElError ElMultiVecUpdate_i
( ElMultiVec_i A, ElInt i, ElInt j, ElInt value );
EL_EXPORT ElError ElMultiVecUpdate_s
( ElMultiVec_s A, ElInt i, ElInt j, float value );
EL_EXPORT ElError ElMultiVecUpdate_d
( ElMultiVec_d A, ElInt i, ElInt j, double value );
EL_EXPORT ElError ElMultiVecUpdate_c
( ElMultiVec_c A, ElInt i, ElInt j, complex_float value );
EL_EXPORT ElError ElMultiVecUpdate_z
( ElMultiVec_z A, ElInt i, ElInt j, complex_double value );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_MULTIVEC_C_H */
