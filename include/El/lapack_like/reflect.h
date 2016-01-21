/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_REFLECT_C_H
#define EL_REFLECT_C_H

#include "El/core/DistMatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Apply packed reflectors
   ======================= */
EL_EXPORT ElError ElApplyPackedReflectors_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, ElInt offset, 
  ElConstMatrix_s H, ElConstMatrix_s t, ElMatrix_s A );
EL_EXPORT ElError ElApplyPackedReflectors_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, ElInt offset, 
  ElConstMatrix_d H, ElConstMatrix_d t, ElMatrix_d A );
EL_EXPORT ElError ElApplyPackedReflectors_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, 
  ElConjugation conjugation, ElInt offset, 
  ElConstMatrix_c H, ElConstMatrix_c t, ElMatrix_c A );
EL_EXPORT ElError ElApplyPackedReflectors_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, 
  ElConjugation conjugation, ElInt offset, 
  ElConstMatrix_z H, ElConstMatrix_z t, ElMatrix_z A );

EL_EXPORT ElError ElApplyPackedReflectorsDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, ElInt offset, 
  ElConstDistMatrix_s H, ElConstDistMatrix_s t, ElDistMatrix_s A );
EL_EXPORT ElError ElApplyPackedReflectorsDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, ElInt offset, 
  ElConstDistMatrix_d H, ElConstDistMatrix_d t, ElDistMatrix_d A );
EL_EXPORT ElError ElApplyPackedReflectorsDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, 
  ElConjugation conjugation, ElInt offset, 
  ElConstDistMatrix_c H, ElConstDistMatrix_c t, ElDistMatrix_c A );
EL_EXPORT ElError ElApplyPackedReflectorsDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, 
  ElConjugation conjugation, ElInt offset, 
  ElConstDistMatrix_z H, ElConstDistMatrix_z t, ElDistMatrix_z A );

/* Expand packed reflectors
   ======================== */
EL_EXPORT ElError ElExpandPackedReflectors_s
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElInt offset, 
  ElMatrix_s H, ElConstMatrix_s t );
EL_EXPORT ElError ElExpandPackedReflectors_d
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElInt offset, 
  ElMatrix_d H, ElConstMatrix_d t );
EL_EXPORT ElError ElExpandPackedReflectors_c
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElConjugation conjugation,
  ElInt offset, ElMatrix_c H, ElConstMatrix_c t );
EL_EXPORT ElError ElExpandPackedReflectors_z
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElConjugation conjugation,
  ElInt offset, ElMatrix_z H, ElConstMatrix_z t );

EL_EXPORT ElError ElExpandPackedReflectorsDist_s
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElInt offset, 
  ElDistMatrix_s H, ElConstDistMatrix_s t );
EL_EXPORT ElError ElExpandPackedReflectorsDist_d
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElInt offset, 
  ElDistMatrix_d H, ElConstDistMatrix_d t );
EL_EXPORT ElError ElExpandPackedReflectorsDist_c
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElConjugation conjugation,
  ElInt offset, ElDistMatrix_c H, ElConstDistMatrix_c t );
EL_EXPORT ElError ElExpandPackedReflectorsDist_z
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElConjugation conjugation,
  ElInt offset, ElDistMatrix_z H, ElConstDistMatrix_z t );

/* Hyperbolic reflector
   ==================== */

/* Left application
   ---------------- */
EL_EXPORT ElError ElLeftHyperbolicReflector_s
( float* chi, ElMatrix_s x, float* tau );
EL_EXPORT ElError ElLeftHyperbolicReflector_d
( double* chi, ElMatrix_d x, double* tau );
EL_EXPORT ElError ElLeftHyperbolicReflector_c
( complex_float* chi, ElMatrix_c x, complex_float* tau );
EL_EXPORT ElError ElLeftHyperbolicReflector_z
( complex_double* chi, ElMatrix_z x, complex_double* tau );

EL_EXPORT ElError ElLeftHyperbolicReflectorDist_s
( float* chi, ElDistMatrix_s x, float* tau );
EL_EXPORT ElError ElLeftHyperbolicReflectorDist_d
( double* chi, ElDistMatrix_d x, double* tau );
EL_EXPORT ElError ElLeftHyperbolicReflectorDist_c
( complex_float* chi, ElDistMatrix_c x, complex_float* tau );
EL_EXPORT ElError ElLeftHyperbolicReflectorDist_z
( complex_double* chi, ElDistMatrix_z x, complex_double* tau );

/* Right application
   ----------------- */
EL_EXPORT ElError ElRightHyperbolicReflector_s
( float* chi, ElMatrix_s x, float* tau );
EL_EXPORT ElError ElRightHyperbolicReflector_d
( double* chi, ElMatrix_d x, double* tau );
EL_EXPORT ElError ElRightHyperbolicReflector_c
( complex_float* chi, ElMatrix_c x, complex_float* tau );
EL_EXPORT ElError ElRightHyperbolicReflector_z
( complex_double* chi, ElMatrix_z x, complex_double* tau );

EL_EXPORT ElError ElRightHyperbolicReflectorDist_s
( float* chi, ElDistMatrix_s x, float* tau );
EL_EXPORT ElError ElRightHyperbolicReflectorDist_d
( double* chi, ElDistMatrix_d x, double* tau );
EL_EXPORT ElError ElRightHyperbolicReflectorDist_c
( complex_float* chi, ElDistMatrix_c x, complex_float* tau );
EL_EXPORT ElError ElRightHyperbolicReflectorDist_z
( complex_double* chi, ElDistMatrix_z x, complex_double* tau );

/* TODO: Wrappers for versions which only communicate within a particular row 
         or column of the process grid */
        
/* Householder reflector
   ===================== */

/* Left application
   ---------------- */
EL_EXPORT ElError ElLeftReflector_s
( float* chi, ElMatrix_s x, float* tau );
EL_EXPORT ElError ElLeftReflector_d
( double* chi, ElMatrix_d x, double* tau );
EL_EXPORT ElError ElLeftReflector_c
( complex_float* chi, ElMatrix_c x, complex_float* tau );
EL_EXPORT ElError ElLeftReflector_z
( complex_double* chi, ElMatrix_z x, complex_double* tau );

EL_EXPORT ElError ElLeftReflectorDist_s
( float* chi, ElDistMatrix_s x, float* tau );
EL_EXPORT ElError ElLeftReflectorDist_d
( double* chi, ElDistMatrix_d x, double* tau );
EL_EXPORT ElError ElLeftReflectorDist_c
( complex_float* chi, ElDistMatrix_c x, complex_float* tau );
EL_EXPORT ElError ElLeftReflectorDist_z
( complex_double* chi, ElDistMatrix_z x, complex_double* tau );

/* Right application
   ----------------- */
EL_EXPORT ElError ElRightReflector_s
( float* chi, ElMatrix_s x, float* tau );
EL_EXPORT ElError ElRightReflector_d
( double* chi, ElMatrix_d x, double* tau );
EL_EXPORT ElError ElRightReflector_c
( complex_float* chi, ElMatrix_c x, complex_float* tau );
EL_EXPORT ElError ElRightReflector_z
( complex_double* chi, ElMatrix_z x, complex_double* tau );

EL_EXPORT ElError ElRightReflectorDist_s
( float* chi, ElDistMatrix_s x, float* tau );
EL_EXPORT ElError ElRightReflectorDist_d
( double* chi, ElDistMatrix_d x, double* tau );
EL_EXPORT ElError ElRightReflectorDist_c
( complex_float* chi, ElDistMatrix_c x, complex_float* tau );
EL_EXPORT ElError ElRightReflectorDist_z
( complex_double* chi, ElDistMatrix_z x, complex_double* tau );

/* TODO: Wrappers for versions which only communicate within a particular row 
         or column of the process grid */

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_REFLECT_C_H */
