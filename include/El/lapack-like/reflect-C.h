/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_REFLECT_C_H
#define EL_REFLECT_C_H

#include "El/core/DistMatrix-C.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Apply packed reflectors
   ======================= */
ElError ElApplyPackedReflectors_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, ElInt offset, 
  ElConstMatrix_s H, ElConstMatrix_s t, ElMatrix_s A );
ElError ElApplyPackedReflectors_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, ElInt offset, 
  ElConstMatrix_d H, ElConstMatrix_d t, ElMatrix_d A );
ElError ElApplyPackedReflectors_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, 
  ElConjugation conjugation, ElInt offset, 
  ElConstMatrix_c H, ElConstMatrix_c t, ElMatrix_c A );
ElError ElApplyPackedReflectors_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, 
  ElConjugation conjugation, ElInt offset, 
  ElConstMatrix_z H, ElConstMatrix_z t, ElMatrix_z A );

ElError ElApplyPackedReflectorsDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, ElInt offset, 
  ElConstDistMatrix_s H, ElConstDistMatrix_s t, ElDistMatrix_s A );
ElError ElApplyPackedReflectorsDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, ElInt offset, 
  ElConstDistMatrix_d H, ElConstDistMatrix_d t, ElDistMatrix_d A );
ElError ElApplyPackedReflectorsDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, 
  ElConjugation conjugation, ElInt offset, 
  ElConstDistMatrix_c H, ElConstDistMatrix_c t, ElDistMatrix_c A );
ElError ElApplyPackedReflectorsDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElVerticalOrHorizontal dir, ElForwardOrBackward order, 
  ElConjugation conjugation, ElInt offset, 
  ElConstDistMatrix_z H, ElConstDistMatrix_z t, ElDistMatrix_z A );

/* Expand packed reflectors
   ======================== */
ElError ElExpandPackedReflectors_s
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElInt offset, 
  ElMatrix_s H, ElConstMatrix_s t );
ElError ElExpandPackedReflectors_d
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElInt offset, 
  ElMatrix_d H, ElConstMatrix_d t );
ElError ElExpandPackedReflectors_c
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElConjugation conjugation,
  ElInt offset, ElMatrix_c H, ElConstMatrix_c t );
ElError ElExpandPackedReflectors_z
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElConjugation conjugation,
  ElInt offset, ElMatrix_z H, ElConstMatrix_z t );

ElError ElExpandPackedReflectorsDist_s
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElInt offset, 
  ElDistMatrix_s H, ElConstDistMatrix_s t );
ElError ElExpandPackedReflectorsDist_d
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElInt offset, 
  ElDistMatrix_d H, ElConstDistMatrix_d t );
ElError ElExpandPackedReflectorsDist_c
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElConjugation conjugation,
  ElInt offset, ElDistMatrix_c H, ElConstDistMatrix_c t );
ElError ElExpandPackedReflectorsDist_z
( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElConjugation conjugation,
  ElInt offset, ElDistMatrix_z H, ElConstDistMatrix_z t );

/* Hyperbolic reflector
   ==================== */

/* Left application
   ---------------- */
ElError ElLeftHyperbolicReflector_s
( float* chi, ElMatrix_s x, float* tau );
ElError ElLeftHyperbolicReflector_d
( double* chi, ElMatrix_d x, double* tau );
ElError ElLeftHyperbolicReflector_c
( complex_float* chi, ElMatrix_c x, complex_float* tau );
ElError ElLeftHyperbolicReflector_z
( complex_double* chi, ElMatrix_z x, complex_double* tau );

ElError ElLeftHyperbolicReflectorDist_s
( float* chi, ElDistMatrix_s x, float* tau );
ElError ElLeftHyperbolicReflectorDist_d
( double* chi, ElDistMatrix_d x, double* tau );
ElError ElLeftHyperbolicReflectorDist_c
( complex_float* chi, ElDistMatrix_c x, complex_float* tau );
ElError ElLeftHyperbolicReflectorDist_z
( complex_double* chi, ElDistMatrix_z x, complex_double* tau );

/* Right application
   ----------------- */
ElError ElRightHyperbolicReflector_s
( float* chi, ElMatrix_s x, float* tau );
ElError ElRightHyperbolicReflector_d
( double* chi, ElMatrix_d x, double* tau );
ElError ElRightHyperbolicReflector_c
( complex_float* chi, ElMatrix_c x, complex_float* tau );
ElError ElRightHyperbolicReflector_z
( complex_double* chi, ElMatrix_z x, complex_double* tau );

ElError ElRightHyperbolicReflectorDist_s
( float* chi, ElDistMatrix_s x, float* tau );
ElError ElRightHyperbolicReflectorDist_d
( double* chi, ElDistMatrix_d x, double* tau );
ElError ElRightHyperbolicReflectorDist_c
( complex_float* chi, ElDistMatrix_c x, complex_float* tau );
ElError ElRightHyperbolicReflectorDist_z
( complex_double* chi, ElDistMatrix_z x, complex_double* tau );

/* TODO: Wrappers for versions which only communicate within a particular row 
         or column of the process grid */
        
/* Householder reflector
   ===================== */

/* Left application
   ---------------- */
ElError ElLeftReflector_s
( float* chi, ElMatrix_s x, float* tau );
ElError ElLeftReflector_d
( double* chi, ElMatrix_d x, double* tau );
ElError ElLeftReflector_c
( complex_float* chi, ElMatrix_c x, complex_float* tau );
ElError ElLeftReflector_z
( complex_double* chi, ElMatrix_z x, complex_double* tau );

ElError ElLeftReflectorDist_s
( float* chi, ElDistMatrix_s x, float* tau );
ElError ElLeftReflectorDist_d
( double* chi, ElDistMatrix_d x, double* tau );
ElError ElLeftReflectorDist_c
( complex_float* chi, ElDistMatrix_c x, complex_float* tau );
ElError ElLeftReflectorDist_z
( complex_double* chi, ElDistMatrix_z x, complex_double* tau );

/* Right application
   ----------------- */
ElError ElRightReflector_s
( float* chi, ElMatrix_s x, float* tau );
ElError ElRightReflector_d
( double* chi, ElMatrix_d x, double* tau );
ElError ElRightReflector_c
( complex_float* chi, ElMatrix_c x, complex_float* tau );
ElError ElRightReflector_z
( complex_double* chi, ElMatrix_z x, complex_double* tau );

ElError ElRightReflectorDist_s
( float* chi, ElDistMatrix_s x, float* tau );
ElError ElRightReflectorDist_d
( double* chi, ElDistMatrix_d x, double* tau );
ElError ElRightReflectorDist_c
( complex_float* chi, ElDistMatrix_c x, complex_float* tau );
ElError ElRightReflectorDist_z
( complex_double* chi, ElDistMatrix_z x, complex_double* tau );

/* TODO: Wrappers for versions which only communicate within a particular row 
         or column of the process grid */

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_REFLECT_C_H */
