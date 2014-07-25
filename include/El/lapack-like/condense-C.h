/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LAPACK_CONDENSE_C_H
#define EL_LAPACK_CONDENSE_C_H

#include "El/core/DistMatrix-C.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Bidiag
   ====== */

/* Return the packed reduction to bidiagonal form, B := Q^H A P
   ------------------------------------------------------------ */
ElError ElBidiag_s( ElMatrix_s A, ElMatrix_s tP, ElMatrix_s tQ );
ElError ElBidiag_d( ElMatrix_d A, ElMatrix_d tP, ElMatrix_d tQ );
ElError ElBidiag_c( ElMatrix_c A, ElMatrix_c tP, ElMatrix_c tQ );
ElError ElBidiag_z( ElMatrix_z A, ElMatrix_z tP, ElMatrix_z tQ );

/* NOTE: 'A' must be in a [MC,MR] distribution, while
         'tP' and 'tQ' must be in a [STAR,STAR] distribution */
ElError ElBidiagDist_s
( ElDistMatrix_s A, ElDistMatrix_s tP, ElDistMatrix_s tQ );
ElError ElBidiagDist_d
( ElDistMatrix_d A, ElDistMatrix_d tP, ElDistMatrix_d tQ );
ElError ElBidiagDist_c
( ElDistMatrix_c A, ElDistMatrix_c tP, ElDistMatrix_c tQ );
ElError ElBidiagDist_z
( ElDistMatrix_z A, ElDistMatrix_z tP, ElDistMatrix_z tQ );

/* Only return the condensed bidiagonal matrix, B := Q^H A P
   --------------------------------------------------------- */
ElError ElBidiagOnly_s( ElMatrix_s A );
ElError ElBidiagOnly_d( ElMatrix_d A );
ElError ElBidiagOnly_c( ElMatrix_c A );
ElError ElBidiagOnly_z( ElMatrix_z A );

ElError ElBidiagOnlyDist_s( ElDistMatrix_s A );
ElError ElBidiagOnlyDist_d( ElDistMatrix_d A );
ElError ElBidiagOnlyDist_c( ElDistMatrix_c A );
ElError ElBidiagOnlyDist_z( ElDistMatrix_z A );

/* Apply Q from B := Q^H A P to a set of vectors
   --------------------------------------------- */
ElError ElApplyQAfterBidiag_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_s A, ElConstMatrix_s t, ElMatrix_s B );
ElError ElApplyQAfterBidiag_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_d A, ElConstMatrix_d t, ElMatrix_d B );
ElError ElApplyQAfterBidiag_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_c A, ElConstMatrix_c t, ElMatrix_c B );
ElError ElApplyQAfterBidiag_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_z A, ElConstMatrix_z t, ElMatrix_z B );

/* NOTE: 'A' and 'B' must be in [MC,MR] distributions, while
         't' must be in a [MD,STAR] or [STAR,STAR] distribution */
ElError ElApplyQAfterBidiagDist_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_s A, ElConstDistMatrix_s t, ElDistMatrix_s B );
ElError ElApplyQAfterBidiagDist_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_d A, ElConstDistMatrix_d t, ElDistMatrix_d B );
ElError ElApplyQAfterBidiagDist_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_c A, ElConstDistMatrix_c t, ElDistMatrix_c B );
ElError ElApplyQAfterBidiagDist_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_z A, ElConstDistMatrix_z t, ElDistMatrix_z B );

/* Apply P from B := Q^H A P to a set of vectors
   --------------------------------------------- */
ElError ElApplyPAfterBidiag_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_s A, ElConstMatrix_s t, ElMatrix_s B );
ElError ElApplyPAfterBidiag_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_d A, ElConstMatrix_d t, ElMatrix_d B );
ElError ElApplyPAfterBidiag_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_c A, ElConstMatrix_c t, ElMatrix_c B );
ElError ElApplyPAfterBidiag_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_z A, ElConstMatrix_z t, ElMatrix_z B );

/* NOTE: 'A' and 'B' must be in [MC,MR] distributions, while
         't' must be in a [MD,STAR] or [STAR,STAR] distribution */
ElError ElApplyPAfterBidiagDist_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_s A, ElConstDistMatrix_s t, ElDistMatrix_s B );
ElError ElApplyPAfterBidiagDist_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_d A, ElConstDistMatrix_d t, ElDistMatrix_d B );
ElError ElApplyPAfterBidiagDist_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_c A, ElConstDistMatrix_c t, ElDistMatrix_c B );
ElError ElApplyPAfterBidiagDist_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_z A, ElConstDistMatrix_z t, ElDistMatrix_z B );

/* HermitianTridiag
   ================ */
typedef enum {
  EL_HERMITIAN_TRIDIAG_NORMAL,
  EL_HERMITIAN_TRIDIAG_SQUARE,
  EL_HERMITIAN_TRIDIAG_DEFAULT
} ElHermitianTridiagApproach;

typedef struct {
  ElHermitianTridiagApproach approach;
  ElGridOrderType order;
} ElHermitianTridiagCtrl;
ElError ElHermitianTridiagCtrlDefault( ElHermitianTridiagCtrl* ctrl );

/* Return packed reduction to real symmetric tridiagonal form, T := Q^H A Q 
   ------------------------------------------------------------------------ */
ElError ElHermitianTridiag_s( ElUpperOrLower uplo, ElMatrix_s A, ElMatrix_s t );
ElError ElHermitianTridiag_d( ElUpperOrLower uplo, ElMatrix_d A, ElMatrix_d t );
ElError ElHermitianTridiag_c( ElUpperOrLower uplo, ElMatrix_c A, ElMatrix_c t );
ElError ElHermitianTridiag_z( ElUpperOrLower uplo, ElMatrix_z A, ElMatrix_z t );

/* NOTE: 'A' must be in a [MC,MR] distribution, while
         't' must be in a [STAR,STAR] distribution */
ElError ElHermitianTridiagDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_s t );
ElError ElHermitianTridiagDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_d t );
ElError ElHermitianTridiagDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_c t );
ElError ElHermitianTridiagDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_z t );

/* Expert version
   ^^^^^^^^^^^^^^ */
/* NOTE: 'A' must be in a [MC,MR] distribution, while
         't' must be in a [STAR,STAR] distribution */
ElError ElHermitianTridiagXDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_s t, 
  ElHermitianTridiagCtrl ctrl );
ElError ElHermitianTridiagXDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_d t,
  ElHermitianTridiagCtrl ctrl );
ElError ElHermitianTridiagXDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_c t,
  ElHermitianTridiagCtrl ctrl );
ElError ElHermitianTridiagXDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_z t,
  ElHermitianTridiagCtrl ctrl );

/* Return only the condensed form, T := Q^H A Q
   -------------------------------------------- */
ElError ElHermitianTridiagOnly_s( ElUpperOrLower uplo, ElMatrix_s A );
ElError ElHermitianTridiagOnly_d( ElUpperOrLower uplo, ElMatrix_d A );
ElError ElHermitianTridiagOnly_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElHermitianTridiagOnly_z( ElUpperOrLower uplo, ElMatrix_z A );

/* NOTE: 'A' must be in a [MC,MR] distribution */
ElError ElHermitianTridiagOnlyDist_s( ElUpperOrLower uplo, ElDistMatrix_s A );
ElError ElHermitianTridiagOnlyDist_d( ElUpperOrLower uplo, ElDistMatrix_d A );
ElError ElHermitianTridiagOnlyDist_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElHermitianTridiagOnlyDist_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* Expert version
   ^^^^^^^^^^^^^^ */
/* NOTE: 'A' must be in a [MC,MR] distribution */
ElError ElHermitianTridiagOnlyXDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElHermitianTridiagCtrl ctrl );
ElError ElHermitianTridiagOnlyXDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElHermitianTridiagCtrl ctrl );
ElError ElHermitianTridiagOnlyXDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElHermitianTridiagCtrl ctrl );
ElError ElHermitianTridiagOnlyXDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElHermitianTridiagCtrl ctrl );

/* ApplyQAfterHermitianTridiag 
   --------------------------- */
ElError ElApplyQAfterHermitianTridiag_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_s A, ElConstMatrix_s t, ElMatrix_s B );
ElError ElApplyQAfterHermitianTridiag_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_d A, ElConstMatrix_d t, ElMatrix_d B );
ElError ElApplyQAfterHermitianTridiag_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_c A, ElConstMatrix_c t, ElMatrix_c B );
ElError ElApplyQAfterHermitianTridiag_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_z A, ElConstMatrix_z t, ElMatrix_z B );

/* NOTE: 'A' and 'B' must be in [MC,MR] distributions, while
         't' must be in a [MD,STAR] or [STAR,STAR] distribution */ 
ElError ElApplyQAfterHermitianTridiagDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_s A, ElConstDistMatrix_s t, ElDistMatrix_s B );
ElError ElApplyQAfterHermitianTridiagDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_d A, ElConstDistMatrix_d t, ElDistMatrix_d B );
ElError ElApplyQAfterHermitianTridiagDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_c A, ElConstDistMatrix_c t, ElDistMatrix_c B );
ElError ElApplyQAfterHermitianTridiagDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_z A, ElConstDistMatrix_z t, ElDistMatrix_z B );

/* Hessenberg
   ========== */

/* Packed reduction to Hessenberg form, H := Q^H A Q
   ------------------------------------------------- */
ElError ElHessenberg_s( ElUpperOrLower uplo, ElMatrix_s A, ElMatrix_s t );
ElError ElHessenberg_d( ElUpperOrLower uplo, ElMatrix_d A, ElMatrix_d t );
ElError ElHessenberg_c( ElUpperOrLower uplo, ElMatrix_c A, ElMatrix_c t );
ElError ElHessenberg_z( ElUpperOrLower uplo, ElMatrix_z A, ElMatrix_z t );

/* NOTE: 'A' must be in a [MC,MR] distribution, while
         't' must be in a [STAR,STAR] distribution */
ElError ElHessenbergDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_s t );
ElError ElHessenbergDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_d t );
ElError ElHessenbergDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_c t );
ElError ElHessenbergDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_z t );

/* Only return the similar Hessenberg matrix, H := Q^H A Q
   ------------------------------------------------------- */
ElError ElHessenbergOnly_s( ElUpperOrLower uplo, ElMatrix_s A );
ElError ElHessenbergOnly_d( ElUpperOrLower uplo, ElMatrix_d A );
ElError ElHessenbergOnly_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElHessenbergOnly_z( ElUpperOrLower uplo, ElMatrix_z A );

/* NOTE: 'A' must be in a [MC,MR] distribution */
ElError ElHessenbergOnlyDist_s( ElUpperOrLower uplo, ElDistMatrix_s A );
ElError ElHessenbergOnlyDist_d( ElUpperOrLower uplo, ElDistMatrix_d A );
ElError ElHessenbergOnlyDist_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElHessenbergOnlyDist_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* Apply Q from a Hessenberg decomposition, H := Q^H A Q
   ----------------------------------------------------- */
ElError ElApplyQAfterHessenberg_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_s A, ElConstMatrix_s t, ElMatrix_s B );
ElError ElApplyQAfterHessenberg_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_d A, ElConstMatrix_d t, ElMatrix_d B );
ElError ElApplyQAfterHessenberg_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_c A, ElConstMatrix_c t, ElMatrix_c B );
ElError ElApplyQAfterHessenberg_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_z A, ElConstMatrix_z t, ElMatrix_z B );

/* NOTE: 'A' and 'H' must be in a [MC,MR] distribution, while
         't' must be in either a [MD,STAR] or [STAR,STAR] distribution */ 
ElError ElApplyQAfterHessenbergDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_s A, ElConstDistMatrix_s t, ElDistMatrix_s B );
ElError ElApplyQAfterHessenbergDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_d A, ElConstDistMatrix_d t, ElDistMatrix_d B );
ElError ElApplyQAfterHessenbergDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_c A, ElConstDistMatrix_c t, ElDistMatrix_c B );
ElError ElApplyQAfterHessenbergDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_z A, ElConstDistMatrix_z t, ElDistMatrix_z B );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_LAPACK_CONDENSE_C_H */
