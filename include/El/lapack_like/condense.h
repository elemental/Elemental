/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LAPACK_CONDENSE_C_H
#define EL_LAPACK_CONDENSE_C_H

#include "El/core/DistMatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Bidiag
   ====== */

/* Return the packed reduction to bidiagonal form, B := Q^H A P
   ------------------------------------------------------------ */
EL_EXPORT ElError ElBidiag_s( ElMatrix_s A, ElMatrix_s tP, ElMatrix_s tQ );
EL_EXPORT ElError ElBidiag_d( ElMatrix_d A, ElMatrix_d tP, ElMatrix_d tQ );
EL_EXPORT ElError ElBidiag_c( ElMatrix_c A, ElMatrix_c tP, ElMatrix_c tQ );
EL_EXPORT ElError ElBidiag_z( ElMatrix_z A, ElMatrix_z tP, ElMatrix_z tQ );

EL_EXPORT ElError ElBidiagDist_s
( ElDistMatrix_s A, ElDistMatrix_s tP, ElDistMatrix_s tQ );
EL_EXPORT ElError ElBidiagDist_d
( ElDistMatrix_d A, ElDistMatrix_d tP, ElDistMatrix_d tQ );
EL_EXPORT ElError ElBidiagDist_c
( ElDistMatrix_c A, ElDistMatrix_c tP, ElDistMatrix_c tQ );
EL_EXPORT ElError ElBidiagDist_z
( ElDistMatrix_z A, ElDistMatrix_z tP, ElDistMatrix_z tQ );

/* Only return the condensed bidiagonal matrix, B := Q^H A P
   --------------------------------------------------------- */
EL_EXPORT ElError ElBidiagOnly_s( ElMatrix_s A );
EL_EXPORT ElError ElBidiagOnly_d( ElMatrix_d A );
EL_EXPORT ElError ElBidiagOnly_c( ElMatrix_c A );
EL_EXPORT ElError ElBidiagOnly_z( ElMatrix_z A );

EL_EXPORT ElError ElBidiagOnlyDist_s( ElDistMatrix_s A );
EL_EXPORT ElError ElBidiagOnlyDist_d( ElDistMatrix_d A );
EL_EXPORT ElError ElBidiagOnlyDist_c( ElDistMatrix_c A );
EL_EXPORT ElError ElBidiagOnlyDist_z( ElDistMatrix_z A );

/* Apply Q from B := Q^H A P to a set of vectors
   --------------------------------------------- */
EL_EXPORT ElError ElApplyQAfterBidiag_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_s A, ElConstMatrix_s t, ElMatrix_s B );
EL_EXPORT ElError ElApplyQAfterBidiag_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_d A, ElConstMatrix_d t, ElMatrix_d B );
EL_EXPORT ElError ElApplyQAfterBidiag_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_c A, ElConstMatrix_c t, ElMatrix_c B );
EL_EXPORT ElError ElApplyQAfterBidiag_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_z A, ElConstMatrix_z t, ElMatrix_z B );

EL_EXPORT ElError ElApplyQAfterBidiagDist_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_s A, ElConstDistMatrix_s t, ElDistMatrix_s B );
EL_EXPORT ElError ElApplyQAfterBidiagDist_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_d A, ElConstDistMatrix_d t, ElDistMatrix_d B );
EL_EXPORT ElError ElApplyQAfterBidiagDist_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_c A, ElConstDistMatrix_c t, ElDistMatrix_c B );
EL_EXPORT ElError ElApplyQAfterBidiagDist_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_z A, ElConstDistMatrix_z t, ElDistMatrix_z B );

/* Apply P from B := Q^H A P to a set of vectors
   --------------------------------------------- */
EL_EXPORT ElError ElApplyPAfterBidiag_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_s A, ElConstMatrix_s t, ElMatrix_s B );
EL_EXPORT ElError ElApplyPAfterBidiag_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_d A, ElConstMatrix_d t, ElMatrix_d B );
EL_EXPORT ElError ElApplyPAfterBidiag_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_c A, ElConstMatrix_c t, ElMatrix_c B );
EL_EXPORT ElError ElApplyPAfterBidiag_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_z A, ElConstMatrix_z t, ElMatrix_z B );

EL_EXPORT ElError ElApplyPAfterBidiagDist_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_s A, ElConstDistMatrix_s t, ElDistMatrix_s B );
EL_EXPORT ElError ElApplyPAfterBidiagDist_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_d A, ElConstDistMatrix_d t, ElDistMatrix_d B );
EL_EXPORT ElError ElApplyPAfterBidiagDist_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_c A, ElConstDistMatrix_c t, ElDistMatrix_c B );
EL_EXPORT ElError ElApplyPAfterBidiagDist_z
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
  ElSymvCtrl symvCtrl;
} ElHermitianTridiagCtrl;
EL_EXPORT ElError 
ElHermitianTridiagCtrlDefault_s( ElHermitianTridiagCtrl* ctrl );
EL_EXPORT ElError 
ElHermitianTridiagCtrlDefault_d( ElHermitianTridiagCtrl* ctrl );
EL_EXPORT ElError 
ElHermitianTridiagCtrlDefault_c( ElHermitianTridiagCtrl* ctrl );
EL_EXPORT ElError 
ElHermitianTridiagCtrlDefault_z( ElHermitianTridiagCtrl* ctrl );

/* Return packed reduction to real symmetric tridiagonal form, T := Q^H A Q 
   ------------------------------------------------------------------------ */
EL_EXPORT ElError ElHermitianTridiag_s
( ElUpperOrLower uplo, ElMatrix_s A, ElMatrix_s t );
EL_EXPORT ElError ElHermitianTridiag_d
( ElUpperOrLower uplo, ElMatrix_d A, ElMatrix_d t );
EL_EXPORT ElError ElHermitianTridiag_c
( ElUpperOrLower uplo, ElMatrix_c A, ElMatrix_c t );
EL_EXPORT ElError ElHermitianTridiag_z
( ElUpperOrLower uplo, ElMatrix_z A, ElMatrix_z t );

EL_EXPORT ElError ElHermitianTridiagDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_s t );
EL_EXPORT ElError ElHermitianTridiagDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_d t );
EL_EXPORT ElError ElHermitianTridiagDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_c t );
EL_EXPORT ElError ElHermitianTridiagDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_z t );

/* Expert version
   ^^^^^^^^^^^^^^ */
EL_EXPORT ElError ElHermitianTridiagXDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_s t, 
  ElHermitianTridiagCtrl ctrl );
EL_EXPORT ElError ElHermitianTridiagXDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_d t,
  ElHermitianTridiagCtrl ctrl );
EL_EXPORT ElError ElHermitianTridiagXDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_c t,
  ElHermitianTridiagCtrl ctrl );
EL_EXPORT ElError ElHermitianTridiagXDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_z t,
  ElHermitianTridiagCtrl ctrl );

/* Return only the condensed form, T := Q^H A Q
   -------------------------------------------- */
EL_EXPORT ElError ElHermitianTridiagOnly_s( ElUpperOrLower uplo, ElMatrix_s A );
EL_EXPORT ElError ElHermitianTridiagOnly_d( ElUpperOrLower uplo, ElMatrix_d A );
EL_EXPORT ElError ElHermitianTridiagOnly_c( ElUpperOrLower uplo, ElMatrix_c A );
EL_EXPORT ElError ElHermitianTridiagOnly_z( ElUpperOrLower uplo, ElMatrix_z A );

EL_EXPORT ElError ElHermitianTridiagOnlyDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A );
EL_EXPORT ElError ElHermitianTridiagOnlyDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A );
EL_EXPORT ElError ElHermitianTridiagOnlyDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A );
EL_EXPORT ElError ElHermitianTridiagOnlyDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A );

/* Expert version
   ^^^^^^^^^^^^^^ */
EL_EXPORT ElError ElHermitianTridiagOnlyXDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElHermitianTridiagCtrl ctrl );
EL_EXPORT ElError ElHermitianTridiagOnlyXDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElHermitianTridiagCtrl ctrl );
EL_EXPORT ElError ElHermitianTridiagOnlyXDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElHermitianTridiagCtrl ctrl );
EL_EXPORT ElError ElHermitianTridiagOnlyXDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElHermitianTridiagCtrl ctrl );

/* ApplyQAfterHermitianTridiag 
   --------------------------- */
EL_EXPORT ElError ElApplyQAfterHermitianTridiag_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_s A, ElConstMatrix_s t, ElMatrix_s B );
EL_EXPORT ElError ElApplyQAfterHermitianTridiag_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_d A, ElConstMatrix_d t, ElMatrix_d B );
EL_EXPORT ElError ElApplyQAfterHermitianTridiag_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_c A, ElConstMatrix_c t, ElMatrix_c B );
EL_EXPORT ElError ElApplyQAfterHermitianTridiag_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_z A, ElConstMatrix_z t, ElMatrix_z B );

EL_EXPORT ElError ElApplyQAfterHermitianTridiagDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_s A, ElConstDistMatrix_s t, ElDistMatrix_s B );
EL_EXPORT ElError ElApplyQAfterHermitianTridiagDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_d A, ElConstDistMatrix_d t, ElDistMatrix_d B );
EL_EXPORT ElError ElApplyQAfterHermitianTridiagDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_c A, ElConstDistMatrix_c t, ElDistMatrix_c B );
EL_EXPORT ElError ElApplyQAfterHermitianTridiagDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_z A, ElConstDistMatrix_z t, ElDistMatrix_z B );

/* Hessenberg
   ========== */

/* Packed reduction to Hessenberg form, H := Q^H A Q
   ------------------------------------------------- */
EL_EXPORT ElError ElHessenberg_s
( ElUpperOrLower uplo, ElMatrix_s A, ElMatrix_s t );
EL_EXPORT ElError ElHessenberg_d
( ElUpperOrLower uplo, ElMatrix_d A, ElMatrix_d t );
EL_EXPORT ElError ElHessenberg_c
( ElUpperOrLower uplo, ElMatrix_c A, ElMatrix_c t );
EL_EXPORT ElError ElHessenberg_z
( ElUpperOrLower uplo, ElMatrix_z A, ElMatrix_z t );

EL_EXPORT ElError ElHessenbergDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_s t );
EL_EXPORT ElError ElHessenbergDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_d t );
EL_EXPORT ElError ElHessenbergDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_c t );
EL_EXPORT ElError ElHessenbergDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_z t );

/* Only return the similar Hessenberg matrix, H := Q^H A Q
   ------------------------------------------------------- */
EL_EXPORT ElError ElHessenbergOnly_s( ElUpperOrLower uplo, ElMatrix_s A );
EL_EXPORT ElError ElHessenbergOnly_d( ElUpperOrLower uplo, ElMatrix_d A );
EL_EXPORT ElError ElHessenbergOnly_c( ElUpperOrLower uplo, ElMatrix_c A );
EL_EXPORT ElError ElHessenbergOnly_z( ElUpperOrLower uplo, ElMatrix_z A );

EL_EXPORT ElError ElHessenbergOnlyDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A );
EL_EXPORT ElError ElHessenbergOnlyDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A );
EL_EXPORT ElError ElHessenbergOnlyDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A );
EL_EXPORT ElError ElHessenbergOnlyDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A );

/* Apply Q from a Hessenberg decomposition, H := Q^H A Q
   ----------------------------------------------------- */
EL_EXPORT ElError ElApplyQAfterHessenberg_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_s A, ElConstMatrix_s t, ElMatrix_s B );
EL_EXPORT ElError ElApplyQAfterHessenberg_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_d A, ElConstMatrix_d t, ElMatrix_d B );
EL_EXPORT ElError ElApplyQAfterHessenberg_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_c A, ElConstMatrix_c t, ElMatrix_c B );
EL_EXPORT ElError ElApplyQAfterHessenberg_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_z A, ElConstMatrix_z t, ElMatrix_z B );

EL_EXPORT ElError ElApplyQAfterHessenbergDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_s A, ElConstDistMatrix_s t, ElDistMatrix_s B );
EL_EXPORT ElError ElApplyQAfterHessenbergDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_d A, ElConstDistMatrix_d t, ElDistMatrix_d B );
EL_EXPORT ElError ElApplyQAfterHessenbergDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_c A, ElConstDistMatrix_c t, ElDistMatrix_c B );
EL_EXPORT ElError ElApplyQAfterHessenbergDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_z A, ElConstDistMatrix_z t, ElDistMatrix_z B );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_LAPACK_CONDENSE_C_H */
