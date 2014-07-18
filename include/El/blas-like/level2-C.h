/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS_LEVEL2_C_H
#define EL_BLAS_LEVEL2_C_H

#include "El/core/DistMatrix-C.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Gemv
   ==== */
ElError ElGemv_i
( ElOrientation orientation, 
  ElInt alpha, ElConstMatrix_i A, ElConstMatrix_i x, 
  ElInt beta,                     ElMatrix_i y );
ElError ElGemv_s
( ElOrientation orientation, 
  float alpha, ElConstMatrix_s A, ElConstMatrix_s x, 
  float beta,                     ElMatrix_s y );
ElError ElGemv_d
( ElOrientation orientation, 
  double alpha, ElConstMatrix_d A, ElConstMatrix_d x, 
  double beta,                     ElMatrix_d y );
ElError ElGemv_c
( ElOrientation orientation, 
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c x, 
  complex_float beta,                     ElMatrix_c y );
ElError ElGemv_z
( ElOrientation orientation, 
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z x, 
  complex_double beta,                     ElMatrix_z y );

ElError ElGemvDist_i
( ElOrientation orientation, 
  ElInt alpha, ElConstDistMatrix_i A, ElConstDistMatrix_i x, 
  ElInt beta,                         ElDistMatrix_i y );
ElError ElGemvDist_s
( ElOrientation orientation, 
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s x, 
  float beta,                         ElDistMatrix_s y );
ElError ElGemvDist_d
( ElOrientation orientation, 
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d x, 
  double beta,                         ElDistMatrix_d y );
ElError ElGemvDist_c
( ElOrientation orientation, 
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c x, 
  complex_float beta,                         ElDistMatrix_c y );
ElError ElGemvDist_z
( ElOrientation orientation, 
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z x, 
  complex_double beta,                         ElDistMatrix_z y );

/* Ger
   === */
ElError ElGer_i
( ElInt alpha, ElConstMatrix_i x, ElConstMatrix_i y, ElMatrix_i A );
ElError ElGer_s
( float alpha, ElConstMatrix_s x, ElConstMatrix_s y, ElMatrix_s A );
ElError ElGer_d
( double alpha, ElConstMatrix_d x, ElConstMatrix_d y, ElMatrix_d A );
ElError ElGer_c
( complex_float alpha, ElConstMatrix_c x, ElConstMatrix_c y, ElMatrix_c A );
ElError ElGer_z
( complex_double alpha, ElConstMatrix_z x, ElConstMatrix_z y, ElMatrix_z A );

ElError ElGerDist_i
( ElInt alpha, ElConstDistMatrix_i x, ElConstDistMatrix_i y, 
               ElDistMatrix_i A );
ElError ElGerDist_s
( float alpha, ElConstDistMatrix_s x, ElConstDistMatrix_s y, 
               ElDistMatrix_s A );
ElError ElGerDist_d
( double alpha, ElConstDistMatrix_d x, ElConstDistMatrix_d y, 
                ElDistMatrix_d A );
ElError ElGerDist_c
( complex_float alpha, ElConstDistMatrix_c x, ElConstDistMatrix_c y, 
                       ElDistMatrix_c A );
ElError ElGerDist_z
( complex_double alpha, ElConstDistMatrix_z x, ElConstDistMatrix_z y, 
                        ElDistMatrix_z A );

/* Geru
   ==== */
ElError ElGeru_i
( ElInt alpha, ElConstMatrix_i x, ElConstMatrix_i y, ElMatrix_i A );
ElError ElGeru_s
( float alpha, ElConstMatrix_s x, ElConstMatrix_s y, ElMatrix_s A );
ElError ElGeru_d
( double alpha, ElConstMatrix_d x, ElConstMatrix_d y, ElMatrix_d A );
ElError ElGeru_c
( complex_float alpha, ElConstMatrix_c x, ElConstMatrix_c y, ElMatrix_c A );
ElError ElGeru_z
( complex_double alpha, ElConstMatrix_z x, ElConstMatrix_z y, ElMatrix_z A );

ElError ElGeruDist_i
( ElInt alpha, ElConstDistMatrix_i x, ElConstDistMatrix_i y, 
               ElDistMatrix_i A );
ElError ElGeruDist_s
( float alpha, ElConstDistMatrix_s x, ElConstDistMatrix_s y, 
               ElDistMatrix_s A );
ElError ElGeruDist_d
( double alpha, ElConstDistMatrix_d x, ElConstDistMatrix_d y, 
                ElDistMatrix_d A );
ElError ElGeruDist_c
( complex_float alpha, ElConstDistMatrix_c x, ElConstDistMatrix_c y, 
                       ElDistMatrix_c A );
ElError ElGeruDist_z
( complex_double alpha, ElConstDistMatrix_z x, ElConstDistMatrix_z y, 
                        ElDistMatrix_z A );

/* Hemv
   ==== */
ElError ElHemv_c
( ElUpperOrLower uplo, 
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c x, 
  complex_float beta,                     ElMatrix_c y );
ElError ElHemv_z
( ElUpperOrLower uplo, 
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z x, 
  complex_double beta,                     ElMatrix_z y );

ElError ElHemvDist_c
( ElUpperOrLower uplo, 
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c x, 
  complex_float beta,                         ElDistMatrix_c y );
ElError ElHemvDist_z
( ElUpperOrLower uplo, 
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z x, 
  complex_double beta,                         ElDistMatrix_z y );

/* Her2 
   ==== */
ElError ElHer2_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstMatrix_c x, ElConstMatrix_c y, ElMatrix_c A );
ElError ElHer2_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstMatrix_z x, ElConstMatrix_z y, ElMatrix_z A );

ElError ElHer2Dist_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstDistMatrix_c x, ElConstDistMatrix_c y, ElDistMatrix_c A );
ElError ElHer2Dist_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstDistMatrix_z x, ElConstDistMatrix_z y, ElDistMatrix_z A );

/* Her 
   === */
ElError ElHer_c
( ElUpperOrLower uplo, complex_float alpha, ElConstMatrix_c x, ElMatrix_c A );
ElError ElHer_z
( ElUpperOrLower uplo, complex_double alpha, ElConstMatrix_z x, ElMatrix_z A );

ElError ElHerDist_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstDistMatrix_c x, ElDistMatrix_c A );
ElError ElHerDist_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstDistMatrix_z x, ElDistMatrix_z A );

/* Symv
   ==== */
ElError ElSymv_i
( ElUpperOrLower uplo, 
  ElInt alpha, ElConstMatrix_i A, ElConstMatrix_i x, 
  ElInt beta,                     ElMatrix_i y );
ElError ElSymv_s
( ElUpperOrLower uplo, 
  float alpha, ElConstMatrix_s A, ElConstMatrix_s x, 
  float beta,                     ElMatrix_s y );
ElError ElSymv_d
( ElUpperOrLower uplo, 
  double alpha, ElConstMatrix_d A, ElConstMatrix_d x, 
  double beta,                     ElMatrix_d y );
ElError ElSymv_c
( ElUpperOrLower uplo, 
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c x, 
  complex_float beta,                     ElMatrix_c y );
ElError ElSymv_z
( ElUpperOrLower uplo, 
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z x, 
  complex_double beta,                     ElMatrix_z y );

ElError ElSymvDist_i
( ElUpperOrLower uplo, 
  ElInt alpha, ElConstDistMatrix_i A, ElConstDistMatrix_i x, 
  ElInt beta,                         ElDistMatrix_i y );
ElError ElSymvDist_s
( ElUpperOrLower uplo, 
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s x, 
  float beta,                         ElDistMatrix_s y );
ElError ElSymvDist_d
( ElUpperOrLower uplo, 
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d x, 
  double beta,                         ElDistMatrix_d y );
ElError ElSymvDist_c
( ElUpperOrLower uplo, 
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c x, 
  complex_float beta,                         ElDistMatrix_c y );
ElError ElSymvDist_z
( ElUpperOrLower uplo, 
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z x, 
  complex_double beta,                         ElDistMatrix_z y );

/* Syr2 
   ==== */
ElError ElSyr2_i
( ElUpperOrLower uplo, ElInt alpha, 
  ElConstMatrix_i x, ElConstMatrix_i y, ElMatrix_i A );
ElError ElSyr2_s
( ElUpperOrLower uplo, float alpha, 
  ElConstMatrix_s x, ElConstMatrix_s y, ElMatrix_s A );
ElError ElSyr2_d
( ElUpperOrLower uplo, double alpha, 
  ElConstMatrix_d x, ElConstMatrix_d y, ElMatrix_d A );
ElError ElSyr2_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstMatrix_c x, ElConstMatrix_c y, ElMatrix_c A );
ElError ElSyr2_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstMatrix_z x, ElConstMatrix_z y, ElMatrix_z A );

ElError ElSyr2Dist_i
( ElUpperOrLower uplo, ElInt alpha, 
  ElConstDistMatrix_i x, ElConstDistMatrix_i y, ElDistMatrix_i A );
ElError ElSyr2Dist_s
( ElUpperOrLower uplo, float alpha, 
  ElConstDistMatrix_s x, ElConstDistMatrix_s y, ElDistMatrix_s A );
ElError ElSyr2Dist_d
( ElUpperOrLower uplo, double alpha, 
  ElConstDistMatrix_d x, ElConstDistMatrix_d y, ElDistMatrix_d A );
ElError ElSyr2Dist_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstDistMatrix_c x, ElConstDistMatrix_c y, ElDistMatrix_c A );
ElError ElSyr2Dist_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstDistMatrix_z x, ElConstDistMatrix_z y, ElDistMatrix_z A );

/* Syr 
   === */
ElError ElSyr_i
( ElUpperOrLower uplo, ElInt alpha, ElConstMatrix_i x, ElMatrix_i A );
ElError ElSyr_s
( ElUpperOrLower uplo, float alpha, ElConstMatrix_s x, ElMatrix_s A );
ElError ElSyr_d
( ElUpperOrLower uplo, double alpha, ElConstMatrix_d x, ElMatrix_d A );
ElError ElSyr_c
( ElUpperOrLower uplo, complex_float alpha, ElConstMatrix_c x, ElMatrix_c A );
ElError ElSyr_z
( ElUpperOrLower uplo, complex_double alpha, ElConstMatrix_z x, ElMatrix_z A );

ElError ElSyrDist_i
( ElUpperOrLower uplo, ElInt alpha, 
  ElConstDistMatrix_i x, ElDistMatrix_i A );
ElError ElSyrDist_s
( ElUpperOrLower uplo, float alpha, 
  ElConstDistMatrix_s x, ElDistMatrix_s A );
ElError ElSyrDist_d
( ElUpperOrLower uplo, double alpha, 
  ElConstDistMatrix_d x, ElDistMatrix_d A );
ElError ElSyrDist_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstDistMatrix_c x, ElDistMatrix_c A );
ElError ElSyrDist_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstDistMatrix_z x, ElDistMatrix_z A );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_BLAS_LEVEL2_C_H */
