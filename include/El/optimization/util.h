/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_OPTIMIZATION_UTIL_C_H
#define EL_OPTIMIZATION_UTIL_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Coherence
   --------- */
EL_EXPORT ElError ElCoherence_s( ElConstMatrix_s A, float* coherence );
EL_EXPORT ElError ElCoherence_d( ElConstMatrix_d A, double* coherence );
EL_EXPORT ElError ElCoherence_c( ElConstMatrix_c A, float* coherence );
EL_EXPORT ElError ElCoherence_z( ElConstMatrix_z A, double* coherence );

EL_EXPORT ElError ElCoherenceDist_s( ElConstDistMatrix_s A, float* coherence );
EL_EXPORT ElError ElCoherenceDist_d( ElConstDistMatrix_d A, double* coherence );
EL_EXPORT ElError ElCoherenceDist_c( ElConstDistMatrix_c A, float* coherence );
EL_EXPORT ElError ElCoherenceDist_z( ElConstDistMatrix_z A, double* coherence );

/* Covariance
   ---------- */
EL_EXPORT ElError ElCovariance_s( ElConstMatrix_s D, ElMatrix_s S );
EL_EXPORT ElError ElCovariance_d( ElConstMatrix_d D, ElMatrix_d S );
EL_EXPORT ElError ElCovariance_c( ElConstMatrix_c D, ElMatrix_c S );
EL_EXPORT ElError ElCovariance_z( ElConstMatrix_z D, ElMatrix_z S );

EL_EXPORT ElError ElCovarianceDist_s( ElConstDistMatrix_s D, ElDistMatrix_s S );
EL_EXPORT ElError ElCovarianceDist_d( ElConstDistMatrix_d D, ElDistMatrix_d S );
EL_EXPORT ElError ElCovarianceDist_c( ElConstDistMatrix_c D, ElDistMatrix_c S );
EL_EXPORT ElError ElCovarianceDist_z( ElConstDistMatrix_z D, ElDistMatrix_z S );

/* Log barrier
   ----------- */
EL_EXPORT ElError ElLogBarrier_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* barrier );
EL_EXPORT ElError ElLogBarrier_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* barrier );
EL_EXPORT ElError ElLogBarrier_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* barrier );
EL_EXPORT ElError ElLogBarrier_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* barrier );

EL_EXPORT ElError ElLogBarrierDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* barrier );
EL_EXPORT ElError ElLogBarrierDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* barrier );
EL_EXPORT ElError ElLogBarrierDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* barrier );
EL_EXPORT ElError ElLogBarrierDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* barrier );

/* TODO: Version which allows overwriting? */

/* Log-det divergence
   ------------------ */
EL_EXPORT ElError ElLogDetDiv_s
( ElUpperOrLower uplo, ElConstMatrix_s A, ElConstMatrix_s B, float* div );
EL_EXPORT ElError ElLogDetDiv_d
( ElUpperOrLower uplo, ElConstMatrix_d A, ElConstMatrix_d B, double* div );
EL_EXPORT ElError ElLogDetDiv_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElConstMatrix_c B, float* div );
EL_EXPORT ElError ElLogDetDiv_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElConstMatrix_z B, double* div );

EL_EXPORT ElError ElLogDetDivDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, ElConstDistMatrix_s B, 
  float* div );
EL_EXPORT ElError ElLogDetDivDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, ElConstDistMatrix_d B, 
  double* div );
EL_EXPORT ElError ElLogDetDivDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElConstDistMatrix_c B, 
  float* div );
EL_EXPORT ElError ElLogDetDivDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElConstDistMatrix_z B, 
  double* div );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_OPTIMIZATION_UTIL_C_H */
