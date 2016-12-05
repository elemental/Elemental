/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
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

/* SOC Identity
   ------------ */
EL_EXPORT ElError ElSOCIdentity_s
( ElMatrix_s x, ElConstMatrix_i orders, ElConstMatrix_i firstInds );
EL_EXPORT ElError ElSOCIdentity_d
( ElMatrix_d x, ElConstMatrix_i orders, ElConstMatrix_i firstInds );

EL_EXPORT ElError ElSOCIdentityDist_s
( ElDistMatrix_s x,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds );
EL_EXPORT ElError ElSOCIdentityDist_d
( ElDistMatrix_d x,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds );

EL_EXPORT ElError ElSOCIdentityDistMultiVec_s
( ElDistMultiVec_s x,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds );
EL_EXPORT ElError ElSOCIdentityDistMultiVec_d
( ElDistMultiVec_d x,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds );

/* SOC Dots
   -------- */
EL_EXPORT ElError ElSOCDots_s
( ElConstMatrix_s x, ElConstMatrix_s y, ElMatrix_s z,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );
EL_EXPORT ElError ElSOCDots_d
( ElConstMatrix_d x, ElConstMatrix_d y, ElMatrix_d z,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );

EL_EXPORT ElError ElSOCDotsDist_s
( ElConstDistMatrix_s x, ElConstDistMatrix_s y, ElDistMatrix_s z,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElSOCDotsDist_d
( ElConstDistMatrix_d x, ElConstDistMatrix_d y, ElDistMatrix_d z,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );

EL_EXPORT ElError ElSOCDotsDistMultiVec_s
( ElConstDistMultiVec_s x, ElConstDistMultiVec_s y, ElDistMultiVec_s z,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElSOCDotsDistMultiVec_d
( ElConstDistMultiVec_d x, ElConstDistMultiVec_d y, ElDistMultiVec_d z,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  ElInt cutoff );

/* Cone Broadcast
   -------------- */
EL_EXPORT ElError ElConeBroadcast_s
( ElMatrix_s x,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );
EL_EXPORT ElError ElConeBroadcast_d
( ElMatrix_d x,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );

EL_EXPORT ElError ElConeBroadcastDist_s
( ElDistMatrix_s x,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElConeBroadcastDist_d
( ElDistMatrix_d x,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );

EL_EXPORT ElError ElConeBroadcastDistMultiVec_s
( ElDistMultiVec_s x,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElConeBroadcastDistMultiVec_d
( ElDistMultiVec_d x,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  ElInt cutoff );

/* SOC Reflection
   -------------- */
EL_EXPORT ElError ElSOCReflect_s
( ElMatrix_s x,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );
EL_EXPORT ElError ElSOCReflect_d
( ElMatrix_d x,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );

EL_EXPORT ElError ElSOCReflectDist_s
( ElDistMatrix_s x,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds );
EL_EXPORT ElError ElSOCReflectDist_d
( ElDistMatrix_d x,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds );

EL_EXPORT ElError ElSOCReflectDistMultiVec_s
( ElDistMultiVec_s x,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds );
EL_EXPORT ElError ElSOCReflectDistMultiVec_d
( ElDistMultiVec_d x,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds );

/* SOC Determinants
   ---------------- */
EL_EXPORT ElError ElSOCDets_s
( ElConstMatrix_s x, ElMatrix_s d,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );
EL_EXPORT ElError ElSOCDets_d
( ElConstMatrix_d x, ElMatrix_d d,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );

EL_EXPORT ElError ElSOCDetsDist_s
( ElConstDistMatrix_s x, ElDistMatrix_s d,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElSOCDetsDist_d
( ElConstDistMatrix_d x, ElDistMatrix_d d,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );

EL_EXPORT ElError ElSOCDetsDistMultiVec_s
( ElConstDistMultiVec_s x, ElDistMultiVec_s d,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElSOCDetsDistMultiVec_d
( ElConstDistMultiVec_d x, ElDistMultiVec_d d,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  ElInt cutoff );

/* Num non-SOC
   ----------- */
EL_EXPORT ElError ElNumNonSOC_s
( ElConstMatrix_s x,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds,
  ElInt* numNonSOC );
EL_EXPORT ElError ElNumNonSOC_d
( ElConstMatrix_d x,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds,
  ElInt* numNonSOC );

EL_EXPORT ElError ElNumNonSOCDist_s
( ElConstDistMatrix_s x,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff, ElInt* numNonSOC );
EL_EXPORT ElError ElNumNonSOCDist_d
( ElConstDistMatrix_d x,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff, ElInt* numNonSOC );

EL_EXPORT ElError ElNumNonSOCDistMultiVec_s
( ElConstDistMultiVec_s x,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  ElInt cutoff, ElInt* numNonSOC );
EL_EXPORT ElError ElNumNonSOCDistMultiVec_d
( ElConstDistMultiVec_d x,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  ElInt cutoff, ElInt* numNonSOC );

/* SOC Apply
   --------- */
EL_EXPORT ElError ElSOCApply_s
( ElConstMatrix_s x, ElConstMatrix_s y, ElMatrix_s z,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );
EL_EXPORT ElError ElSOCApply_d
( ElConstMatrix_d x, ElConstMatrix_d y, ElMatrix_d z,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );

EL_EXPORT ElError ElSOCApplyDist_s
( ElConstDistMatrix_s x, ElConstDistMatrix_s y, ElDistMatrix_s z,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElSOCApplyDist_d
( ElConstDistMatrix_d x, ElConstDistMatrix_d y, ElDistMatrix_d z,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );

EL_EXPORT ElError ElSOCApplyDist_s
( ElConstDistMatrix_s x, ElConstDistMatrix_s y, ElDistMatrix_s z,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElSOCApplyDist_d
( ElConstDistMatrix_d x, ElConstDistMatrix_d y, ElDistMatrix_d z,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );

/* SOC Apply quadratic
   ------------------- */
EL_EXPORT ElError ElSOCApplyQuadratic_s
( ElConstMatrix_s x, ElConstMatrix_s y, ElMatrix_s z,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );
EL_EXPORT ElError ElSOCApplyQuadratic_d
( ElConstMatrix_d x, ElConstMatrix_d y, ElMatrix_d z,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );

EL_EXPORT ElError ElSOCApplyQuadraticDist_s
( ElConstDistMatrix_s x, ElConstDistMatrix_s y, ElDistMatrix_s z,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElSOCApplyQuadraticDist_d
( ElConstDistMatrix_d x, ElConstDistMatrix_d y, ElDistMatrix_d z,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );

EL_EXPORT ElError ElSOCApplyQuadraticDist_s
( ElConstDistMatrix_s x, ElConstDistMatrix_s y, ElDistMatrix_s z,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElSOCApplyQuadraticDist_d
( ElConstDistMatrix_d x, ElConstDistMatrix_d y, ElDistMatrix_d z,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );

/* SOC Inverse
   ----------- */
EL_EXPORT ElError ElSOCInverse_s
( ElConstMatrix_s x, ElMatrix_s xInv,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );
EL_EXPORT ElError ElSOCInverse_d
( ElConstMatrix_d x, ElMatrix_d xInv,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );

EL_EXPORT ElError ElSOCInverseDist_s
( ElConstDistMatrix_s x, ElDistMatrix_s xInv,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElSOCInverseDist_d
( ElConstDistMatrix_d x, ElDistMatrix_d xInv,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );

EL_EXPORT ElError ElSOCInverseDistMultiVec_s
( ElConstDistMultiVec_s x, ElDistMultiVec_s xInv,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElSOCInverseDistMultiVec_d
( ElConstDistMultiVec_d x, ElDistMultiVec_d xInv,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  ElInt cutoff );

/* SOC Square-root
   --------------- */
EL_EXPORT ElError ElSOCSquareRoot_s
( ElConstMatrix_s x, ElMatrix_s xRoot,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );
EL_EXPORT ElError ElSOCSquareRoot_d
( ElConstMatrix_d x, ElMatrix_d xRoot,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );

EL_EXPORT ElError ElSOCSquareRootDist_s
( ElConstDistMatrix_s x, ElDistMatrix_s xRoot,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElSOCSquareRootDist_d
( ElConstDistMatrix_d x, ElDistMatrix_d xRoot,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );

EL_EXPORT ElError ElSOCSquareRootDistMultiVec_s
( ElConstDistMultiVec_s x, ElDistMultiVec_s xRoot,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElSOCSquareRootDistMultiVec_d
( ElConstDistMultiVec_d x, ElDistMultiVec_d xRoot,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  ElInt cutoff );

/* SOC Nesterov-Todd
   ----------------- */
EL_EXPORT ElError ElSOCNesterovTodd_s
( ElConstMatrix_s x, ElConstMatrix_s y, ElMatrix_s z,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );
EL_EXPORT ElError ElSOCNesterovTodd_d
( ElConstMatrix_d x, ElConstMatrix_d y, ElMatrix_d z,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds );

EL_EXPORT ElError ElSOCNesterovToddDist_s
( ElConstDistMatrix_s x, ElConstDistMatrix_s y, ElDistMatrix_s z,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElSOCNesterovToddDist_d
( ElConstDistMatrix_d x, ElConstDistMatrix_d y, ElDistMatrix_d z,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  ElInt cutoff );

EL_EXPORT ElError ElSOCNesterovToddDistMultiVec_s
( ElConstDistMultiVec_s x, ElConstDistMultiVec_s y, ElDistMultiVec_s z,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  ElInt cutoff );
EL_EXPORT ElError ElSOCNesterovToddDistMultiVec_d
( ElConstDistMultiVec_d x, ElConstDistMultiVec_d y, ElDistMultiVec_d z,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  ElInt cutoff );

/* Max step in SOC
   --------------- */
EL_EXPORT ElError ElMaxStepInSOC_s
( ElConstMatrix_s x, ElConstMatrix_s y,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds,
  float upperBound, float* alpha );
EL_EXPORT ElError ElMaxStepInSOC_d
( ElConstMatrix_d x, ElConstMatrix_d y,
  ElConstMatrix_i orders, ElConstMatrix_i firstInds,
  double upperBound, double* alpha );

EL_EXPORT ElError ElMaxStepInSOCDist_s
( ElConstDistMatrix_s x, ElConstDistMatrix_s y,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  float upperBound, ElInt cutoff, float* alpha );
EL_EXPORT ElError ElMaxStepInSOCDist_d
( ElConstDistMatrix_d x, ElConstDistMatrix_d y,
  ElConstDistMatrix_i orders, ElConstDistMatrix_i firstInds,
  double upperBound, ElInt cutoff, double* alpha );

EL_EXPORT ElError ElMaxStepInSOCDistMultiVec_s
( ElConstDistMultiVec_s x, ElConstDistMultiVec_s y,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  float upperBound, ElInt cutoff, float* alpha );
EL_EXPORT ElError ElMaxStepInSOCDistMultiVec_d
( ElConstDistMultiVec_d x, ElConstDistMultiVec_d y,
  ElConstDistMultiVec_i orders, ElConstDistMultiVec_i firstInds,
  double upperBound, ElInt cutoff, double* alpha );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_OPTIMIZATION_UTIL_C_H */
