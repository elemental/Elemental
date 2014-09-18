/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_OPTIMIZATION_C_H
#define EL_OPTIMIZATION_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Basis pursuit
   ============= */
/* min || z ||_1 such that A z = b */

ElError ElBasisPursuit_s
( ElConstMatrix_s A, ElConstMatrix_s b, ElMatrix_s z, ElInt* numIts );
ElError ElBasisPursuit_d
( ElConstMatrix_d A, ElConstMatrix_d b, ElMatrix_d z, ElInt* numIts );
ElError ElBasisPursuit_c
( ElConstMatrix_c A, ElConstMatrix_c b, ElMatrix_c z, ElInt* numIts );
ElError ElBasisPursuit_z
( ElConstMatrix_z A, ElConstMatrix_z b, ElMatrix_z z, ElInt* numIts );

ElError ElBasisPursuitDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s z, 
  ElInt* numIts );
ElError ElBasisPursuitDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d z, 
  ElInt* numIts );
ElError ElBasisPursuitDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c b, ElDistMatrix_c z, 
  ElInt* numIts );
ElError ElBasisPursuitDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z b, ElDistMatrix_z z, 
  ElInt* numIts );

/* TODO: Expert versions */

/* Least Absolute Shrinkage and Selection Operator (LASSO)
   ======================================================= */
ElError ElLasso_s
( ElConstMatrix_s A, ElConstMatrix_s b, float lambda, 
  ElMatrix_s z, ElInt* numIts );
ElError ElLasso_d
( ElConstMatrix_d A, ElConstMatrix_d b, double lambda, 
  ElMatrix_d z, ElInt* numIts );
ElError ElLasso_c
( ElConstMatrix_c A, ElConstMatrix_c b, float lambda, 
  ElMatrix_c z, ElInt* numIts );
ElError ElLasso_z
( ElConstMatrix_z A, ElConstMatrix_z b, double lambda, 
  ElMatrix_z z, ElInt* numIts );

ElError ElLassoDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, float lambda, 
  ElDistMatrix_s z, ElInt* numIts );
ElError ElLassoDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, double lambda, 
  ElDistMatrix_d z, ElInt* numIts );
ElError ElLassoDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c b, float lambda, 
  ElDistMatrix_c z, ElInt* numIts );
ElError ElLassoDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z b, double lambda, 
  ElDistMatrix_z z, ElInt* numIts );

/* Linear program
   ============== */
ElError ElLinearProgram_s
( ElConstMatrix_s A, ElConstMatrix_s b, ElConstMatrix_s c, 
  ElMatrix_s z, ElInt* numIts );
ElError ElLinearProgram_d
( ElConstMatrix_d A, ElConstMatrix_d b, ElConstMatrix_d c, 
  ElMatrix_d z, ElInt* numIts );

ElError ElLinearProgramDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElDistMatrix_s z, ElInt* numIts );
ElError ElLinearProgramDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElDistMatrix_d z, ElInt* numIts );

/* Logistic regression
   =================== */
typedef enum {
  EL_NO_PENALTY,
  EL_L1_PENALTY,
  EL_L2_PENALTY
} ElRegularization;

ElError ElLogisticRegression_s
( ElConstMatrix_s G, ElConstMatrix_s q, ElMatrix_s z, float gamma, 
  ElRegularization penalty, ElInt* numIts );
ElError ElLogisticRegression_d
( ElConstMatrix_d G, ElConstMatrix_d q, ElMatrix_d z, double gamma, 
  ElRegularization penalty, ElInt* numIts );

ElError ElLogisticRegressionDist_s
( ElConstDistMatrix_s G, ElConstDistMatrix_s q, ElDistMatrix_s z, float gamma, 
  ElRegularization penalty, ElInt* numIts );
ElError ElLogisticRegressionDist_d
( ElConstDistMatrix_d G, ElConstDistMatrix_d q, ElDistMatrix_d z, double gamma, 
  ElRegularization penalty, ElInt* numIts );

/* Model fit
   ========= */
ElError ElModelFit_s
( void (*lossProx)(ElMatrix_s,float), 
  void (*regProx)(ElMatrix_s,float),
  ElConstMatrix_s A, ElConstMatrix_s b, ElMatrix_s w, float rho, 
  ElInt* numIts );
ElError ElModelFit_d
( void (*lossProx)(ElMatrix_d,double), 
  void (*regProx)(ElMatrix_d,double),
  ElConstMatrix_d A, ElConstMatrix_d b, ElMatrix_d w, double rho,
  ElInt* numIts );

ElError ElModelFitDist_s
( void (*lossProx)(ElDistMatrix_s,float), 
  void (*regProx)(ElDistMatrix_s,float),
  ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s w, float rho,
  ElInt* numIts );
ElError ElModelFitDist_d
( void (*lossProx)(ElDistMatrix_d,double), 
  void (*regProx)(ElDistMatrix_d,double),
  ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d w, double rho,
  ElInt* numIts );

/* TODO: Expert versions */

/* Non-negative matrix factorization
   ================================= */
ElError ElNMF_s( ElConstMatrix_s A, ElMatrix_s X, ElMatrix_s Y );
ElError ElNMF_d( ElConstMatrix_d A, ElMatrix_d X, ElMatrix_d Y );

ElError ElNMFDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s X, ElDistMatrix_s Y );
ElError ElNMFDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d X, ElDistMatrix_d Y );

/* Non-negative least squares
   ========================== */
ElError ElNonNegativeLeastSquares_s
( ElConstMatrix_s A, ElConstMatrix_s Y, ElMatrix_s Z, ElInt* numIts );
ElError ElNonNegativeLeastSquares_d
( ElConstMatrix_d A, ElConstMatrix_d Y, ElMatrix_d Z, ElInt* numIts );

ElError ElNonNegativeLeastSquaresDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s Y, ElDistMatrix_s Z, 
  ElInt* numIts );
ElError ElNonNegativeLeastSquaresDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d Y, ElDistMatrix_d Z, 
  ElInt* numIts );

/* TODO: Expert versions */

/* Quadratic program
   ================= */
ElError ElQuadraticProgram_s
( ElConstMatrix_s P, ElConstMatrix_s S, float lb, float ub, 
  ElMatrix_s Z, ElInt* numIts );
ElError ElQuadraticProgram_d
( ElConstMatrix_d P, ElConstMatrix_d S, double lb, double ub, 
  ElMatrix_d Z, ElInt* numIts );

ElError ElQuadraticProgramDist_s
( ElConstDistMatrix_s P, ElConstDistMatrix_s S, float lb, float ub, 
  ElDistMatrix_s Z, ElInt* numIts );
ElError ElQuadraticProgramDist_d
( ElConstDistMatrix_d P, ElConstDistMatrix_d S, double lb, double ub, 
  ElDistMatrix_d Z, ElInt* numIts );

/* TODO: Expert versions */

/* Robust Principal Component Analysis
   =================================== */
ElError ElRPCA_s( ElConstMatrix_s M, ElMatrix_s L, ElMatrix_s S );
ElError ElRPCA_d( ElConstMatrix_d M, ElMatrix_d L, ElMatrix_d S );
ElError ElRPCA_c( ElConstMatrix_c M, ElMatrix_c L, ElMatrix_c S );
ElError ElRPCA_z( ElConstMatrix_z M, ElMatrix_z L, ElMatrix_z S );

ElError ElRPCADist_s
( ElConstDistMatrix_s M, ElDistMatrix_s L, ElDistMatrix_s S );
ElError ElRPCADist_d
( ElConstDistMatrix_d M, ElDistMatrix_d L, ElDistMatrix_d S );
ElError ElRPCADist_c
( ElConstDistMatrix_c M, ElDistMatrix_c L, ElDistMatrix_c S );
ElError ElRPCADist_z
( ElConstDistMatrix_z M, ElDistMatrix_z L, ElDistMatrix_z S );

/* Sparse inverse covariance selection
   =================================== */
ElError ElSparseInvCov_s
( ElConstMatrix_s D, float lambda, ElMatrix_s Z, ElInt* numIts );
ElError ElSparseInvCov_d
( ElConstMatrix_d D, double lambda, ElMatrix_d Z, ElInt* numIts );
ElError ElSparseInvCov_c
( ElConstMatrix_c D, float lambda, ElMatrix_c Z, ElInt* numIts );
ElError ElSparseInvCov_z
( ElConstMatrix_z D, double lambda, ElMatrix_z Z, ElInt* numIts );

ElError ElSparseInvCovDist_s
( ElConstDistMatrix_s D, float lambda, ElDistMatrix_s Z, ElInt* numIts );
ElError ElSparseInvCovDist_d
( ElConstDistMatrix_d D, double lambda, ElDistMatrix_d Z, ElInt* numIts );
ElError ElSparseInvCovDist_c
( ElConstDistMatrix_c D, float lambda, ElDistMatrix_c Z, ElInt* numIts );
ElError ElSparseInvCovDist_z
( ElConstDistMatrix_z D, double lambda, ElDistMatrix_z Z, ElInt* numIts );

/* TODO: Expert versions */

/* Support Vector Machine
   ====================== */
ElError ElSVM_s
( ElConstMatrix_s G, ElConstMatrix_s q, ElMatrix_s z, float gamma, 
  ElInt* numIts );
ElError ElSVM_d
( ElConstMatrix_d G, ElConstMatrix_d q, ElMatrix_d z, double gamma, 
  ElInt* numIts );

ElError ElSVMDist_s
( ElConstDistMatrix_s G, ElConstDistMatrix_s q, ElDistMatrix_s z, float gamma, 
  ElInt* numIts );
ElError ElSVMDist_d
( ElConstDistMatrix_d G, ElConstDistMatrix_d q, ElDistMatrix_d z, double gamma, 
  ElInt* numIts );

/* TODO: Expert versions */

/* Utilities
   ========= */

/* Clipping
   -------- */
ElError ElLowerClip_s( ElMatrix_s X, float lowerBound );
ElError ElLowerClip_d( ElMatrix_d X, double lowerBound );

ElError ElLowerClipDist_s( ElDistMatrix_s X, float lowerBound );
ElError ElLowerClipDist_d( ElDistMatrix_d X, double lowerBound );

ElError ElUpperClip_s( ElMatrix_s X, float upperBound );
ElError ElUpperClip_d( ElMatrix_d X, double upperBound );

ElError ElUpperClipDist_s( ElDistMatrix_s X, float upperBound );
ElError ElUpperClipDist_d( ElDistMatrix_d X, double upperBound );

ElError ElClip_s( ElMatrix_s X, float lowerBound, float upperBound );
ElError ElClip_d( ElMatrix_d X, double lowerBound, double upperBound );

ElError ElClipDist_s( ElDistMatrix_s X, float lowerBound, float upperBound );
ElError ElClipDist_d( ElDistMatrix_d X, double lowerBound, double upperBound );

/* Covariance
   ---------- */
ElError ElCovariance_s( ElConstMatrix_s D, ElMatrix_s S );
ElError ElCovariance_d( ElConstMatrix_d D, ElMatrix_d S );
ElError ElCovariance_c( ElConstMatrix_c D, ElMatrix_c S );
ElError ElCovariance_z( ElConstMatrix_z D, ElMatrix_z S );

ElError ElCovarianceDist_s( ElConstDistMatrix_s D, ElDistMatrix_s S );
ElError ElCovarianceDist_d( ElConstDistMatrix_d D, ElDistMatrix_d S );
ElError ElCovarianceDist_c( ElConstDistMatrix_c D, ElDistMatrix_c S );
ElError ElCovarianceDist_z( ElConstDistMatrix_z D, ElDistMatrix_z S );

/* Frobenius-norm proximal map
   --------------------------- */
/* arg min || A ||_F + rho/2 || A - A0 ||_F^2 
      A                                       */
ElError ElFrobeniusProx_s( ElMatrix_s A, float rho );
ElError ElFrobeniusProx_d( ElMatrix_d A, double rho );
ElError ElFrobeniusProx_c( ElMatrix_c A, float rho );
ElError ElFrobeniusProx_z( ElMatrix_z A, double rho );

ElError ElFrobeniusProxDist_s( ElDistMatrix_s A, float rho );
ElError ElFrobeniusProxDist_d( ElDistMatrix_d A, double rho );
ElError ElFrobeniusProxDist_c( ElDistMatrix_c A, float rho );
ElError ElFrobeniusProxDist_z( ElDistMatrix_z A, double rho );

/* Hinge-loss proximal map
   ----------------------- */
ElError ElHingeLossProx_s( ElMatrix_s A, float rho );
ElError ElHingeLossProx_d( ElMatrix_d A, double rho );

ElError ElHingeLossProxDist_s( ElDistMatrix_s A, float rho );
ElError ElHingeLossProxDist_d( ElDistMatrix_d A, double rho );

/* Log barrier
   ----------- */
ElError ElLogBarrier_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* barrier );
ElError ElLogBarrier_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* barrier );
ElError ElLogBarrier_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* barrier );
ElError ElLogBarrier_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* barrier );

ElError ElLogBarrierDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* barrier );
ElError ElLogBarrierDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* barrier );
ElError ElLogBarrierDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* barrier );
ElError ElLogBarrierDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* barrier );

/* TODO: Version which allows overwriting? */

/* Log-det divergence
   ------------------ */
ElError ElLogDetDiv_s
( ElUpperOrLower uplo, ElConstMatrix_s A, ElConstMatrix_s B, float* div );
ElError ElLogDetDiv_d
( ElUpperOrLower uplo, ElConstMatrix_d A, ElConstMatrix_d B, double* div );
ElError ElLogDetDiv_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElConstMatrix_c B, float* div );
ElError ElLogDetDiv_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElConstMatrix_z B, double* div );

ElError ElLogDetDivDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, ElConstDistMatrix_s B, 
  float* div );
ElError ElLogDetDivDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, ElConstDistMatrix_d B, 
  double* div );
ElError ElLogDetDivDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElConstDistMatrix_c B, 
  float* div );
ElError ElLogDetDivDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElConstDistMatrix_z B, 
  double* div );

/* Logistic proximal map
   --------------------- */
ElError ElLogisticProx_s( ElMatrix_s A, float rho );
ElError ElLogisticProx_d( ElMatrix_d A, double rho );

ElError ElLogisticProxDist_s( ElDistMatrix_s A, float rho );
ElError ElLogisticProxDist_d( ElDistMatrix_d A, double rho );

/* Singular-value soft thresholding
   -------------------------------- */
ElError ElSVT_s( ElMatrix_s A, float rho, bool relative );
ElError ElSVT_d( ElMatrix_d A, double rho, bool relative );
ElError ElSVT_c( ElMatrix_c A, float rho, bool relative );
ElError ElSVT_z( ElMatrix_z A, double rho, bool relative );

ElError ElSVTDist_s( ElDistMatrix_s A, float rho, bool relative );
ElError ElSVTDist_d( ElDistMatrix_d A, double rho, bool relative );
ElError ElSVTDist_c( ElDistMatrix_c A, float rho, bool relative );
ElError ElSVTDist_z( ElDistMatrix_z A, double rho, bool relative );

/* TODO: Expert versions */

/* Soft-thresholding
   ----------------- */
ElError ElSoftThreshold_s( ElMatrix_s A, float rho, bool relative );
ElError ElSoftThreshold_d( ElMatrix_d A, double rho, bool relative );
ElError ElSoftThreshold_c( ElMatrix_c A, float rho, bool relative );
ElError ElSoftThreshold_z( ElMatrix_z A, double rho, bool relative );

ElError ElSoftThresholdDist_s( ElDistMatrix_s A, float rho, bool relative );
ElError ElSoftThresholdDist_d( ElDistMatrix_d A, double rho, bool relative );
ElError ElSoftThresholdDist_c( ElDistMatrix_c A, float rho, bool relative );
ElError ElSoftThresholdDist_z( ElDistMatrix_z A, double rho, bool relative );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_OPTIMIZATION_C_H */
