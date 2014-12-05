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

EL_EXPORT ElError ElBasisPursuit_s
( ElConstMatrix_s A, ElConstMatrix_s b, ElMatrix_s z, ElInt* numIts );
EL_EXPORT ElError ElBasisPursuit_d
( ElConstMatrix_d A, ElConstMatrix_d b, ElMatrix_d z, ElInt* numIts );
EL_EXPORT ElError ElBasisPursuit_c
( ElConstMatrix_c A, ElConstMatrix_c b, ElMatrix_c z, ElInt* numIts );
EL_EXPORT ElError ElBasisPursuit_z
( ElConstMatrix_z A, ElConstMatrix_z b, ElMatrix_z z, ElInt* numIts );

EL_EXPORT ElError ElBasisPursuitDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s z, 
  ElInt* numIts );
EL_EXPORT ElError ElBasisPursuitDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d z, 
  ElInt* numIts );
EL_EXPORT ElError ElBasisPursuitDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c b, ElDistMatrix_c z, 
  ElInt* numIts );
EL_EXPORT ElError ElBasisPursuitDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z b, ElDistMatrix_z z, 
  ElInt* numIts );

/* TODO: Expert versions */

/* Least Absolute Shrinkage and Selection Operator (LASSO)
   ======================================================= */
EL_EXPORT ElError ElLasso_s
( ElConstMatrix_s A, ElConstMatrix_s b, float lambda, 
  ElMatrix_s z, ElInt* numIts );
EL_EXPORT ElError ElLasso_d
( ElConstMatrix_d A, ElConstMatrix_d b, double lambda, 
  ElMatrix_d z, ElInt* numIts );
EL_EXPORT ElError ElLasso_c
( ElConstMatrix_c A, ElConstMatrix_c b, float lambda, 
  ElMatrix_c z, ElInt* numIts );
EL_EXPORT ElError ElLasso_z
( ElConstMatrix_z A, ElConstMatrix_z b, double lambda, 
  ElMatrix_z z, ElInt* numIts );

EL_EXPORT ElError ElLassoDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, float lambda, 
  ElDistMatrix_s z, ElInt* numIts );
EL_EXPORT ElError ElLassoDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, double lambda, 
  ElDistMatrix_d z, ElInt* numIts );
EL_EXPORT ElError ElLassoDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c b, float lambda, 
  ElDistMatrix_c z, ElInt* numIts );
EL_EXPORT ElError ElLassoDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z b, double lambda, 
  ElDistMatrix_z z, ElInt* numIts );

/* Linear program
   ============== */
EL_EXPORT ElError ElLinearProgram_s
( ElConstMatrix_s A, ElConstMatrix_s b, ElConstMatrix_s c, ElMatrix_s x );
EL_EXPORT ElError ElLinearProgram_d
( ElConstMatrix_d A, ElConstMatrix_d b, ElConstMatrix_d c, ElMatrix_d x );

EL_EXPORT ElError ElLinearProgramDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElDistMatrix_s x );
EL_EXPORT ElError ElLinearProgramDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElDistMatrix_d x );

EL_EXPORT ElError ElLinearProgramSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s b, ElConstMatrix_s c, ElMatrix_s x );
EL_EXPORT ElError ElLinearProgramSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d b, ElConstMatrix_d c, ElMatrix_d x );

EL_EXPORT ElError ElLinearProgramDistSparse_s
( ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s b, ElConstDistMultiVec_s c, ElDistMultiVec_s x );
EL_EXPORT ElError ElLinearProgramDistSparse_d
( ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d b, ElConstDistMultiVec_d c, ElDistMultiVec_d x );

/* Infeasible Path-following IPM
   ----------------------------- */
EL_EXPORT ElError ElLinearProgramIPF_s
( ElConstMatrix_s A,
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s s, ElMatrix_s x, ElMatrix_s l );
EL_EXPORT ElError ElLinearProgramIPF_d
( ElConstMatrix_d A,
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d s, ElMatrix_d x, ElMatrix_d l );
EL_EXPORT ElError ElLinearProgramIPFSparse_s
( ElConstSparseMatrix_s A,
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s s, ElMatrix_s x, ElMatrix_s l );
EL_EXPORT ElError ElLinearProgramIPFSparse_d
( ElConstSparseMatrix_d A,
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d s, ElMatrix_d x, ElMatrix_d l );

EL_EXPORT ElError ElLinearProgramIPFDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b, ElConstDistMatrix_s c,
  ElDistMatrix_s s, ElDistMatrix_s x, ElDistMatrix_s l );
EL_EXPORT ElError ElLinearProgramIPFDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b, ElConstDistMatrix_d c,
  ElDistMatrix_d s, ElDistMatrix_d x, ElDistMatrix_d l );
EL_EXPORT ElError ElLinearProgramIPFDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b, ElConstDistMultiVec_s c,
  ElDistMultiVec_s s, ElDistMultiVec_s x, ElDistMultiVec_s l );
EL_EXPORT ElError ElLinearProgramIPFDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b, ElConstDistMultiVec_d c,
  ElDistMultiVec_d s, ElDistMultiVec_d x, ElDistMultiVec_d l );

/* Mehotra's Predictor-Corrector IPM
   --------------------------------- */
EL_EXPORT ElError ElLinearProgramMehrotra_s
( ElConstMatrix_s A,
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s s, ElMatrix_s x, ElMatrix_s l );
EL_EXPORT ElError ElLinearProgramMehrotra_d
( ElConstMatrix_d A,
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d s, ElMatrix_d x, ElMatrix_d l );

EL_EXPORT ElError ElLinearProgramMehrotraDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b, ElConstDistMatrix_s c,
  ElDistMatrix_s s, ElDistMatrix_s x, ElDistMatrix_s l );
EL_EXPORT ElError ElLinearProgramMehrotraDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b, ElConstDistMatrix_d c,
  ElDistMatrix_d s, ElDistMatrix_d x, ElDistMatrix_d l );

EL_EXPORT ElError ElLinearProgramMehrotraSparse_s
( ElConstSparseMatrix_s A,
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s s, ElMatrix_s x, ElMatrix_s l );
EL_EXPORT ElError ElLinearProgramMehrotraSparse_d
( ElConstSparseMatrix_d A,
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d s, ElMatrix_d x, ElMatrix_d l );

EL_EXPORT ElError ElLinearProgramMehrotraDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b, ElConstDistMultiVec_s c,
  ElDistMultiVec_s s, ElDistMultiVec_s x, ElDistMultiVec_s l );
EL_EXPORT ElError ElLinearProgramMehrotraDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b, ElConstDistMultiVec_d c,
  ElDistMultiVec_d s, ElDistMultiVec_d x, ElDistMultiVec_d l );

/* ADMM
   ---- */
EL_EXPORT ElError ElLinearProgramADMM_s
( ElConstMatrix_s A, ElConstMatrix_s b, ElConstMatrix_s c, 
  ElMatrix_s z, ElInt* numIts );
EL_EXPORT ElError ElLinearProgramADMM_d
( ElConstMatrix_d A, ElConstMatrix_d b, ElConstMatrix_d c, 
  ElMatrix_d z, ElInt* numIts );

EL_EXPORT ElError ElLinearProgramADMMDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElDistMatrix_s z, ElInt* numIts );
EL_EXPORT ElError ElLinearProgramADMMDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElDistMatrix_d z, ElInt* numIts );

/* Logistic regression
   =================== */
typedef enum {
  EL_NO_PENALTY,
  EL_L1_PENALTY,
  EL_L2_PENALTY
} ElRegularization;

EL_EXPORT ElError ElLogisticRegression_s
( ElConstMatrix_s G, ElConstMatrix_s q, ElMatrix_s z, float gamma, 
  ElRegularization penalty, ElInt* numIts );
EL_EXPORT ElError ElLogisticRegression_d
( ElConstMatrix_d G, ElConstMatrix_d q, ElMatrix_d z, double gamma, 
  ElRegularization penalty, ElInt* numIts );

EL_EXPORT ElError ElLogisticRegressionDist_s
( ElConstDistMatrix_s G, ElConstDistMatrix_s q, ElDistMatrix_s z, float gamma, 
  ElRegularization penalty, ElInt* numIts );
EL_EXPORT ElError ElLogisticRegressionDist_d
( ElConstDistMatrix_d G, ElConstDistMatrix_d q, ElDistMatrix_d z, double gamma, 
  ElRegularization penalty, ElInt* numIts );

/* Model fit
   ========= */
EL_EXPORT ElError ElModelFit_s
( void (*lossProx)(ElMatrix_s,float), 
  void (*regProx)(ElMatrix_s,float),
  ElConstMatrix_s A, ElConstMatrix_s b, ElMatrix_s w, float rho, 
  ElInt* numIts );
EL_EXPORT ElError ElModelFit_d
( void (*lossProx)(ElMatrix_d,double), 
  void (*regProx)(ElMatrix_d,double),
  ElConstMatrix_d A, ElConstMatrix_d b, ElMatrix_d w, double rho,
  ElInt* numIts );

EL_EXPORT ElError ElModelFitDist_s
( void (*lossProx)(ElDistMatrix_s,float), 
  void (*regProx)(ElDistMatrix_s,float),
  ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s w, float rho,
  ElInt* numIts );
EL_EXPORT ElError ElModelFitDist_d
( void (*lossProx)(ElDistMatrix_d,double), 
  void (*regProx)(ElDistMatrix_d,double),
  ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d w, double rho,
  ElInt* numIts );

/* TODO: Expert versions */

/* Non-negative matrix factorization
   ================================= */
EL_EXPORT ElError ElNMF_s( ElConstMatrix_s A, ElMatrix_s X, ElMatrix_s Y );
EL_EXPORT ElError ElNMF_d( ElConstMatrix_d A, ElMatrix_d X, ElMatrix_d Y );

EL_EXPORT ElError ElNMFDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s X, ElDistMatrix_s Y );
EL_EXPORT ElError ElNMFDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d X, ElDistMatrix_d Y );

/* Non-negative least squares
   ========================== */
EL_EXPORT ElError ElNonNegativeLeastSquares_s
( ElConstMatrix_s A, ElConstMatrix_s Y, ElMatrix_s Z, ElInt* numIts );
EL_EXPORT ElError ElNonNegativeLeastSquares_d
( ElConstMatrix_d A, ElConstMatrix_d Y, ElMatrix_d Z, ElInt* numIts );

EL_EXPORT ElError ElNonNegativeLeastSquaresDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s Y, ElDistMatrix_s Z, 
  ElInt* numIts );
EL_EXPORT ElError ElNonNegativeLeastSquaresDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d Y, ElDistMatrix_d Z, 
  ElInt* numIts );

/* TODO: Expert versions */

/* Quadratic program
   ================= */
EL_EXPORT ElError ElQuadraticProgram_s
( ElConstMatrix_s Q, ElConstMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s x );
EL_EXPORT ElError ElQuadraticProgram_d
( ElConstMatrix_d Q, ElConstMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d x );

EL_EXPORT ElError ElQuadraticProgramDist_s
( ElConstDistMatrix_s Q, ElConstDistMatrix_s A, 
  ElConstDistMatrix_s b, ElConstDistMatrix_s c,
  ElDistMatrix_s x );
EL_EXPORT ElError ElQuadraticProgramDist_d
( ElConstDistMatrix_d Q, ElConstDistMatrix_d A, 
  ElConstDistMatrix_d b, ElConstDistMatrix_d c,
  ElDistMatrix_d x );

EL_EXPORT ElError ElQuadraticProgramSparse_s
( ElConstSparseMatrix_s Q, ElConstSparseMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s x );
EL_EXPORT ElError ElQuadraticProgramSparse_d
( ElConstSparseMatrix_d Q, ElConstSparseMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d x );

EL_EXPORT ElError ElQuadraticProgramDistSparse_s
( ElConstDistSparseMatrix_s Q, ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s b, ElConstDistMultiVec_s c,
  ElDistMultiVec_s x );
EL_EXPORT ElError ElQuadraticProgramDistSparse_d
( ElConstDistSparseMatrix_d Q, ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d b, ElConstDistMultiVec_d c,
  ElDistMultiVec_d x );

/* IPF
   --- */
EL_EXPORT ElError ElQuadraticProgramIPF_s
( ElConstMatrix_s Q, ElConstMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s s, ElMatrix_s x, ElMatrix_s l );
EL_EXPORT ElError ElQuadraticProgramIPF_d
( ElConstMatrix_d Q, ElConstMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d s, ElMatrix_d x, ElMatrix_d l );

EL_EXPORT ElError ElQuadraticProgramIPFDist_s
( ElConstDistMatrix_s Q, ElConstDistMatrix_s A, 
  ElConstDistMatrix_s b, ElConstDistMatrix_s c,
  ElDistMatrix_s s, ElDistMatrix_s x, ElDistMatrix_s l );
EL_EXPORT ElError ElQuadraticProgramIPFDist_d
( ElConstDistMatrix_d Q, ElConstDistMatrix_d A, 
  ElConstDistMatrix_d b, ElConstDistMatrix_d c,
  ElDistMatrix_d s, ElDistMatrix_d x, ElDistMatrix_d l );

EL_EXPORT ElError ElQuadraticProgramIPFSparse_s
( ElConstSparseMatrix_s Q, ElConstSparseMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s s, ElMatrix_s x, ElMatrix_s l );
EL_EXPORT ElError ElQuadraticProgramIPFSparse_d
( ElConstSparseMatrix_d Q, ElConstSparseMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d s, ElMatrix_d x, ElMatrix_d l );

EL_EXPORT ElError ElQuadraticProgramIPFDistSparse_s
( ElConstDistSparseMatrix_s Q, ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s b, ElConstDistMultiVec_s c,
  ElDistMultiVec_s s, ElDistMultiVec_s x, ElDistMultiVec_s l );
EL_EXPORT ElError ElQuadraticProgramIPFDistSparse_d
( ElConstDistSparseMatrix_d Q, ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d b, ElConstDistMultiVec_d c,
  ElDistMultiVec_d s, ElDistMultiVec_d x, ElDistMultiVec_d l );

/* Mehrotra
   -------- */
EL_EXPORT ElError ElQuadraticProgramMehrotra_s
( ElConstMatrix_s Q, ElConstMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s s, ElMatrix_s x, ElMatrix_s l );
EL_EXPORT ElError ElQuadraticProgramMehrotra_d
( ElConstMatrix_d Q, ElConstMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d s, ElMatrix_d x, ElMatrix_d l );

EL_EXPORT ElError ElQuadraticProgramMehrotraDist_s
( ElConstDistMatrix_s Q, ElConstDistMatrix_s A, 
  ElConstDistMatrix_s b, ElConstDistMatrix_s c,
  ElDistMatrix_s s, ElDistMatrix_s x, ElDistMatrix_s l );
EL_EXPORT ElError ElQuadraticProgramMehrotraDist_d
( ElConstDistMatrix_d Q, ElConstDistMatrix_d A, 
  ElConstDistMatrix_d b, ElConstDistMatrix_d c,
  ElDistMatrix_d s, ElDistMatrix_d x, ElDistMatrix_d l );

EL_EXPORT ElError ElQuadraticProgramMehrotraSparse_s
( ElConstSparseMatrix_s Q, ElConstSparseMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s s, ElMatrix_s x, ElMatrix_s l );
EL_EXPORT ElError ElQuadraticProgramMehrotraSparse_d
( ElConstSparseMatrix_d Q, ElConstSparseMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d s, ElMatrix_d x, ElMatrix_d l );

EL_EXPORT ElError ElQuadraticProgramMehrotraDistSparse_s
( ElConstDistSparseMatrix_s Q, ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s b, ElConstDistMultiVec_s c,
  ElDistMultiVec_s s, ElDistMultiVec_s x, ElDistMultiVec_s l );
EL_EXPORT ElError ElQuadraticProgramMehrotraDistSparse_d
( ElConstDistSparseMatrix_d Q, ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d b, ElConstDistMultiVec_d c,
  ElDistMultiVec_d s, ElDistMultiVec_d x, ElDistMultiVec_d l );

/* ADMM (non-conic form)
   --------------------- */
EL_EXPORT ElError ElQuadraticProgramADMM_s
( ElConstMatrix_s Q, ElConstMatrix_s C, float lb, float ub, 
  ElMatrix_s Z, ElInt* numIts );
EL_EXPORT ElError ElQuadraticProgramADMM_d
( ElConstMatrix_d Q, ElConstMatrix_d C, double lb, double ub, 
  ElMatrix_d Z, ElInt* numIts );

EL_EXPORT ElError ElQuadraticProgramADMMDist_s
( ElConstDistMatrix_s Q, ElConstDistMatrix_s C, float lb, float ub, 
  ElDistMatrix_s Z, ElInt* numIts );
EL_EXPORT ElError ElQuadraticProgramADMMDist_d
( ElConstDistMatrix_d Q, ElConstDistMatrix_d C, double lb, double ub, 
  ElDistMatrix_d Z, ElInt* numIts );

/* TODO: Expert versions */

/* Robust Principal Component Analysis
   =================================== */
EL_EXPORT ElError ElRPCA_s( ElConstMatrix_s M, ElMatrix_s L, ElMatrix_s S );
EL_EXPORT ElError ElRPCA_d( ElConstMatrix_d M, ElMatrix_d L, ElMatrix_d S );
EL_EXPORT ElError ElRPCA_c( ElConstMatrix_c M, ElMatrix_c L, ElMatrix_c S );
EL_EXPORT ElError ElRPCA_z( ElConstMatrix_z M, ElMatrix_z L, ElMatrix_z S );

EL_EXPORT ElError ElRPCADist_s
( ElConstDistMatrix_s M, ElDistMatrix_s L, ElDistMatrix_s S );
EL_EXPORT ElError ElRPCADist_d
( ElConstDistMatrix_d M, ElDistMatrix_d L, ElDistMatrix_d S );
EL_EXPORT ElError ElRPCADist_c
( ElConstDistMatrix_c M, ElDistMatrix_c L, ElDistMatrix_c S );
EL_EXPORT ElError ElRPCADist_z
( ElConstDistMatrix_z M, ElDistMatrix_z L, ElDistMatrix_z S );

/* Sparse inverse covariance selection
   =================================== */
EL_EXPORT ElError ElSparseInvCov_s
( ElConstMatrix_s D, float lambda, ElMatrix_s Z, ElInt* numIts );
EL_EXPORT ElError ElSparseInvCov_d
( ElConstMatrix_d D, double lambda, ElMatrix_d Z, ElInt* numIts );
EL_EXPORT ElError ElSparseInvCov_c
( ElConstMatrix_c D, float lambda, ElMatrix_c Z, ElInt* numIts );
EL_EXPORT ElError ElSparseInvCov_z
( ElConstMatrix_z D, double lambda, ElMatrix_z Z, ElInt* numIts );

EL_EXPORT ElError ElSparseInvCovDist_s
( ElConstDistMatrix_s D, float lambda, ElDistMatrix_s Z, ElInt* numIts );
EL_EXPORT ElError ElSparseInvCovDist_d
( ElConstDistMatrix_d D, double lambda, ElDistMatrix_d Z, ElInt* numIts );
EL_EXPORT ElError ElSparseInvCovDist_c
( ElConstDistMatrix_c D, float lambda, ElDistMatrix_c Z, ElInt* numIts );
EL_EXPORT ElError ElSparseInvCovDist_z
( ElConstDistMatrix_z D, double lambda, ElDistMatrix_z Z, ElInt* numIts );

/* TODO: Expert versions */

/* Support Vector Machine
   ====================== */
EL_EXPORT ElError ElSVM_s
( ElConstMatrix_s G, ElConstMatrix_s q, ElMatrix_s z, float gamma, 
  ElInt* numIts );
EL_EXPORT ElError ElSVM_d
( ElConstMatrix_d G, ElConstMatrix_d q, ElMatrix_d z, double gamma, 
  ElInt* numIts );

EL_EXPORT ElError ElSVMDist_s
( ElConstDistMatrix_s G, ElConstDistMatrix_s q, ElDistMatrix_s z, float gamma, 
  ElInt* numIts );
EL_EXPORT ElError ElSVMDist_d
( ElConstDistMatrix_d G, ElConstDistMatrix_d q, ElDistMatrix_d z, double gamma, 
  ElInt* numIts );

/* TODO: Expert versions */

/* Utilities
   ========= */

/* Clipping
   -------- */
EL_EXPORT ElError ElLowerClip_s( ElMatrix_s X, float lowerBound );
EL_EXPORT ElError ElLowerClip_d( ElMatrix_d X, double lowerBound );

EL_EXPORT ElError ElLowerClipDist_s( ElDistMatrix_s X, float lowerBound );
EL_EXPORT ElError ElLowerClipDist_d( ElDistMatrix_d X, double lowerBound );

EL_EXPORT ElError ElUpperClip_s( ElMatrix_s X, float upperBound );
EL_EXPORT ElError ElUpperClip_d( ElMatrix_d X, double upperBound );

EL_EXPORT ElError ElUpperClipDist_s( ElDistMatrix_s X, float upperBound );
EL_EXPORT ElError ElUpperClipDist_d( ElDistMatrix_d X, double upperBound );

EL_EXPORT ElError ElClip_s
( ElMatrix_s X, float lowerBound, float upperBound );
EL_EXPORT ElError ElClip_d
( ElMatrix_d X, double lowerBound, double upperBound );

EL_EXPORT ElError ElClipDist_s
( ElDistMatrix_s X, float lowerBound, float upperBound );
EL_EXPORT ElError ElClipDist_d
( ElDistMatrix_d X, double lowerBound, double upperBound );

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

/* Frobenius-norm proximal map
   --------------------------- */
/* arg min || A ||_F + rho/2 || A - A0 ||_F^2 
      A                                       */
EL_EXPORT ElError ElFrobeniusProx_s( ElMatrix_s A, float rho );
EL_EXPORT ElError ElFrobeniusProx_d( ElMatrix_d A, double rho );
EL_EXPORT ElError ElFrobeniusProx_c( ElMatrix_c A, float rho );
EL_EXPORT ElError ElFrobeniusProx_z( ElMatrix_z A, double rho );

EL_EXPORT ElError ElFrobeniusProxDist_s( ElDistMatrix_s A, float rho );
EL_EXPORT ElError ElFrobeniusProxDist_d( ElDistMatrix_d A, double rho );
EL_EXPORT ElError ElFrobeniusProxDist_c( ElDistMatrix_c A, float rho );
EL_EXPORT ElError ElFrobeniusProxDist_z( ElDistMatrix_z A, double rho );

/* Hinge-loss proximal map
   ----------------------- */
EL_EXPORT ElError ElHingeLossProx_s( ElMatrix_s A, float rho );
EL_EXPORT ElError ElHingeLossProx_d( ElMatrix_d A, double rho );

EL_EXPORT ElError ElHingeLossProxDist_s( ElDistMatrix_s A, float rho );
EL_EXPORT ElError ElHingeLossProxDist_d( ElDistMatrix_d A, double rho );

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

/* Logistic proximal map
   --------------------- */
EL_EXPORT ElError ElLogisticProx_s( ElMatrix_s A, float rho );
EL_EXPORT ElError ElLogisticProx_d( ElMatrix_d A, double rho );

EL_EXPORT ElError ElLogisticProxDist_s( ElDistMatrix_s A, float rho );
EL_EXPORT ElError ElLogisticProxDist_d( ElDistMatrix_d A, double rho );

/* Singular-value soft thresholding
   -------------------------------- */
EL_EXPORT ElError ElSVT_s( ElMatrix_s A, float rho, bool relative );
EL_EXPORT ElError ElSVT_d( ElMatrix_d A, double rho, bool relative );
EL_EXPORT ElError ElSVT_c( ElMatrix_c A, float rho, bool relative );
EL_EXPORT ElError ElSVT_z( ElMatrix_z A, double rho, bool relative );

EL_EXPORT ElError ElSVTDist_s( ElDistMatrix_s A, float rho, bool relative );
EL_EXPORT ElError ElSVTDist_d( ElDistMatrix_d A, double rho, bool relative );
EL_EXPORT ElError ElSVTDist_c( ElDistMatrix_c A, float rho, bool relative );
EL_EXPORT ElError ElSVTDist_z( ElDistMatrix_z A, double rho, bool relative );

/* TODO: Expert versions */

/* Soft-thresholding
   ----------------- */
EL_EXPORT ElError ElSoftThreshold_s( ElMatrix_s A, float rho, bool relative );
EL_EXPORT ElError ElSoftThreshold_d( ElMatrix_d A, double rho, bool relative );
EL_EXPORT ElError ElSoftThreshold_c( ElMatrix_c A, float rho, bool relative );
EL_EXPORT ElError ElSoftThreshold_z( ElMatrix_z A, double rho, bool relative );

EL_EXPORT ElError ElSoftThresholdDist_s
( ElDistMatrix_s A, float rho, bool relative );
EL_EXPORT ElError ElSoftThresholdDist_d
( ElDistMatrix_d A, double rho, bool relative );
EL_EXPORT ElError ElSoftThresholdDist_c
( ElDistMatrix_c A, float rho, bool relative );
EL_EXPORT ElError ElSoftThresholdDist_z
( ElDistMatrix_z A, double rho, bool relative );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_OPTIMIZATION_C_H */
