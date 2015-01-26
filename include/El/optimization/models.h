/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_OPTIMIZATION_MODELS_C_H
#define EL_OPTIMIZATION_MODELS_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Basis pursuit
   ============= */
EL_EXPORT ElError ElBP_s
( ElConstMatrix_s A, ElConstMatrix_s b, ElMatrix_s x );
EL_EXPORT ElError ElBP_d
( ElConstMatrix_d A, ElConstMatrix_d b, ElMatrix_d x );

EL_EXPORT ElError ElBPDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s x );
EL_EXPORT ElError ElBPDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d x );

EL_EXPORT ElError ElBPSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s b, ElMatrix_s x );
EL_EXPORT ElError ElBPSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d b, ElMatrix_d x );

EL_EXPORT ElError ElBPDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s b, ElDistMultiVec_s x );
EL_EXPORT ElError ElBPDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d b, ElDistMultiVec_d x );

/* Expert verions
   -------------- */
EL_EXPORT ElError ElBPX_s
( ElConstMatrix_s A, ElConstMatrix_s b, ElMatrix_s x,
  ElLPDirectCtrl_s ctrl );
EL_EXPORT ElError ElBPX_d
( ElConstMatrix_d A, ElConstMatrix_d b, ElMatrix_d x,
  ElLPDirectCtrl_d ctrl );

EL_EXPORT ElError ElBPXDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s x,
  ElLPDirectCtrl_s ctrl );
EL_EXPORT ElError ElBPXDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d x,
  ElLPDirectCtrl_d ctrl );

EL_EXPORT ElError ElBPXSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s b, ElMatrix_s x,
  ElLPDirectCtrl_s ctrl );
EL_EXPORT ElError ElBPXSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d b, ElMatrix_d x,
  ElLPDirectCtrl_d ctrl );

EL_EXPORT ElError ElBPXDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s b, ElDistMultiVec_s x,
  ElLPDirectCtrl_s ctrl );
EL_EXPORT ElError ElBPXDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d b, ElDistMultiVec_d x,
  ElLPDirectCtrl_d ctrl );

/* ADMM
   ---- */
EL_EXPORT ElError ElBPADMM_s
( ElConstMatrix_s A, ElConstMatrix_s b, ElMatrix_s z, ElInt* numIts );
EL_EXPORT ElError ElBPADMM_d
( ElConstMatrix_d A, ElConstMatrix_d b, ElMatrix_d z, ElInt* numIts );
EL_EXPORT ElError ElBPADMM_c
( ElConstMatrix_c A, ElConstMatrix_c b, ElMatrix_c z, ElInt* numIts );
EL_EXPORT ElError ElBPADMM_z
( ElConstMatrix_z A, ElConstMatrix_z b, ElMatrix_z z, ElInt* numIts );

EL_EXPORT ElError ElBPADMMDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s z, 
  ElInt* numIts );
EL_EXPORT ElError ElBPADMMDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d z, 
  ElInt* numIts );
EL_EXPORT ElError ElBPADMMDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c b, ElDistMatrix_c z, 
  ElInt* numIts );
EL_EXPORT ElError ElBPADMMDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z b, ElDistMatrix_z z, 
  ElInt* numIts );

/* TODO: Expert versions */

/* Chebyshev point
   =============== */
EL_EXPORT ElError ElCP_s
( ElConstMatrix_s A, ElConstMatrix_s b, ElMatrix_s x );
EL_EXPORT ElError ElCP_d
( ElConstMatrix_d A, ElConstMatrix_d b, ElMatrix_d x );

EL_EXPORT ElError ElCPDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s x );
EL_EXPORT ElError ElCPDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d x );

EL_EXPORT ElError ElCPSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s b, ElMatrix_s x );
EL_EXPORT ElError ElCPSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d b, ElMatrix_d x );

EL_EXPORT ElError ElCPDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s b, ElDistMultiVec_s x );
EL_EXPORT ElError ElCPDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d b, ElDistMultiVec_d x );

/* Expert verions
   -------------- */
EL_EXPORT ElError ElCPX_s
( ElConstMatrix_s A, ElConstMatrix_s b, ElMatrix_s x,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElCPX_d
( ElConstMatrix_d A, ElConstMatrix_d b, ElMatrix_d x,
  ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElCPXDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s x,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElCPXDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d x,
  ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElCPXSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s b, ElMatrix_s x,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElCPXSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d b, ElMatrix_d x,
  ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElCPXDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s b, ElDistMultiVec_s x,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElCPXDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d b, ElDistMultiVec_d x,
  ElLPAffineCtrl_d ctrl );

/* Dantzig selector
   ================ */
EL_EXPORT ElError ElDS_s
( ElConstMatrix_s A, ElConstMatrix_s b, 
  float lambda, ElMatrix_s x );
EL_EXPORT ElError ElDS_d
( ElConstMatrix_d A, ElConstMatrix_d b, 
  double lambda, ElMatrix_d x );

EL_EXPORT ElError ElDSDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, 
  float lambda, ElDistMatrix_s x );
EL_EXPORT ElError ElDSDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, 
  double lambda, ElDistMatrix_d x );

EL_EXPORT ElError ElDSSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s b, 
  float lambda, ElMatrix_s x );
EL_EXPORT ElError ElDSSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d b, 
  double lambda, ElMatrix_d x );

EL_EXPORT ElError ElDSDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s b, 
  float lambda, ElDistMultiVec_s x );
EL_EXPORT ElError ElDSDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d b, 
  double lambda, ElDistMultiVec_d x );

/* Expert verions
   -------------- */
EL_EXPORT ElError ElDSX_s
( ElConstMatrix_s A, ElConstMatrix_s b, 
  float lambda, ElMatrix_s x, ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElDSX_d
( ElConstMatrix_d A, ElConstMatrix_d b, 
  double lambda, ElMatrix_d x, ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElDSXDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, 
  float lambda, ElDistMatrix_s x, ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElDSXDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, 
  double lambda, ElDistMatrix_d x, ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElDSXSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s b, 
  float lambda, ElMatrix_s x, ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElDSXSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d b, 
  double lambda, ElMatrix_d x, ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElDSXDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s b, 
  float lambda, ElDistMultiVec_s x, ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElDSXDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d b, 
  double lambda, ElDistMultiVec_d x, ElLPAffineCtrl_d ctrl );

/* Least Absolute Value regression
   =============================== */
EL_EXPORT ElError ElLAV_s
( ElConstMatrix_s A, ElConstMatrix_s b, ElMatrix_s x );
EL_EXPORT ElError ElLAV_d
( ElConstMatrix_d A, ElConstMatrix_d b, ElMatrix_d x );

EL_EXPORT ElError ElLAVDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s x );
EL_EXPORT ElError ElLAVDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d x );

EL_EXPORT ElError ElLAVSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s b, ElMatrix_s x );
EL_EXPORT ElError ElLAVSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d b, ElMatrix_d x );

EL_EXPORT ElError ElLAVDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s b, ElDistMultiVec_s x );
EL_EXPORT ElError ElLAVDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d b, ElDistMultiVec_d x );

/* Expert verions
   -------------- */
EL_EXPORT ElError ElLAVX_s
( ElConstMatrix_s A, ElConstMatrix_s b, ElMatrix_s x,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElLAVX_d
( ElConstMatrix_d A, ElConstMatrix_d b, ElMatrix_d x,
  ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElLAVXDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s x,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElLAVXDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d x,
  ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElLAVXSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s b, ElMatrix_s x,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElLAVXSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d b, ElMatrix_d x,
  ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElLAVXDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s b, ElDistMultiVec_s x,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElLAVXDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d b, ElDistMultiVec_d x,
  ElLPAffineCtrl_d ctrl );

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

/* Expert versions */
EL_EXPORT ElError ElNMFX_s
( ElConstMatrix_s A, ElMatrix_s X, ElMatrix_s Y, 
  ElQPDirectCtrl_s ctrl );
EL_EXPORT ElError ElNMFX_d
( ElConstMatrix_d A, ElMatrix_d X, ElMatrix_d Y, 
  ElQPDirectCtrl_d ctrl );

EL_EXPORT ElError ElNMFXDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s X, ElDistMatrix_s Y,
  ElQPDirectCtrl_s ctrl );
EL_EXPORT ElError ElNMFXDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d X, ElDistMatrix_d Y,
  ElQPDirectCtrl_d ctrl );

/* Non-negative least squares
   ========================== */
EL_EXPORT ElError ElNNLS_s
( ElConstMatrix_s A, ElConstMatrix_s b, ElMatrix_s x );
EL_EXPORT ElError ElNNLS_d
( ElConstMatrix_d A, ElConstMatrix_d b, ElMatrix_d x );

EL_EXPORT ElError ElNNLSDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s x );
EL_EXPORT ElError ElNNLSDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d x );

EL_EXPORT ElError ElNNLSDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s x );
EL_EXPORT ElError ElNNLSDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d x );

EL_EXPORT ElError ElNNLSDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s b, ElDistMultiVec_s x );
EL_EXPORT ElError ElNNLSDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d b, ElDistMultiVec_d x );

/* Expert versions */
EL_EXPORT ElError ElNNLSX_s
( ElConstMatrix_s A, ElConstMatrix_s b, ElMatrix_s x, 
  ElQPDirectCtrl_s ctrl );
EL_EXPORT ElError ElNNLSX_d
( ElConstMatrix_d A, ElConstMatrix_d b, ElMatrix_d x, 
  ElQPDirectCtrl_d ctrl );

EL_EXPORT ElError ElNNLSXDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s x,
  ElQPDirectCtrl_s ctrl );
EL_EXPORT ElError ElNNLSXDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d x,
  ElQPDirectCtrl_d ctrl );

EL_EXPORT ElError ElNNLSXDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, ElDistMatrix_s x,
  ElQPDirectCtrl_s ctrl );
EL_EXPORT ElError ElNNLSXDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, ElDistMatrix_d x,
  ElQPDirectCtrl_d ctrl );

EL_EXPORT ElError ElNNLSXDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s b, ElDistMultiVec_s x,
  ElQPDirectCtrl_s ctrl );
EL_EXPORT ElError ElNNLSXDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d b, ElDistMultiVec_d x,
  ElQPDirectCtrl_d ctrl );

/* ADMM
   ---- */
EL_EXPORT ElError ElNNLSADMM_s
( ElConstMatrix_s A, ElConstMatrix_s B, ElMatrix_s X, ElInt* numIts );
EL_EXPORT ElError ElNNLSADMM_d
( ElConstMatrix_d A, ElConstMatrix_d B, ElMatrix_d X, ElInt* numIts );

EL_EXPORT ElError ElNNLSADMMDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, ElDistMatrix_s X, 
  ElInt* numIts );
EL_EXPORT ElError ElNNLSADMMDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, ElDistMatrix_d X, 
  ElInt* numIts );

/* TODO: Expert versions */

/* Basis pursuit denoising
   ======================= */
EL_EXPORT ElError ElBPDN_s
( ElConstMatrix_s A, ElConstMatrix_s b, float lambda, 
  ElMatrix_s x );
EL_EXPORT ElError ElBPDN_d
( ElConstMatrix_d A, ElConstMatrix_d b, double lambda, 
  ElMatrix_d x );

EL_EXPORT ElError ElBPDNDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, float lambda, 
  ElDistMatrix_s x );
EL_EXPORT ElError ElBPDNDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, double lambda, 
  ElDistMatrix_d x );

EL_EXPORT ElError ElBPDNSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s b, float lambda,
  ElMatrix_s x );
EL_EXPORT ElError ElBPDNSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d b, double lambda,
  ElMatrix_d x );

EL_EXPORT ElError ElBPDNDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s b, float lambda,
  ElDistMultiVec_s x );
EL_EXPORT ElError ElBPDNDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d b, double lambda,
  ElDistMultiVec_d x );

/* Expert verions
   -------------- */
EL_EXPORT ElError ElBPDNX_s
( ElConstMatrix_s A, ElConstMatrix_s b, float lambda,
  ElMatrix_s x, ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElBPDNX_d
( ElConstMatrix_d A, ElConstMatrix_d b, double lambda,
  ElMatrix_d x, ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElBPDNXDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, float lambda,
  ElDistMatrix_s x, ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElBPDNXDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, double lambda,
  ElDistMatrix_d x, ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElBPDNXSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s b, float lambda,
  ElMatrix_s x, ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElBPDNXSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d b, double lambda,
  ElMatrix_d x, ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElBPDNXDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s b, float lambda,
  ElDistMultiVec_s x, ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElBPDNXDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d b, double lambda,
  ElDistMultiVec_d x, ElQPAffineCtrl_d ctrl );

/* ADMM
   ---- */
EL_EXPORT ElError ElBPDNADMM_s
( ElConstMatrix_s A, ElConstMatrix_s b, float lambda, 
  ElMatrix_s z, ElInt* numIts );
EL_EXPORT ElError ElBPDNADMM_d
( ElConstMatrix_d A, ElConstMatrix_d b, double lambda, 
  ElMatrix_d z, ElInt* numIts );
EL_EXPORT ElError ElBPDNADMM_c
( ElConstMatrix_c A, ElConstMatrix_c b, float lambda, 
  ElMatrix_c z, ElInt* numIts );
EL_EXPORT ElError ElBPDNADMM_z
( ElConstMatrix_z A, ElConstMatrix_z b, double lambda, 
  ElMatrix_z z, ElInt* numIts );

EL_EXPORT ElError ElBPDNADMMDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, float lambda, 
  ElDistMatrix_s z, ElInt* numIts );
EL_EXPORT ElError ElBPDNADMMDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, double lambda, 
  ElDistMatrix_d z, ElInt* numIts );
EL_EXPORT ElError ElBPDNADMMDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c b, float lambda, 
  ElDistMatrix_c z, ElInt* numIts );
EL_EXPORT ElError ElBPDNADMMDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z b, double lambda, 
  ElDistMatrix_z z, ElInt* numIts );

// Elastic net (EN): 
//   min || b - A x ||_2^2 + lambda_1 || x ||_1 + lambda_2 || x ||_2^2
// ===================================================================
EL_EXPORT ElError ElEN_s
( ElConstMatrix_s A, ElConstMatrix_s b, 
  float lambda1, float lambda2,
  ElMatrix_s x );
EL_EXPORT ElError ElEN_d
( ElConstMatrix_d A, ElConstMatrix_d b, 
  double lambda1, double lambda2,
  ElMatrix_d x );

EL_EXPORT ElError ElENDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, 
  float lambda1, float lambda2,
  ElDistMatrix_s x );
EL_EXPORT ElError ElENDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, 
  double lambda1, double lambda2,
  ElDistMatrix_d x );

EL_EXPORT ElError ElENSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s b, 
  float lambda1, float lambda2,
  ElMatrix_s x );
EL_EXPORT ElError ElENSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d b, 
  double lambda1, double lambda2,
  ElMatrix_d x );

EL_EXPORT ElError ElENDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s b, 
  float lambda1, float lambda2,
  ElDistMultiVec_s x );
EL_EXPORT ElError ElENDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d b, 
  double lambda1, double lambda2,
  ElDistMultiVec_d x );

/* Expert verions
   -------------- */
EL_EXPORT ElError ElENX_s
( ElConstMatrix_s A, ElConstMatrix_s b, 
  float lambda1, float lambda2,
  ElMatrix_s x, ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElENX_d
( ElConstMatrix_d A, ElConstMatrix_d b, 
  double lambda1, double lambda2,
  ElMatrix_d x, ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElENXDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s b, 
  float lambda1, float lambda2,
  ElDistMatrix_s x, ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElENXDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d b, 
  double lambda1, double lambda2,
  ElDistMatrix_d x, ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElENXSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s b, 
  float lambda1, float lambda2,
  ElMatrix_s x, ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElENXSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d b, 
  double lambda1, double lambda2,
  ElMatrix_d x, ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElENXDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s b, 
  float lambda1, float lambda2,
  ElDistMultiVec_s x, ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElENXDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d b, 
  double lambda1, double lambda2,
  ElDistMultiVec_d x, ElQPAffineCtrl_d ctrl );

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
( ElConstMatrix_s A, ElConstMatrix_s d, float lambda, 
  ElMatrix_s x );
EL_EXPORT ElError ElSVM_d
( ElConstMatrix_d A, ElConstMatrix_d d, double lambda, 
  ElMatrix_d x );

EL_EXPORT ElError ElSVMDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s d, float lambda, 
  ElDistMatrix_s x );
EL_EXPORT ElError ElSVMDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d d, double lambda, 
  ElDistMatrix_d x );

EL_EXPORT ElError ElSVMSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s d, float lambda,
  ElMatrix_s x );
EL_EXPORT ElError ElSVMSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d d, double lambda,
  ElMatrix_d x );

EL_EXPORT ElError ElSVMDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s d, float lambda,
  ElDistMultiVec_s x );
EL_EXPORT ElError ElSVMDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d d, double lambda,
  ElDistMultiVec_d x );

/* Expert verions
   -------------- */
EL_EXPORT ElError ElSVMX_s
( ElConstMatrix_s A, ElConstMatrix_s d, float lambda,
  ElMatrix_s x, ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElSVMX_d
( ElConstMatrix_d A, ElConstMatrix_d d, double lambda,
  ElMatrix_d x, ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElSVMXDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s d, float lambda,
  ElDistMatrix_s x, ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElSVMXDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d d, double lambda,
  ElDistMatrix_d x, ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElSVMXSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s d, float lambda,
  ElMatrix_s x, ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElSVMXSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d d, double lambda,
  ElMatrix_d x, ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElSVMXDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s d, float lambda,
  ElDistMultiVec_s x, ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElSVMXDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d d, double lambda,
  ElDistMultiVec_d x, ElQPAffineCtrl_d ctrl );

/* ADMM
   ---- */
EL_EXPORT ElError ElSVMADMM_s
( ElConstMatrix_s G, ElConstMatrix_s q, 
  float gamma,       ElMatrix_s z, 
  ElInt* numIts );
EL_EXPORT ElError ElSVMADMM_d
( ElConstMatrix_d G, ElConstMatrix_d q, 
  double gamma,      ElMatrix_d z, 
  ElInt* numIts );

EL_EXPORT ElError ElSVMADMMDist_s
( ElConstDistMatrix_s G, ElConstDistMatrix_s q, 
  float gamma,           ElDistMatrix_s z, 
  ElInt* numIts );
EL_EXPORT ElError ElSVMADMMDist_d
( ElConstDistMatrix_d G, ElConstDistMatrix_d q, 
  double gamma,          ElDistMatrix_d z,
  ElInt* numIts );

/* TODO: Expert versions */

/* Total variation denoising
   ========================= */
EL_EXPORT ElError ElTV_s
( ElConstMatrix_s b, float lambda, ElMatrix_s x );
EL_EXPORT ElError ElTV_d
( ElConstMatrix_d b, double lambda, 
  ElMatrix_d x );

EL_EXPORT ElError ElTVDist_s
( ElConstDistMatrix_s b, float lambda, ElDistMatrix_s x );
EL_EXPORT ElError ElTVDist_d
( ElConstDistMatrix_d b, double lambda, ElDistMatrix_d x );

EL_EXPORT ElError ElTVSparse_s
( ElConstMatrix_s b, float lambda, ElMatrix_s x );
EL_EXPORT ElError ElTVSparse_d
( ElConstMatrix_d b, double lambda, ElMatrix_d x );

EL_EXPORT ElError ElTVDistSparse_s
( ElConstDistMultiVec_s b, float lambda, ElDistMultiVec_s x );
EL_EXPORT ElError ElTVDistSparse_d
( ElConstDistMultiVec_d b, double lambda, ElDistMultiVec_d x );

/* Expert verions
   -------------- */
EL_EXPORT ElError ElTVX_s
( ElConstMatrix_s b, float lambda, ElMatrix_s x, 
  ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElTVX_d
( ElConstMatrix_d b, double lambda, ElMatrix_d x, 
  ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElTVXDist_s
( ElConstDistMatrix_s b, float lambda, ElDistMatrix_s x, 
  ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElTVXDist_d
( ElConstDistMatrix_d b, double lambda, ElDistMatrix_d x, 
  ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElTVXSparse_s
( ElConstMatrix_s b, float lambda, ElMatrix_s x, 
  ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElTVXSparse_d
( ElConstMatrix_d b, double lambda, ElMatrix_d x, 
  ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElTVXDistSparse_s
( ElConstDistMultiVec_s b, float lambda, ElDistMultiVec_s x, 
  ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElTVXDistSparse_d
( ElConstDistMultiVec_d b, double lambda, ElDistMultiVec_d x, 
  ElQPAffineCtrl_d ctrl );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_OPTIMIZATION_MODELS_C_H */
