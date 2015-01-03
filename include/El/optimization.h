/*
   Copyright (c) 2009-2015, Jack Poulson
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

/* Linear Program control structures
   ================================= */
typedef enum {
  EL_LP_ADMM,
  EL_LP_IPF,
  EL_LP_IPF_SELFDUAL,
  EL_LP_MEHROTRA,
  EL_LP_MEHROTRA_SELFDUAL
} ElLPApproach;

typedef struct {
  float gamma;
  float beta;
  float psi;
  float stepRatio;
  bool print;
} ElLPIPFLineSearchCtrl_s;

typedef struct {
  double gamma;
  double beta;
  double psi;
  double stepRatio;
  bool print;
} ElLPIPFLineSearchCtrl_d;

ElError ElLPIPFLineSearchCtrlDefault_s( ElLPIPFLineSearchCtrl_s* ctrl );
ElError ElLPIPFLineSearchCtrlDefault_d( ElLPIPFLineSearchCtrl_d* ctrl );

/* Direct conic form
   ----------------- */
typedef enum {
  EL_LP_PRIMAL_FULL_KKT,
  EL_LP_PRIMAL_AUGMENTED_KKT,
  EL_LP_PRIMAL_NORMAL_KKT
} ElLPDirectKKTSystem;

typedef struct {
  float rho;  
  float alpha;
  ElInt maxIter;
  float absTol;
  float relTol;
  bool inv;
  bool print;
} ElLPDirectADMMCtrl_s;

typedef struct {
  double rho;  
  double alpha;
  ElInt maxIter;
  double absTol;
  double relTol;
  bool inv;
  bool print;
} ElLPDirectADMMCtrl_d;

ElError ElLPDirectADMMCtrlDefault_s( ElLPDirectADMMCtrl_s* ctrl );
ElError ElLPDirectADMMCtrlDefault_d( ElLPDirectADMMCtrl_d* ctrl );

typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  float tol;
  ElInt maxIts;
  float centering;
  ElLPDirectKKTSystem system;

  ElLPIPFLineSearchCtrl_s lineSearchCtrl;
  bool print;
} ElLPDirectIPFCtrl_s;

typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  double tol;
  ElInt maxIts;
  double centering;
  ElLPDirectKKTSystem system;

  ElLPIPFLineSearchCtrl_d lineSearchCtrl;
  bool print;
} ElLPDirectIPFCtrl_d;

ElError ElLPDirectIPFCtrlDefault_s( ElLPDirectIPFCtrl_s* ctrl, bool isSparse );
ElError ElLPDirectIPFCtrlDefault_d( ElLPDirectIPFCtrl_d* ctrl, bool isSparse );

typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  float tol;
  ElInt maxIts;
  float maxStepRatio;
  ElLPDirectKKTSystem system;
  bool print;
} ElLPDirectMehrotraCtrl_s;

typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  double tol;
  ElInt maxIts;
  double maxStepRatio;
  ElLPDirectKKTSystem system;
  bool print;
} ElLPDirectMehrotraCtrl_d;

ElError ElLPDirectMehrotraCtrlDefault_s
( ElLPDirectMehrotraCtrl_s* ctrl, bool isSparse );
ElError ElLPDirectMehrotraCtrlDefault_d
( ElLPDirectMehrotraCtrl_d* ctrl, bool isSparse );

typedef struct {
  ElLPApproach approach;  
  ElLPDirectADMMCtrl_s admmCtrl;
  ElLPDirectIPFCtrl_s ipfCtrl;
  ElLPDirectMehrotraCtrl_s mehrotraCtrl;
} ElLPDirectCtrl_s;
typedef struct {
  ElLPApproach approach;  
  ElLPDirectADMMCtrl_d admmCtrl;
  ElLPDirectIPFCtrl_d ipfCtrl;
  ElLPDirectMehrotraCtrl_d mehrotraCtrl;
} ElLPDirectCtrl_d;

ElError ElLPDirectCtrlDefault_s( ElLPDirectCtrl_s* ctrl, bool isSparse );
ElError ElLPDirectCtrlDefault_d( ElLPDirectCtrl_d* ctrl, bool isSparse );

/* Affine conic form
   ----------------- */
typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  float tol;
  ElInt maxIts;
  float centering;

  ElLPIPFLineSearchCtrl_s lineSearchCtrl;
  bool print;
} ElLPAffineIPFCtrl_s;

typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  double tol;
  ElInt maxIts;
  double centering;

  ElLPIPFLineSearchCtrl_d lineSearchCtrl;
  bool print;
} ElLPAffineIPFCtrl_d;

ElError ElLPAffineIPFCtrlDefault_s( ElLPAffineIPFCtrl_s* ctrl );
ElError ElLPAffineIPFCtrlDefault_d( ElLPAffineIPFCtrl_d* ctrl );

typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  float tol;
  ElInt maxIts;
  float maxStepRatio;
  bool print;
} ElLPAffineMehrotraCtrl_s;

typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  double tol;
  ElInt maxIts;
  double maxStepRatio;
  bool print;
} ElLPAffineMehrotraCtrl_d;

ElError ElLPAffineMehrotraCtrlDefault_s( ElLPAffineMehrotraCtrl_s* ctrl );
ElError ElLPAffineMehrotraCtrlDefault_d( ElLPAffineMehrotraCtrl_d* ctrl );

typedef struct {
  ElLPApproach approach;  
  ElLPAffineIPFCtrl_s ipfCtrl;
  ElLPAffineMehrotraCtrl_s mehrotraCtrl;
} ElLPAffineCtrl_s;
typedef struct {
  ElLPApproach approach;  
  ElLPAffineIPFCtrl_d ipfCtrl;
  ElLPAffineMehrotraCtrl_d mehrotraCtrl;
} ElLPAffineCtrl_d;

ElError ElLPAffineCtrlDefault_s( ElLPAffineCtrl_s* ctrl );
ElError ElLPAffineCtrlDefault_d( ElLPAffineCtrl_d* ctrl );

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
/* Direct conic form
   ----------------- */ 
EL_EXPORT ElError ElLPDirect_s
( ElConstMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c, 
  ElMatrix_s x,      ElMatrix_s y, 
  ElMatrix_s z );
EL_EXPORT ElError ElLPDirect_d
( ElConstMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c, 
  ElMatrix_d x,      ElMatrix_d y, 
  ElMatrix_d z );

EL_EXPORT ElError ElLPDirectDist_s
( ElConstDistMatrix_s A, 
  ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElDistMatrix_s x,      ElDistMatrix_s y, 
  ElDistMatrix_s z );
EL_EXPORT ElError ElLPDirectDist_d
( ElConstDistMatrix_d A, 
  ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElDistMatrix_d x,      ElDistMatrix_d y, 
  ElDistMatrix_d z );

EL_EXPORT ElError ElLPDirectSparse_s
( ElConstSparseMatrix_s A, 
  ElConstMatrix_s b,       ElConstMatrix_s c, 
  ElMatrix_s x,            ElMatrix_s y, 
  ElMatrix_s z );
EL_EXPORT ElError ElLPDirectSparse_d
( ElConstSparseMatrix_d A, 
  ElConstMatrix_d b,       ElConstMatrix_d c, 
  ElMatrix_d x,            ElMatrix_d y, 
  ElMatrix_d z );

EL_EXPORT ElError ElLPDirectDistSparse_s
( ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s b,     ElConstDistMultiVec_s c, 
  ElDistMultiVec_s x,          ElDistMultiVec_s y, 
  ElDistMultiVec_s z );
EL_EXPORT ElError ElLPDirectDistSparse_d
( ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d b,     ElConstDistMultiVec_d c, 
  ElDistMultiVec_d x,          ElDistMultiVec_d y, 
  ElDistMultiVec_d z );

/* Expert versions
   ^^^^^^^^^^^^^^^ */
EL_EXPORT ElError ElLPDirectX_s
( ElConstMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c, 
  ElMatrix_s x,      ElMatrix_s y, 
  ElMatrix_s z, 
  ElLPDirectCtrl_s ctrl );
EL_EXPORT ElError ElLPDirectX_d
( ElConstMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c, 
  ElMatrix_d x,      ElMatrix_d y, 
  ElMatrix_d z, 
  ElLPDirectCtrl_d ctrl );

EL_EXPORT ElError ElLPDirectXDist_s
( ElConstDistMatrix_s A, 
  ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElDistMatrix_s x,      ElDistMatrix_s y, 
  ElDistMatrix_s z, 
  ElLPDirectCtrl_s ctrl );
EL_EXPORT ElError ElLPDirectXDist_d
( ElConstDistMatrix_d A, 
  ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElDistMatrix_d x,      ElDistMatrix_d y, 
  ElDistMatrix_d z, 
  ElLPDirectCtrl_d ctrl );

EL_EXPORT ElError ElLPDirectXSparse_s
( ElConstSparseMatrix_s A, 
  ElConstMatrix_s b,       ElConstMatrix_s c, 
  ElMatrix_s x,            ElMatrix_s y, 
  ElMatrix_s z, 
  ElLPDirectCtrl_s ctrl );
EL_EXPORT ElError ElLPDirectXSparse_d
( ElConstSparseMatrix_d A, 
  ElConstMatrix_d b,       ElConstMatrix_d c, 
  ElMatrix_d x,            ElMatrix_d y, 
  ElMatrix_d z, 
  ElLPDirectCtrl_d ctrl );

EL_EXPORT ElError ElLPDirectXDistSparse_s
( ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s b,     ElConstDistMultiVec_s c, 
  ElDistMultiVec_s x,          ElDistMultiVec_s y, 
  ElDistMultiVec_s z, 
  ElLPDirectCtrl_s ctrl );
EL_EXPORT ElError ElLPDirectXDistSparse_d
( ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d b,     ElConstDistMultiVec_d c, 
  ElDistMultiVec_d x,          ElDistMultiVec_d y, 
  ElDistMultiVec_d z,
  ElLPDirectCtrl_d ctrl );

/* Affine conic form
   ----------------- */
EL_EXPORT ElError ElLPAffine_s
( ElConstMatrix_s A, ElConstMatrix_s G,
  ElConstMatrix_s b, ElConstMatrix_s c, 
  ElConstMatrix_s h,
  ElMatrix_s x,      ElMatrix_s y, 
  ElMatrix_s z,      ElMatrix_s s );
EL_EXPORT ElError ElLPAffine_d
( ElConstMatrix_d A, ElConstMatrix_d G,
  ElConstMatrix_d b, ElConstMatrix_d c, 
  ElConstMatrix_d h,
  ElMatrix_d x,      ElMatrix_d y, 
  ElMatrix_d z,      ElMatrix_d s );

EL_EXPORT ElError ElLPAffineDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s G,
  ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElConstDistMatrix_s h,
  ElDistMatrix_s x,      ElDistMatrix_s y, 
  ElDistMatrix_s z,      ElDistMatrix_s s );
EL_EXPORT ElError ElLPAffineDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d G,
  ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElConstDistMatrix_d h,
  ElDistMatrix_d x,      ElDistMatrix_d y, 
  ElDistMatrix_d z,      ElDistMatrix_d s );

EL_EXPORT ElError ElLPAffineSparse_s
( ElConstSparseMatrix_s A, ElConstSparseMatrix_s G,
  ElConstMatrix_s b,       ElConstMatrix_s c, 
  ElConstMatrix_s h,
  ElMatrix_s x,            ElMatrix_s y, 
  ElMatrix_s z,            ElMatrix_s s );
EL_EXPORT ElError ElLPAffineSparse_d
( ElConstSparseMatrix_d A, ElConstSparseMatrix_d G,
  ElConstMatrix_d b,       ElConstMatrix_d c, 
  ElConstMatrix_d h,
  ElMatrix_d x,            ElMatrix_d y, 
  ElMatrix_d z,            ElMatrix_d s );

EL_EXPORT ElError ElLPAffineDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistSparseMatrix_s G,
  ElConstDistMultiVec_s b,     ElConstDistMultiVec_s c, 
  ElConstDistMultiVec_s h,
  ElDistMultiVec_s x,          ElDistMultiVec_s y, 
  ElDistMultiVec_s z,          ElDistMultiVec_s s );
EL_EXPORT ElError ElLPAffineDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistSparseMatrix_d G,
  ElConstDistMultiVec_d b,     ElConstDistMultiVec_d c, 
  ElConstDistMultiVec_d h,
  ElDistMultiVec_d x,          ElDistMultiVec_d y, 
  ElDistMultiVec_d z,          ElDistMultiVec_d s );

/* Expert versions
   ^^^^^^^^^^^^^^^ */
EL_EXPORT ElError ElLPAffineX_s
( ElConstMatrix_s A, ElConstMatrix_s G,
  ElConstMatrix_s b, ElConstMatrix_s c, 
  ElConstMatrix_s h,
  ElMatrix_s x,      ElMatrix_s y, 
  ElMatrix_s z,      ElMatrix_s s,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElLPAffineX_d
( ElConstMatrix_d A, ElConstMatrix_d G,
  ElConstMatrix_d b, ElConstMatrix_d c, 
  ElConstMatrix_d h,
  ElMatrix_d x,      ElMatrix_d y, 
  ElMatrix_d z,      ElMatrix_d s,
  ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElLPAffineXDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s G,
  ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElConstDistMatrix_s h,
  ElDistMatrix_s x,      ElDistMatrix_s y, 
  ElDistMatrix_s z,      ElDistMatrix_s s,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElLPAffineXDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d G,
  ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElConstDistMatrix_d h,
  ElDistMatrix_d x,      ElDistMatrix_d y, 
  ElDistMatrix_d z,      ElDistMatrix_d s,
  ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElLPAffineXSparse_s
( ElConstSparseMatrix_s A, ElConstSparseMatrix_s G,
  ElConstMatrix_s b,       ElConstMatrix_s c, 
  ElConstMatrix_s h,
  ElMatrix_s x,            ElMatrix_s y, 
  ElMatrix_s z,            ElMatrix_s s,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElLPAffineXSparse_d
( ElConstSparseMatrix_d A, ElConstSparseMatrix_d G,
  ElConstMatrix_d b,       ElConstMatrix_d c, 
  ElConstMatrix_d h,
  ElMatrix_d x,            ElMatrix_d y, 
  ElMatrix_d z,            ElMatrix_d s,
  ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElLPAffineXDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistSparseMatrix_s G,
  ElConstDistMultiVec_s b,     ElConstDistMultiVec_s c, 
  ElConstDistMultiVec_s h,
  ElDistMultiVec_s x,          ElDistMultiVec_s y, 
  ElDistMultiVec_s z,          ElDistMultiVec_s s,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElLPAffineXDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistSparseMatrix_d G,
  ElConstDistMultiVec_d b,     ElConstDistMultiVec_d c, 
  ElConstDistMultiVec_d h,
  ElDistMultiVec_d x,          ElDistMultiVec_d y, 
  ElDistMultiVec_d z,          ElDistMultiVec_d s,
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
typedef enum {
  EL_QP_ADMM,
  EL_QP_IPF,
  EL_QP_IPF_SELFDUAL,
  EL_QP_MEHROTRA,
  EL_QP_MEHROTRA_SELFDUAL
} ElQPApproach;

typedef struct {
  float gamma;
  float beta;
  float psi;
  float stepRatio;
  bool print;
} ElQPIPFLineSearchCtrl_s;

typedef struct {
  double gamma;
  double beta;
  double psi;
  double stepRatio;
  bool print;
} ElQPIPFLineSearchCtrl_d;

ElError ElQPIPFLineSearchCtrlDefault_s( ElQPIPFLineSearchCtrl_s* ctrl );
ElError ElQPIPFLineSearchCtrlDefault_d( ElQPIPFLineSearchCtrl_d* ctrl );

/* Direct conic form
   ----------------- */ 
EL_EXPORT ElError ElQPDirect_s
( ElConstMatrix_s Q, ElConstMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c, 
  ElMatrix_s x,      ElMatrix_s y, 
  ElMatrix_s z );
EL_EXPORT ElError ElQPDirect_d
( ElConstMatrix_d Q, ElConstMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c, 
  ElMatrix_d x,      ElMatrix_d y, 
  ElMatrix_d z );

EL_EXPORT ElError ElQPDirectDist_s
( ElConstDistMatrix_s Q, ElConstDistMatrix_s A, 
  ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElDistMatrix_s x,      ElDistMatrix_s y, 
  ElDistMatrix_s z );
EL_EXPORT ElError ElQPDirectDist_d
( ElConstDistMatrix_d Q, ElConstDistMatrix_d A, 
  ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElDistMatrix_d x,      ElDistMatrix_d y, 
  ElDistMatrix_d z );

EL_EXPORT ElError ElQPDirectSparse_s
( ElConstSparseMatrix_s Q, ElConstSparseMatrix_s A, 
  ElConstMatrix_s b,       ElConstMatrix_s c, 
  ElMatrix_s x,            ElMatrix_s y, 
  ElMatrix_s z );
EL_EXPORT ElError ElQPDirectSparse_d
( ElConstSparseMatrix_d Q, ElConstSparseMatrix_d A, 
  ElConstMatrix_d b,       ElConstMatrix_d c, 
  ElMatrix_d x,            ElMatrix_d y, 
  ElMatrix_d z );

EL_EXPORT ElError ElQPDirectDistSparse_s
( ElConstDistSparseMatrix_s Q, ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s b,     ElConstDistMultiVec_s c, 
  ElDistMultiVec_s x,          ElDistMultiVec_s y, 
  ElDistMultiVec_s z );
EL_EXPORT ElError ElQPDirectDistSparse_d
( ElConstDistSparseMatrix_d Q, ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d b,     ElConstDistMultiVec_d c, 
  ElDistMultiVec_d x,          ElDistMultiVec_d y, 
  ElDistMultiVec_d z );

/* Expert versions
   ^^^^^^^^^^^^^^^ */
typedef enum {
  EL_QP_PRIMAL_FULL_KKT,
  EL_QP_PRIMAL_AUGMENTED_KKT
} ElQPDirectKKTSystem;

typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  float tol;
  ElInt maxIts;
  float centering;
  ElQPDirectKKTSystem system;

  ElQPIPFLineSearchCtrl_s lineSearchCtrl;
  bool print;
} ElQPDirectIPFCtrl_s;

typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  double tol;
  ElInt maxIts;
  double centering;
  ElQPDirectKKTSystem system;

  ElQPIPFLineSearchCtrl_d lineSearchCtrl;
  bool print;
} ElQPDirectIPFCtrl_d;

ElError ElQPDirectIPFCtrlDefault_s( ElQPDirectIPFCtrl_s* ctrl );
ElError ElQPDirectIPFCtrlDefault_d( ElQPDirectIPFCtrl_d* ctrl );

typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  float tol;
  ElInt maxIts;
  float maxStepRatio;
  ElQPDirectKKTSystem system;
  bool print;
} ElQPDirectMehrotraCtrl_s;

typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  double tol;
  ElInt maxIts;
  double maxStepRatio;
  ElQPDirectKKTSystem system;
  bool print;
} ElQPDirectMehrotraCtrl_d;

ElError ElQPDirectMehrotraCtrlDefault_s( ElQPDirectMehrotraCtrl_s* ctrl );
ElError ElQPDirectMehrotraCtrlDefault_d( ElQPDirectMehrotraCtrl_d* ctrl );

typedef struct {
  ElQPApproach approach;  
  ElQPDirectIPFCtrl_s ipfCtrl;
  ElQPDirectMehrotraCtrl_s mehrotraCtrl;
} ElQPDirectCtrl_s;
typedef struct {
  ElQPApproach approach;  
  ElQPDirectIPFCtrl_d ipfCtrl;
  ElQPDirectMehrotraCtrl_d mehrotraCtrl;
} ElQPDirectCtrl_d;

ElError ElQPDirectCtrlDefault_s( ElQPDirectCtrl_s* ctrl );
ElError ElQPDirectCtrlDefault_d( ElQPDirectCtrl_d* ctrl );

EL_EXPORT ElError ElQPDirectX_s
( ElConstMatrix_s Q, ElConstMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c, 
  ElMatrix_s x,      ElMatrix_s y, 
  ElMatrix_s z, 
  ElQPDirectCtrl_s ctrl );
EL_EXPORT ElError ElQPDirectX_d
( ElConstMatrix_d Q, ElConstMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c, 
  ElMatrix_d x,      ElMatrix_d y, 
  ElMatrix_d z, 
  ElQPDirectCtrl_d ctrl );

EL_EXPORT ElError ElQPDirectXDist_s
( ElConstDistMatrix_s Q, ElConstDistMatrix_s A, 
  ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElDistMatrix_s x,      ElDistMatrix_s y, 
  ElDistMatrix_s z, 
  ElQPDirectCtrl_s ctrl );
EL_EXPORT ElError ElQPDirectXDist_d
( ElConstDistMatrix_d Q, ElConstDistMatrix_d A, 
  ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElDistMatrix_d x,      ElDistMatrix_d y, 
  ElDistMatrix_d z, 
  ElQPDirectCtrl_d ctrl );

EL_EXPORT ElError ElQPDirectXSparse_s
( ElConstSparseMatrix_s Q, ElConstSparseMatrix_s A, 
  ElConstMatrix_s b,       ElConstMatrix_s c, 
  ElMatrix_s x,            ElMatrix_s y, 
  ElMatrix_s z, 
  ElQPDirectCtrl_s ctrl );
EL_EXPORT ElError ElQPDirectXSparse_d
( ElConstSparseMatrix_d Q, ElConstSparseMatrix_d A, 
  ElConstMatrix_d b,       ElConstMatrix_d c, 
  ElMatrix_d x,            ElMatrix_d y, 
  ElMatrix_d z, 
  ElQPDirectCtrl_d ctrl );

EL_EXPORT ElError ElQPDirectXDistSparse_s
( ElConstDistSparseMatrix_s Q, ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s b,     ElConstDistMultiVec_s c, 
  ElDistMultiVec_s x,          ElDistMultiVec_s y, 
  ElDistMultiVec_s z, 
  ElQPDirectCtrl_s ctrl );
EL_EXPORT ElError ElQPDirectXDistSparse_d
( ElConstDistSparseMatrix_d Q, ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d b,     ElConstDistMultiVec_d c, 
  ElDistMultiVec_d x,          ElDistMultiVec_d y, 
  ElDistMultiVec_d z,
  ElQPDirectCtrl_d ctrl );

/* Affine conic form
   ----------------- */
EL_EXPORT ElError ElQPAffine_s
( ElConstMatrix_s Q,
  ElConstMatrix_s A, ElConstMatrix_s G,
  ElConstMatrix_s b, ElConstMatrix_s c, 
  ElConstMatrix_s h,
  ElMatrix_s x,      ElMatrix_s y, 
  ElMatrix_s z,      ElMatrix_s s );
EL_EXPORT ElError ElQPAffine_d
( ElConstMatrix_d Q,
  ElConstMatrix_d A, ElConstMatrix_d G,
  ElConstMatrix_d b, ElConstMatrix_d c, 
  ElConstMatrix_d h,
  ElMatrix_d x,      ElMatrix_d y, 
  ElMatrix_d z,      ElMatrix_d s );

EL_EXPORT ElError ElQPAffineDist_s
( ElConstDistMatrix_s Q,
  ElConstDistMatrix_s A, ElConstDistMatrix_s G,
  ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElConstDistMatrix_s h,
  ElDistMatrix_s x,      ElDistMatrix_s y, 
  ElDistMatrix_s z,      ElDistMatrix_s s );
EL_EXPORT ElError ElQPAffineDist_d
( ElConstDistMatrix_d Q,
  ElConstDistMatrix_d A, ElConstDistMatrix_d G,
  ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElConstDistMatrix_d h,
  ElDistMatrix_d x,      ElDistMatrix_d y, 
  ElDistMatrix_d z,      ElDistMatrix_d s );

EL_EXPORT ElError ElQPAffineSparse_s
( ElConstSparseMatrix_s Q,
  ElConstSparseMatrix_s A, ElConstSparseMatrix_s G,
  ElConstMatrix_s b,       ElConstMatrix_s c, 
  ElConstMatrix_s h,
  ElMatrix_s x,            ElMatrix_s y, 
  ElMatrix_s z,            ElMatrix_s s );
EL_EXPORT ElError ElQPAffineSparse_d
( ElConstSparseMatrix_d Q,
  ElConstSparseMatrix_d A, ElConstSparseMatrix_d G,
  ElConstMatrix_d b,       ElConstMatrix_d c, 
  ElConstMatrix_d h,
  ElMatrix_d x,            ElMatrix_d y, 
  ElMatrix_d z,            ElMatrix_d s );

EL_EXPORT ElError ElQPAffineDistSparse_s
( ElConstDistSparseMatrix_s Q,
  ElConstDistSparseMatrix_s A, ElConstDistSparseMatrix_s G,
  ElConstDistMultiVec_s b,     ElConstDistMultiVec_s c, 
  ElConstDistMultiVec_s h,
  ElDistMultiVec_s x,          ElDistMultiVec_s y, 
  ElDistMultiVec_s z,          ElDistMultiVec_s s );
EL_EXPORT ElError ElQPAffineDistSparse_d
( ElConstDistSparseMatrix_d Q,
  ElConstDistSparseMatrix_d A, ElConstDistSparseMatrix_d G,
  ElConstDistMultiVec_d b,     ElConstDistMultiVec_d c, 
  ElConstDistMultiVec_d h,
  ElDistMultiVec_d x,          ElDistMultiVec_d y, 
  ElDistMultiVec_d z,          ElDistMultiVec_d s );

/* Expert versions
   ^^^^^^^^^^^^^^^ */
typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  float tol;
  ElInt maxIts;
  float centering;

  ElQPIPFLineSearchCtrl_s lineSearchCtrl;
  bool print;
} ElQPAffineIPFCtrl_s;

typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  double tol;
  ElInt maxIts;
  double centering;

  ElQPIPFLineSearchCtrl_d lineSearchCtrl;
  bool print;
} ElQPAffineIPFCtrl_d;

ElError ElQPAffineIPFCtrlDefault_s( ElQPAffineIPFCtrl_s* ctrl );
ElError ElQPAffineIPFCtrlDefault_d( ElQPAffineIPFCtrl_d* ctrl );

typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  float tol;
  ElInt maxIts;
  float maxStepRatio;
  bool print;
} ElQPAffineMehrotraCtrl_s;

typedef struct {
  bool primalInitialized;
  bool dualInitialized;
  double tol;
  ElInt maxIts;
  double maxStepRatio;
  bool print;
} ElQPAffineMehrotraCtrl_d;

ElError ElQPAffineMehrotraCtrlDefault_s( ElQPAffineMehrotraCtrl_s* ctrl );
ElError ElQPAffineMehrotraCtrlDefault_d( ElQPAffineMehrotraCtrl_d* ctrl );

typedef struct {
  ElQPApproach approach;  
  ElQPAffineIPFCtrl_s ipfCtrl;
  ElQPAffineMehrotraCtrl_s mehrotraCtrl;
} ElQPAffineCtrl_s;
typedef struct {
  ElQPApproach approach;  
  ElQPAffineIPFCtrl_d ipfCtrl;
  ElQPAffineMehrotraCtrl_d mehrotraCtrl;
} ElQPAffineCtrl_d;

ElError ElQPAffineCtrlDefault_s( ElQPAffineCtrl_s* ctrl );
ElError ElQPAffineCtrlDefault_d( ElQPAffineCtrl_d* ctrl );

EL_EXPORT ElError ElQPAffineX_s
( ElConstMatrix_s Q,
  ElConstMatrix_s A, ElConstMatrix_s G,
  ElConstMatrix_s b, ElConstMatrix_s c, 
  ElConstMatrix_s h,
  ElMatrix_s x,      ElMatrix_s y, 
  ElMatrix_s z,      ElMatrix_s s,
  ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElQPAffineX_d
( ElConstMatrix_d Q,
  ElConstMatrix_d A, ElConstMatrix_d G,
  ElConstMatrix_d b, ElConstMatrix_d c, 
  ElConstMatrix_d h,
  ElMatrix_d x,      ElMatrix_d y, 
  ElMatrix_d z,      ElMatrix_d s,
  ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElQPAffineXDist_s
( ElConstDistMatrix_s Q,
  ElConstDistMatrix_s A, ElConstDistMatrix_s G,
  ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElConstDistMatrix_s h,
  ElDistMatrix_s x,      ElDistMatrix_s y, 
  ElDistMatrix_s z,      ElDistMatrix_s s,
  ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElQPAffineXDist_d
( ElConstDistMatrix_d Q,
  ElConstDistMatrix_d A, ElConstDistMatrix_d G,
  ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElConstDistMatrix_d h,
  ElDistMatrix_d x,      ElDistMatrix_d y, 
  ElDistMatrix_d z,      ElDistMatrix_d s,
  ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElQPAffineXSparse_s
( ElConstSparseMatrix_s Q,
  ElConstSparseMatrix_s A, ElConstSparseMatrix_s G,
  ElConstMatrix_s b,       ElConstMatrix_s c, 
  ElConstMatrix_s h,
  ElMatrix_s x,            ElMatrix_s y, 
  ElMatrix_s z,            ElMatrix_s s,
  ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElQPAffineXSparse_d
( ElConstSparseMatrix_d Q,
  ElConstSparseMatrix_d A, ElConstSparseMatrix_d G,
  ElConstMatrix_d b,       ElConstMatrix_d c, 
  ElConstMatrix_d h,
  ElMatrix_d x,            ElMatrix_d y, 
  ElMatrix_d z,            ElMatrix_d s,
  ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElQPAffineXDistSparse_s
( ElConstDistSparseMatrix_s Q,
  ElConstDistSparseMatrix_s A, ElConstDistSparseMatrix_s G,
  ElConstDistMultiVec_s b,     ElConstDistMultiVec_s c, 
  ElConstDistMultiVec_s h,
  ElDistMultiVec_s x,          ElDistMultiVec_s y, 
  ElDistMultiVec_s z,          ElDistMultiVec_s s,
  ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElQPAffineXDistSparse_d
( ElConstDistSparseMatrix_d Q,
  ElConstDistSparseMatrix_d A, ElConstDistSparseMatrix_d G,
  ElConstDistMultiVec_d b,     ElConstDistMultiVec_d c, 
  ElConstDistMultiVec_d h,
  ElDistMultiVec_d x,          ElDistMultiVec_d y, 
  ElDistMultiVec_d z,          ElDistMultiVec_d s,
  ElQPAffineCtrl_d ctrl );

/* Box form
   -------- */

/* Alternating Direction Method of Multipliers
   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
EL_EXPORT ElError ElQPBoxADMM_s
( ElConstMatrix_s Q, ElConstMatrix_s C, float lb, float ub, 
  ElMatrix_s Z, ElInt* numIts );
EL_EXPORT ElError ElQPBoxADMM_d
( ElConstMatrix_d Q, ElConstMatrix_d C, double lb, double ub, 
  ElMatrix_d Z, ElInt* numIts );

EL_EXPORT ElError ElQPBoxADMMDist_s
( ElConstDistMatrix_s Q, ElConstDistMatrix_s C, float lb, float ub, 
  ElDistMatrix_s Z, ElInt* numIts );
EL_EXPORT ElError ElQPBoxADMMDist_d
( ElConstDistMatrix_d Q, ElConstDistMatrix_d C, double lb, double ub, 
  ElDistMatrix_d Z, ElInt* numIts );

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
