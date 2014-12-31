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
typedef enum {
  EL_LP_ADMM,
  EL_LP_IPF,
  EL_LP_IPF_SELFDUAL,
  EL_LP_MEHROTRA,
  EL_LP_MEHROTRA_SELFDUAL
} ElLPApproach;

/* Primal conic form
   ----------------- */ 
EL_EXPORT ElError ElLPPrimal_s
( ElConstMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c, 
  ElMatrix_s x,      ElMatrix_s y, 
  ElMatrix_s z );
EL_EXPORT ElError ElLPPrimal_d
( ElConstMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c, 
  ElMatrix_d x,      ElMatrix_d y, 
  ElMatrix_d z );

EL_EXPORT ElError ElLPPrimalDist_s
( ElConstDistMatrix_s A, 
  ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElDistMatrix_s x,      ElDistMatrix_s y, 
  ElDistMatrix_s z );
EL_EXPORT ElError ElLPPrimalDist_d
( ElConstDistMatrix_d A, 
  ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElDistMatrix_d x,      ElDistMatrix_d y, 
  ElDistMatrix_d z );

EL_EXPORT ElError ElLPPrimalSparse_s
( ElConstSparseMatrix_s A, 
  ElConstMatrix_s b,       ElConstMatrix_s c, 
  ElMatrix_s x,            ElMatrix_s y, 
  ElMatrix_s z );
EL_EXPORT ElError ElLPPrimalSparse_d
( ElConstSparseMatrix_d A, 
  ElConstMatrix_d b,       ElConstMatrix_d c, 
  ElMatrix_d x,            ElMatrix_d y, 
  ElMatrix_d z );

EL_EXPORT ElError ElLPPrimalDistSparse_s
( ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s b,     ElConstDistMultiVec_s c, 
  ElDistMultiVec_s x,          ElDistMultiVec_s y, 
  ElDistMultiVec_s z );
EL_EXPORT ElError ElLPPrimalDistSparse_d
( ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d b,     ElConstDistMultiVec_d c, 
  ElDistMultiVec_d x,          ElDistMultiVec_d y, 
  ElDistMultiVec_d z );

/* Expert versions
   ^^^^^^^^^^^^^^^ */
typedef enum {
  EL_LP_PRIMAL_FULL_KKT,
  EL_LP_PRIMAL_AUGMENTED_KKT,
  EL_LP_PRIMAL_NORMAL_KKT
} ElLPPrimalKKTSystem;

typedef struct {
  float rho;  
  float alpha;
  ElInt maxIter;
  float absTol;
  float relTol;
  bool inv;
  bool print;
} ElLPPrimalADMMCtrl_s;

typedef struct {
  double rho;  
  double alpha;
  ElInt maxIter;
  double absTol;
  double relTol;
  bool inv;
  bool print;
} ElLPPrimalADMMCtrl_d;

ElError ElLPPrimalADMMCtrlDefault_s( ElLPPrimalADMMCtrl_s* ctrl );
ElError ElLPPrimalADMMCtrlDefault_d( ElLPPrimalADMMCtrl_d* ctrl );

typedef struct {
  float gamma;
  float beta;
  float psi;
  bool print;
} ElLPPrimalIPFLineSearchCtrl_s;

typedef struct {
  double gamma;
  double beta;
  double psi;
  bool print;
} ElLPPrimalIPFLineSearchCtrl_d;

ElError ElLPPrimalIPFLineSearchCtrlDefault_s
( ElLPPrimalIPFLineSearchCtrl_s* ctrl );
ElError ElLPPrimalIPFLineSearchCtrlDefault_d
( ElLPPrimalIPFLineSearchCtrl_d* ctrl );

typedef struct {
  bool initialized;
  float tol;
  ElInt maxIts;
  float centering;
  ElLPPrimalKKTSystem system;

  ElLPPrimalIPFLineSearchCtrl_s lineSearchCtrl;
  bool print;
} ElLPPrimalIPFCtrl_s;

typedef struct {
  bool initialized;
  double tol;
  ElInt maxIts;
  double centering;
  ElLPPrimalKKTSystem system;

  ElLPPrimalIPFLineSearchCtrl_d lineSearchCtrl;
  bool print;
} ElLPPrimalIPFCtrl_d;

ElError ElLPPrimalIPFCtrlDefault_s( ElLPPrimalIPFCtrl_s* ctrl, bool isSparse );
ElError ElLPPrimalIPFCtrlDefault_d( ElLPPrimalIPFCtrl_d* ctrl, bool isSparse );

typedef struct {
  bool initialized;
  float tol;
  ElInt maxIts;
  float maxStepRatio;
  ElLPPrimalKKTSystem system;
  bool print;
} ElLPPrimalMehrotraCtrl_s;

typedef struct {
  bool initialized;
  double tol;
  ElInt maxIts;
  double maxStepRatio;
  ElLPPrimalKKTSystem system;
  bool print;
} ElLPPrimalMehrotraCtrl_d;

ElError ElLPPrimalMehrotraCtrlDefault_s
( ElLPPrimalMehrotraCtrl_s* ctrl, bool isSparse );
ElError ElLPPrimalMehrotraCtrlDefault_d
( ElLPPrimalMehrotraCtrl_d* ctrl, bool isSparse );

typedef struct {
  ElLPApproach approach;  
  ElLPPrimalADMMCtrl_s admmCtrl;
  ElLPPrimalIPFCtrl_s ipfCtrl;
  ElLPPrimalMehrotraCtrl_s mehrotraCtrl;
} ElLPPrimalCtrl_s;
typedef struct {
  ElLPApproach approach;  
  ElLPPrimalADMMCtrl_d admmCtrl;
  ElLPPrimalIPFCtrl_d ipfCtrl;
  ElLPPrimalMehrotraCtrl_d mehrotraCtrl;
} ElLPPrimalCtrl_d;

ElError ElLPPrimalCtrlDefault_s( ElLPPrimalCtrl_s* ctrl, bool isSparse );
ElError ElLPPrimalCtrlDefault_d( ElLPPrimalCtrl_d* ctrl, bool isSparse );

EL_EXPORT ElError ElLPPrimalX_s
( ElConstMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c, 
  ElMatrix_s x,      ElMatrix_s y, 
  ElMatrix_s z, 
  ElLPPrimalCtrl_s ctrl );
EL_EXPORT ElError ElLPPrimalX_d
( ElConstMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c, 
  ElMatrix_d x,      ElMatrix_d y, 
  ElMatrix_d z, 
  ElLPPrimalCtrl_d ctrl );

EL_EXPORT ElError ElLPPrimalXDist_s
( ElConstDistMatrix_s A, 
  ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElDistMatrix_s x,      ElDistMatrix_s y, 
  ElDistMatrix_s z, 
  ElLPPrimalCtrl_s ctrl );
EL_EXPORT ElError ElLPPrimalXDist_d
( ElConstDistMatrix_d A, 
  ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElDistMatrix_d x,      ElDistMatrix_d y, 
  ElDistMatrix_d z, 
  ElLPPrimalCtrl_d ctrl );

EL_EXPORT ElError ElLPPrimalXSparse_s
( ElConstSparseMatrix_s A, 
  ElConstMatrix_s b,       ElConstMatrix_s c, 
  ElMatrix_s x,            ElMatrix_s y, 
  ElMatrix_s z, 
  ElLPPrimalCtrl_s ctrl );
EL_EXPORT ElError ElLPPrimalXSparse_d
( ElConstSparseMatrix_d A, 
  ElConstMatrix_d b,       ElConstMatrix_d c, 
  ElMatrix_d x,            ElMatrix_d y, 
  ElMatrix_d z, 
  ElLPPrimalCtrl_d ctrl );

EL_EXPORT ElError ElLPPrimalXDistSparse_s
( ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s b,     ElConstDistMultiVec_s c, 
  ElDistMultiVec_s x,          ElDistMultiVec_s y, 
  ElDistMultiVec_s z, 
  ElLPPrimalCtrl_s ctrl );
EL_EXPORT ElError ElLPPrimalXDistSparse_d
( ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d b,     ElConstDistMultiVec_d c, 
  ElDistMultiVec_d x,          ElDistMultiVec_d y, 
  ElDistMultiVec_d z,
  ElLPPrimalCtrl_d ctrl );

/* Dual conic form
   --------------- */
EL_EXPORT ElError ElLPDual_s
( ElConstMatrix_s A, ElConstMatrix_s G,
  ElConstMatrix_s b, ElConstMatrix_s c, 
  ElConstMatrix_s h,
  ElMatrix_s x,      ElMatrix_s y, 
  ElMatrix_s z,      ElMatrix_s s );
EL_EXPORT ElError ElLPDual_d
( ElConstMatrix_d A, ElConstMatrix_d G,
  ElConstMatrix_d b, ElConstMatrix_d c, 
  ElConstMatrix_d h,
  ElMatrix_d x,      ElMatrix_d y, 
  ElMatrix_d z,      ElMatrix_d s );

EL_EXPORT ElError ElLPDualDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s G,
  ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElConstDistMatrix_s h,
  ElDistMatrix_s x,      ElDistMatrix_s y, 
  ElDistMatrix_s z,      ElDistMatrix_s s );
EL_EXPORT ElError ElLPDualDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d G,
  ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElConstDistMatrix_d h,
  ElDistMatrix_d x,      ElDistMatrix_d y, 
  ElDistMatrix_d z,      ElDistMatrix_d s );

EL_EXPORT ElError ElLPDualSparse_s
( ElConstSparseMatrix_s A, ElConstSparseMatrix_s G,
  ElConstMatrix_s b,       ElConstMatrix_s c, 
  ElConstMatrix_s h,
  ElMatrix_s x,            ElMatrix_s y, 
  ElMatrix_s z,            ElMatrix_s s );
EL_EXPORT ElError ElLPDualSparse_d
( ElConstSparseMatrix_d A, ElConstSparseMatrix_d G,
  ElConstMatrix_d b,       ElConstMatrix_d c, 
  ElConstMatrix_d h,
  ElMatrix_d x,            ElMatrix_d y, 
  ElMatrix_d z,            ElMatrix_d s );

EL_EXPORT ElError ElLPDualDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistSparseMatrix_s G,
  ElConstDistMultiVec_s b,     ElConstDistMultiVec_s c, 
  ElConstDistMultiVec_s h,
  ElDistMultiVec_s x,          ElDistMultiVec_s y, 
  ElDistMultiVec_s z,          ElDistMultiVec_s s );
EL_EXPORT ElError ElLPDualDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistSparseMatrix_d G,
  ElConstDistMultiVec_d b,     ElConstDistMultiVec_d c, 
  ElConstDistMultiVec_d h,
  ElDistMultiVec_d x,          ElDistMultiVec_d y, 
  ElDistMultiVec_d z,          ElDistMultiVec_d s );

/* Expert versions
   ^^^^^^^^^^^^^^^ */
typedef struct {
  float gamma;
  float beta;
  float psi;
  bool print;
} ElLPDualIPFLineSearchCtrl_s;

typedef struct {
  double gamma;
  double beta;
  double psi;
  bool print;
} ElLPDualIPFLineSearchCtrl_d;

ElError ElLPDualIPFLineSearchCtrlDefault_s( ElLPDualIPFLineSearchCtrl_s* ctrl );
ElError ElLPDualIPFLineSearchCtrlDefault_d( ElLPDualIPFLineSearchCtrl_d* ctrl );

typedef struct {
  bool initialized;
  float tol;
  ElInt maxIts;
  float centering;

  ElLPDualIPFLineSearchCtrl_s lineSearchCtrl;
  bool print;
} ElLPDualIPFCtrl_s;

typedef struct {
  bool initialized;
  double tol;
  ElInt maxIts;
  double centering;

  ElLPDualIPFLineSearchCtrl_d lineSearchCtrl;
  bool print;
} ElLPDualIPFCtrl_d;

ElError ElLPDualIPFCtrlDefault_s( ElLPDualIPFCtrl_s* ctrl );
ElError ElLPDualIPFCtrlDefault_d( ElLPDualIPFCtrl_d* ctrl );

typedef struct {
  bool initialized;
  float tol;
  ElInt maxIts;
  float maxStepRatio;
  bool print;
} ElLPDualMehrotraCtrl_s;

typedef struct {
  bool initialized;
  double tol;
  ElInt maxIts;
  double maxStepRatio;
  bool print;
} ElLPDualMehrotraCtrl_d;

ElError ElLPDualMehrotraCtrlDefault_s( ElLPDualMehrotraCtrl_s* ctrl );
ElError ElLPDualMehrotraCtrlDefault_d( ElLPDualMehrotraCtrl_d* ctrl );

typedef struct {
  ElLPApproach approach;  
  ElLPDualIPFCtrl_s ipfCtrl;
  ElLPDualMehrotraCtrl_s mehrotraCtrl;
} ElLPDualCtrl_s;
typedef struct {
  ElLPApproach approach;  
  ElLPDualIPFCtrl_d ipfCtrl;
  ElLPDualMehrotraCtrl_d mehrotraCtrl;
} ElLPDualCtrl_d;

ElError ElLPDualCtrlDefault_s( ElLPDualCtrl_s* ctrl );
ElError ElLPDualCtrlDefault_d( ElLPDualCtrl_d* ctrl );

EL_EXPORT ElError ElLPDualX_s
( ElConstMatrix_s A, ElConstMatrix_s G,
  ElConstMatrix_s b, ElConstMatrix_s c, 
  ElConstMatrix_s h,
  ElMatrix_s x,      ElMatrix_s y, 
  ElMatrix_s z,      ElMatrix_s s,
  ElLPDualCtrl_s ctrl );
EL_EXPORT ElError ElLPDualX_d
( ElConstMatrix_d A, ElConstMatrix_d G,
  ElConstMatrix_d b, ElConstMatrix_d c, 
  ElConstMatrix_d h,
  ElMatrix_d x,      ElMatrix_d y, 
  ElMatrix_d z,      ElMatrix_d s,
  ElLPDualCtrl_d ctrl );

EL_EXPORT ElError ElLPDualXDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s G,
  ElConstDistMatrix_s b, ElConstDistMatrix_s c, 
  ElConstDistMatrix_s h,
  ElDistMatrix_s x,      ElDistMatrix_s y, 
  ElDistMatrix_s z,      ElDistMatrix_s s,
  ElLPDualCtrl_s ctrl );
EL_EXPORT ElError ElLPDualXDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d G,
  ElConstDistMatrix_d b, ElConstDistMatrix_d c, 
  ElConstDistMatrix_d h,
  ElDistMatrix_d x,      ElDistMatrix_d y, 
  ElDistMatrix_d z,      ElDistMatrix_d s,
  ElLPDualCtrl_d ctrl );

EL_EXPORT ElError ElLPDualXSparse_s
( ElConstSparseMatrix_s A, ElConstSparseMatrix_s G,
  ElConstMatrix_s b,       ElConstMatrix_s c, 
  ElConstMatrix_s h,
  ElMatrix_s x,            ElMatrix_s y, 
  ElMatrix_s z,            ElMatrix_s s,
  ElLPDualCtrl_s ctrl );
EL_EXPORT ElError ElLPDualXSparse_d
( ElConstSparseMatrix_d A, ElConstSparseMatrix_d G,
  ElConstMatrix_d b,       ElConstMatrix_d c, 
  ElConstMatrix_d h,
  ElMatrix_d x,            ElMatrix_d y, 
  ElMatrix_d z,            ElMatrix_d s,
  ElLPDualCtrl_d ctrl );

EL_EXPORT ElError ElLPDualXDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistSparseMatrix_s G,
  ElConstDistMultiVec_s b,     ElConstDistMultiVec_s c, 
  ElConstDistMultiVec_s h,
  ElDistMultiVec_s x,          ElDistMultiVec_s y, 
  ElDistMultiVec_s z,          ElDistMultiVec_s s,
  ElLPDualCtrl_s ctrl );
EL_EXPORT ElError ElLPDualXDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistSparseMatrix_d G,
  ElConstDistMultiVec_d b,     ElConstDistMultiVec_d c, 
  ElConstDistMultiVec_d h,
  ElDistMultiVec_d x,          ElDistMultiVec_d y, 
  ElDistMultiVec_d z,          ElDistMultiVec_d s,
  ElLPDualCtrl_d ctrl );

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
/* Primal conic form
   ----------------- */
EL_EXPORT ElError ElQPPrimal_s
( ElConstMatrix_s Q, ElConstMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s x );
EL_EXPORT ElError ElQPPrimal_d
( ElConstMatrix_d Q, ElConstMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d x );

EL_EXPORT ElError ElQPPrimalDist_s
( ElConstDistMatrix_s Q, ElConstDistMatrix_s A, 
  ElConstDistMatrix_s b, ElConstDistMatrix_s c,
  ElDistMatrix_s x );
EL_EXPORT ElError ElQPPrimalDist_d
( ElConstDistMatrix_d Q, ElConstDistMatrix_d A, 
  ElConstDistMatrix_d b, ElConstDistMatrix_d c,
  ElDistMatrix_d x );

EL_EXPORT ElError ElQPPrimalSparse_s
( ElConstSparseMatrix_s Q, ElConstSparseMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s x );
EL_EXPORT ElError ElQPPrimalSparse_d
( ElConstSparseMatrix_d Q, ElConstSparseMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d x );

EL_EXPORT ElError ElQPPrimalDistSparse_s
( ElConstDistSparseMatrix_s Q, ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s b, ElConstDistMultiVec_s c,
  ElDistMultiVec_s x );
EL_EXPORT ElError ElQPPrimalDistSparse_d
( ElConstDistSparseMatrix_d Q, ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d b, ElConstDistMultiVec_d c,
  ElDistMultiVec_d x );

/* Infeasible Path-Following
   ^^^^^^^^^^^^^^^^^^^^^^^^^ */
EL_EXPORT ElError ElQPPrimalIPF_s
( ElConstMatrix_s Q, ElConstMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s x, ElMatrix_s y, ElMatrix_s z );
EL_EXPORT ElError ElQPPrimalIPF_d
( ElConstMatrix_d Q, ElConstMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d x, ElMatrix_d y, ElMatrix_d z );

EL_EXPORT ElError ElQPPrimalIPFDist_s
( ElConstDistMatrix_s Q, ElConstDistMatrix_s A, 
  ElConstDistMatrix_s b, ElConstDistMatrix_s c,
  ElDistMatrix_s x, ElDistMatrix_s y, ElDistMatrix_s z );
EL_EXPORT ElError ElQPPrimalIPFDist_d
( ElConstDistMatrix_d Q, ElConstDistMatrix_d A, 
  ElConstDistMatrix_d b, ElConstDistMatrix_d c,
  ElDistMatrix_d x, ElDistMatrix_d y, ElDistMatrix_d z );

EL_EXPORT ElError ElQPPrimalIPFSparse_s
( ElConstSparseMatrix_s Q, ElConstSparseMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s x, ElMatrix_s y, ElMatrix_s z );
EL_EXPORT ElError ElQPPrimalIPFSparse_d
( ElConstSparseMatrix_d Q, ElConstSparseMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d x, ElMatrix_d y, ElMatrix_d z );

EL_EXPORT ElError ElQPPrimalIPFDistSparse_s
( ElConstDistSparseMatrix_s Q, ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s b, ElConstDistMultiVec_s c,
  ElDistMultiVec_s x, ElDistMultiVec_s y, ElDistMultiVec_s z );
EL_EXPORT ElError ElQPPrimalIPFDistSparse_d
( ElConstDistSparseMatrix_d Q, ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d b, ElConstDistMultiVec_d c,
  ElDistMultiVec_d x, ElDistMultiVec_d y, ElDistMultiVec_d z );

/* Mehrotra Predictor-Corrector
   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
EL_EXPORT ElError ElQPPrimalMehrotra_s
( ElConstMatrix_s Q, ElConstMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s x, ElMatrix_s y, ElMatrix_s z );
EL_EXPORT ElError ElQPPrimalMehrotra_d
( ElConstMatrix_d Q, ElConstMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d x, ElMatrix_d y, ElMatrix_d z );

EL_EXPORT ElError ElQPPrimalMehrotraDist_s
( ElConstDistMatrix_s Q, ElConstDistMatrix_s A, 
  ElConstDistMatrix_s b, ElConstDistMatrix_s c,
  ElDistMatrix_s x, ElDistMatrix_s y, ElDistMatrix_s z );
EL_EXPORT ElError ElQPPrimalMehrotraDist_d
( ElConstDistMatrix_d Q, ElConstDistMatrix_d A, 
  ElConstDistMatrix_d b, ElConstDistMatrix_d c,
  ElDistMatrix_d x, ElDistMatrix_d y, ElDistMatrix_d z );

EL_EXPORT ElError ElQPPrimalMehrotraSparse_s
( ElConstSparseMatrix_s Q, ElConstSparseMatrix_s A, 
  ElConstMatrix_s b, ElConstMatrix_s c,
  ElMatrix_s x, ElMatrix_s y, ElMatrix_s z );
EL_EXPORT ElError ElQPPrimalMehrotraSparse_d
( ElConstSparseMatrix_d Q, ElConstSparseMatrix_d A, 
  ElConstMatrix_d b, ElConstMatrix_d c,
  ElMatrix_d x, ElMatrix_d y, ElMatrix_d z );

EL_EXPORT ElError ElQPPrimalMehrotraDistSparse_s
( ElConstDistSparseMatrix_s Q, ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s b, ElConstDistMultiVec_s c,
  ElDistMultiVec_s x, ElDistMultiVec_s y, ElDistMultiVec_s z );
EL_EXPORT ElError ElQPPrimalMehrotraDistSparse_d
( ElConstDistSparseMatrix_d Q, ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d b, ElConstDistMultiVec_d c,
  ElDistMultiVec_d x, ElDistMultiVec_d y, ElDistMultiVec_d z );

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
