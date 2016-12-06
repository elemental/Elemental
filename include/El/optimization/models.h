/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_MODELS_C_H
#define EL_OPTIMIZATION_MODELS_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Basis pursuit
   ============= */
EL_EXPORT ElError ElBP_s
( ElConstMatrix_s A,
  ElConstMatrix_s b,
  ElMatrix_s x );
EL_EXPORT ElError ElBP_d
( ElConstMatrix_d A,
  ElConstMatrix_d b,
  ElMatrix_d x );

EL_EXPORT ElError ElBPDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b,
  ElDistMatrix_s x );
EL_EXPORT ElError ElBPDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b,
  ElDistMatrix_d x );

EL_EXPORT ElError ElBPSparse_s
( ElConstSparseMatrix_s A,
  ElConstMatrix_s b,
  ElMatrix_s x );
EL_EXPORT ElError ElBPSparse_d
( ElConstSparseMatrix_d A,
  ElConstMatrix_d b,
  ElMatrix_d x );

EL_EXPORT ElError ElBPDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b,
  ElDistMultiVec_s x );
EL_EXPORT ElError ElBPDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b,
  ElDistMultiVec_d x );

/* Expert verions
   -------------- */
typedef struct {
  float rho;
  float alpha;
  ElInt maxIter;
  float absTol;
  float relTol;
  bool usePinv;
  float pinvTol;
  bool progress;
} ElBPADMMCtrl_s;

typedef struct {
  double rho;
  double alpha;
  ElInt maxIter;
  double absTol;
  double relTol;
  bool usePinv;
  double pinvTol;
  bool progress;
} ElBPADMMCtrl_d;

EL_EXPORT ElError ElBPADMMCtrlDefault_s( ElBPADMMCtrl_s* ctrl );
EL_EXPORT ElError ElBPADMMCtrlDefault_d( ElBPADMMCtrl_d* ctrl );

typedef struct {
  bool useIPM;
  bool useSOCP;
  ElBPADMMCtrl_s admmCtrl;
  ElLPDirectCtrl_s lpIPMCtrl;
  ElSOCPDirectCtrl_s socpIPMCtrl;
} ElBPCtrl_s;

typedef struct {
  bool useIPM;
  bool useSOCP;
  ElBPADMMCtrl_d admmCtrl;
  ElLPDirectCtrl_d lpIPMCtrl;
  ElSOCPDirectCtrl_d socpIPMCtrl;
} ElBPCtrl_d;

typedef struct {
  ElSOCPDirectCtrl_s ipmCtrl;
} ElBPCtrl_c;

typedef struct {
  ElSOCPDirectCtrl_d ipmCtrl;
} ElBPCtrl_z;

EL_EXPORT ElError ElBPCtrlDefault_s( ElBPCtrl_s* ctrl, bool isSparse );
EL_EXPORT ElError ElBPCtrlDefault_d( ElBPCtrl_d* ctrl, bool isSparse );
EL_EXPORT ElError ElBPCtrlDefault_c( ElBPCtrl_c* ctrl );
EL_EXPORT ElError ElBPCtrlDefault_z( ElBPCtrl_z* ctrl );

EL_EXPORT ElError ElBPX_s
( ElConstMatrix_s A,
  ElConstMatrix_s b,
  ElMatrix_s x,
  ElBPCtrl_s ctrl );
EL_EXPORT ElError ElBPX_d
( ElConstMatrix_d A,
  ElConstMatrix_d b,
  ElMatrix_d x,
  ElBPCtrl_d ctrl );
EL_EXPORT ElError ElBPX_c
( ElConstMatrix_c A,
  ElConstMatrix_c b,
  ElMatrix_c x,
  ElBPCtrl_c ctrl );
EL_EXPORT ElError ElBPX_z
( ElConstMatrix_z A,
  ElConstMatrix_z b,
  ElMatrix_z x,
  ElBPCtrl_z ctrl );

EL_EXPORT ElError ElBPXDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b,
  ElDistMatrix_s x,
  ElBPCtrl_s ctrl );
EL_EXPORT ElError ElBPXDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b,
  ElDistMatrix_d x,
  ElBPCtrl_d ctrl );
EL_EXPORT ElError ElBPXDist_c
( ElConstDistMatrix_c A,
  ElConstDistMatrix_c b,
  ElDistMatrix_c x,
  ElBPCtrl_c ctrl );
EL_EXPORT ElError ElBPXDist_z
( ElConstDistMatrix_z A,
  ElConstDistMatrix_z b,
  ElDistMatrix_z x,
  ElBPCtrl_z ctrl );

EL_EXPORT ElError ElBPXSparse_s
( ElConstSparseMatrix_s A,
  ElConstMatrix_s b,
  ElMatrix_s x,
  ElBPCtrl_s ctrl );
EL_EXPORT ElError ElBPXSparse_d
( ElConstSparseMatrix_d A,
  ElConstMatrix_d b,
  ElMatrix_d x,
  ElBPCtrl_d ctrl );
EL_EXPORT ElError ElBPXSparse_c
( ElConstSparseMatrix_c A,
  ElConstMatrix_c b,
  ElMatrix_c x,
  ElBPCtrl_c ctrl );
EL_EXPORT ElError ElBPXSparse_z
( ElConstSparseMatrix_z A,
  ElConstMatrix_z b,
  ElMatrix_z x,
  ElBPCtrl_z ctrl );

EL_EXPORT ElError ElBPXDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b,
  ElDistMultiVec_s x,
  ElBPCtrl_s ctrl );
EL_EXPORT ElError ElBPXDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b,
  ElDistMultiVec_d x,
  ElBPCtrl_d ctrl );
EL_EXPORT ElError ElBPXDistSparse_c
( ElConstDistSparseMatrix_c A,
  ElConstDistMultiVec_c b,
  ElDistMultiVec_c x,
  ElBPCtrl_c ctrl );
EL_EXPORT ElError ElBPXDistSparse_z
( ElConstDistSparseMatrix_z A,
  ElConstDistMultiVec_z b,
  ElDistMultiVec_z x,
  ElBPCtrl_z ctrl );

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

/* Robust least squares
   ==================== */
EL_EXPORT ElError ElRLS_s
( ElConstMatrix_s A,
  ElConstMatrix_s b, float rho,
  ElMatrix_s x );
EL_EXPORT ElError ElRLS_d
( ElConstMatrix_d A,
  ElConstMatrix_d b, double rho,
  ElMatrix_d x );

EL_EXPORT ElError ElRLSSparse_s
( ElConstSparseMatrix_s A,
  ElConstMatrix_s b, float rho,
  ElMatrix_s x );
EL_EXPORT ElError ElRLSSparse_d
( ElConstSparseMatrix_d A,
  ElConstMatrix_d b, double rho,
  ElMatrix_d x );

EL_EXPORT ElError ElRLSDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b, float rho,
  ElDistMatrix_s x );
EL_EXPORT ElError ElRLSDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b, double rho,
  ElDistMatrix_d x );

EL_EXPORT ElError ElRLSDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b, float rho,
  ElDistMultiVec_s x );
EL_EXPORT ElError ElRLSDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b, double rho,
  ElDistMultiVec_d x );

/* Expert versions
   --------------- */
EL_EXPORT ElError ElRLSX_s
( ElConstMatrix_s A,
  ElConstMatrix_s b, float rho,
  ElMatrix_s x,
  ElSOCPAffineCtrl_s ctrl );
EL_EXPORT ElError ElRLSX_d
( ElConstMatrix_d A,
  ElConstMatrix_d b, double rho,
  ElMatrix_d x,
  ElSOCPAffineCtrl_d ctrl );

EL_EXPORT ElError ElRLSXSparse_s
( ElConstSparseMatrix_s A,
  ElConstMatrix_s b, float rho,
  ElMatrix_s x,
  ElSOCPAffineCtrl_s ctrl );
EL_EXPORT ElError ElRLSXSparse_d
( ElConstSparseMatrix_d A,
  ElConstMatrix_d b, double rho,
  ElMatrix_d x,
  ElSOCPAffineCtrl_d ctrl );

EL_EXPORT ElError ElRLSXDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b, float rho,
  ElDistMatrix_s x,
  ElSOCPAffineCtrl_s ctrl );
EL_EXPORT ElError ElRLSXDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b, double rho,
  ElDistMatrix_d x,
  ElSOCPAffineCtrl_d ctrl );

EL_EXPORT ElError ElRLSXDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b, float rho,
  ElDistMultiVec_s x,
  ElSOCPAffineCtrl_s ctrl );
EL_EXPORT ElError ElRLSXDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b, double rho,
  ElDistMultiVec_d x,
  ElSOCPAffineCtrl_d ctrl );

/* Robust non-negative least squares
   ================================= */
EL_EXPORT ElError ElRNNLS_s
( ElConstMatrix_s A,
  ElConstMatrix_s b, float rho,
  ElMatrix_s x );
EL_EXPORT ElError ElRNNLS_d
( ElConstMatrix_d A,
  ElConstMatrix_d b, double rho,
  ElMatrix_d x );

EL_EXPORT ElError ElRNNLSSparse_s
( ElConstSparseMatrix_s A,
  ElConstMatrix_s b, float rho,
  ElMatrix_s x );
EL_EXPORT ElError ElRNNLSSparse_d
( ElConstSparseMatrix_d A,
  ElConstMatrix_d b, double rho,
  ElMatrix_d x );

EL_EXPORT ElError ElRNNLSDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b, float rho,
  ElDistMatrix_s x );
EL_EXPORT ElError ElRNNLSDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b, double rho,
  ElDistMatrix_d x );

EL_EXPORT ElError ElRNNLSDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b, float rho,
  ElDistMultiVec_s x );
EL_EXPORT ElError ElRNNLSDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b, double rho,
  ElDistMultiVec_d x );

/* Expert versions
   --------------- */
EL_EXPORT ElError ElRNNLSX_s
( ElConstMatrix_s A,
  ElConstMatrix_s b, float rho,
  ElMatrix_s x,
  ElSOCPAffineCtrl_s ctrl );
EL_EXPORT ElError ElRNNLSX_d
( ElConstMatrix_d A,
  ElConstMatrix_d b, double rho,
  ElMatrix_d x,
  ElSOCPAffineCtrl_d ctrl );

EL_EXPORT ElError ElRNNLSXSparse_s
( ElConstSparseMatrix_s A,
  ElConstMatrix_s b, float rho,
  ElMatrix_s x,
  ElSOCPAffineCtrl_s ctrl );
EL_EXPORT ElError ElRNNLSXSparse_d
( ElConstSparseMatrix_d A,
  ElConstMatrix_d b, double rho,
  ElMatrix_d x,
  ElSOCPAffineCtrl_d ctrl );

EL_EXPORT ElError ElRNNLSXDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b, float rho,
  ElDistMatrix_s x,
  ElSOCPAffineCtrl_s ctrl );
EL_EXPORT ElError ElRNNLSXDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b, double rho,
  ElDistMatrix_d x,
  ElSOCPAffineCtrl_d ctrl );

EL_EXPORT ElError ElRNNLSXDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b, float rho,
  ElDistMultiVec_s x,
  ElSOCPAffineCtrl_s ctrl );
EL_EXPORT ElError ElRNNLSXDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b, double rho,
  ElDistMultiVec_d x,
  ElSOCPAffineCtrl_d ctrl );

/* Non-negative least squares
   ========================== */
EL_EXPORT ElError ElNNLS_s
( ElConstMatrix_s A,
  ElConstMatrix_s b,
  ElMatrix_s x );
EL_EXPORT ElError ElNNLS_d
( ElConstMatrix_d A,
  ElConstMatrix_d b,
  ElMatrix_d x );

EL_EXPORT ElError ElNNLSDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b,
  ElDistMatrix_s x );
EL_EXPORT ElError ElNNLSDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b,
  ElDistMatrix_d x );

EL_EXPORT ElError ElNNLSDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b,
  ElDistMatrix_s x );
EL_EXPORT ElError ElNNLSDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b,
  ElDistMatrix_d x );

EL_EXPORT ElError ElNNLSDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b,
  ElDistMultiVec_s x );
EL_EXPORT ElError ElNNLSDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b,
  ElDistMultiVec_d x );

/* Expert versions
   --------------- */
typedef enum {
  EL_NNLS_ADMM,
  EL_NNLS_QP,
  EL_NNLS_SOCP
} ElNNLSApproach;

typedef struct {
  ElNNLSApproach approach;
  ElADMMCtrl_s admmCtrl;
  ElQPDirectCtrl_s qpCtrl;
  ElSOCPAffineCtrl_s socpCtrl;
} ElNNLSCtrl_s;

typedef struct {
  ElNNLSApproach approach;
  ElADMMCtrl_d admmCtrl;
  ElQPDirectCtrl_d qpCtrl;
  ElSOCPAffineCtrl_d socpCtrl;
} ElNNLSCtrl_d;

EL_EXPORT ElError ElNNLSCtrlDefault_s( ElNNLSCtrl_s* ctrl );
EL_EXPORT ElError ElNNLSCtrlDefault_d( ElNNLSCtrl_d* ctrl );

EL_EXPORT ElError ElNNLSX_s
( ElConstMatrix_s A,
  ElConstMatrix_s b,
  ElMatrix_s x,
  ElNNLSCtrl_s ctrl );
EL_EXPORT ElError ElNNLSX_d
( ElConstMatrix_d A,
  ElConstMatrix_d b,
  ElMatrix_d x,
  ElNNLSCtrl_d ctrl );

EL_EXPORT ElError ElNNLSXDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b,
  ElDistMatrix_s x,
  ElNNLSCtrl_s ctrl );
EL_EXPORT ElError ElNNLSXDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b,
  ElDistMatrix_d x,
  ElNNLSCtrl_d ctrl );

EL_EXPORT ElError ElNNLSXDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b,
  ElDistMatrix_s x,
  ElNNLSCtrl_s ctrl );
EL_EXPORT ElError ElNNLSXDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b,
  ElDistMatrix_d x,
  ElNNLSCtrl_d ctrl );

EL_EXPORT ElError ElNNLSXDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b,
  ElDistMultiVec_s x,
  ElNNLSCtrl_s ctrl );
EL_EXPORT ElError ElNNLSXDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b,
  ElDistMultiVec_d x,
  ElNNLSCtrl_d ctrl );

/* Non-negative matrix factorization
   ================================= */
EL_EXPORT ElError ElNMF_s( ElConstMatrix_s A, ElMatrix_s X, ElMatrix_s Y );
EL_EXPORT ElError ElNMF_d( ElConstMatrix_d A, ElMatrix_d X, ElMatrix_d Y );

EL_EXPORT ElError ElNMFDist_s
( ElConstDistMatrix_s A,
  ElDistMatrix_s X,
  ElDistMatrix_s Y );
EL_EXPORT ElError ElNMFDist_d
( ElConstDistMatrix_d A,
  ElDistMatrix_d X,
  ElDistMatrix_d Y );

/* Expert versions */
typedef struct {
  ElNNLSCtrl_s nnlsCtrl;
  ElInt maxIter;
} ElNMFCtrl_s;

typedef struct {
  ElNNLSCtrl_d nnlsCtrl;
  ElInt maxIter;
} ElNMFCtrl_d;

EL_EXPORT ElError ElNMFCtrlDefault_s( ElNMFCtrl_s* ctrl );
EL_EXPORT ElError ElNMFCtrlDefault_d( ElNMFCtrl_d* ctrl );

EL_EXPORT ElError ElNMFX_s
( ElConstMatrix_s A,
  ElMatrix_s X,
  ElMatrix_s Y,
  ElNMFCtrl_s ctrl );
EL_EXPORT ElError ElNMFX_d
( ElConstMatrix_d A,
  ElMatrix_d X,
  ElMatrix_d Y,
  ElNMFCtrl_d ctrl );

EL_EXPORT ElError ElNMFXDist_s
( ElConstDistMatrix_s A,
  ElDistMatrix_s X,
  ElDistMatrix_s Y,
  ElNMFCtrl_s ctrl );
EL_EXPORT ElError ElNMFXDist_d
( ElConstDistMatrix_d A,
  ElDistMatrix_d X,
  ElDistMatrix_d Y,
  ElNMFCtrl_d ctrl );

/* Basis pursuit denoising
   ======================= */
EL_EXPORT ElError ElBPDN_s
( ElConstMatrix_s A,
  ElConstMatrix_s b,
  float lambda,
  ElMatrix_s x );
EL_EXPORT ElError ElBPDN_d
( ElConstMatrix_d A,
  ElConstMatrix_d b,
  double lambda,
  ElMatrix_d x );

EL_EXPORT ElError ElBPDNDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b,
  float lambda,
  ElDistMatrix_s x );
EL_EXPORT ElError ElBPDNDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b,
  double lambda,
  ElDistMatrix_d x );

EL_EXPORT ElError ElBPDNSparse_s
( ElConstSparseMatrix_s A,
  ElConstMatrix_s b,
  float lambda,
  ElMatrix_s x );
EL_EXPORT ElError ElBPDNSparse_d
( ElConstSparseMatrix_d A,
  ElConstMatrix_d b,
  double lambda,
  ElMatrix_d x );

EL_EXPORT ElError ElBPDNDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b,
  float lambda,
  ElDistMultiVec_s x );
EL_EXPORT ElError ElBPDNDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b,
  double lambda,
  ElDistMultiVec_d x );

/* Expert versions
   --------------- */
typedef struct {
  float rho;
  float alpha;
  ElInt maxIter;
  float absTol;
  float relTol;
  bool inv;
  bool progress;
} ElBPDNADMMCtrl_s;

typedef struct {
  double rho;
  double alpha;
  ElInt maxIter;
  double absTol;
  double relTol;
  bool inv;
  bool progress;
} ElBPDNADMMCtrl_d;

EL_EXPORT ElError ElBPDNADMMCtrlDefault_s( ElBPDNADMMCtrl_s* ctrl );
EL_EXPORT ElError ElBPDNADMMCtrlDefault_d( ElBPDNADMMCtrl_d* ctrl );

typedef struct {
  bool useIPM;
  ElBPDNADMMCtrl_s admmCtrl;
  ElQPAffineCtrl_s ipmCtrl;
} ElBPDNCtrl_s;

typedef struct {
  bool useIPM;
  ElBPDNADMMCtrl_d admmCtrl;
  ElQPAffineCtrl_d ipmCtrl;
} ElBPDNCtrl_d;

EL_EXPORT ElError ElBPDNCtrlDefault_s( ElBPDNCtrl_s* ctrl );
EL_EXPORT ElError ElBPDNCtrlDefault_d( ElBPDNCtrl_d* ctrl );

EL_EXPORT ElError ElBPDNX_s
( ElConstMatrix_s A,
  ElConstMatrix_s b,
  float lambda,
  ElMatrix_s x,
  ElBPDNCtrl_s ctrl );
EL_EXPORT ElError ElBPDNX_d
( ElConstMatrix_d A,
  ElConstMatrix_d b,
  double lambda,
  ElMatrix_d x,
  ElBPDNCtrl_d ctrl );

EL_EXPORT ElError ElBPDNXDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b,
  float lambda,
  ElDistMatrix_s x,
  ElBPDNCtrl_s ctrl );
EL_EXPORT ElError ElBPDNXDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b,
  double lambda,
  ElDistMatrix_d x,
  ElBPDNCtrl_d ctrl );

EL_EXPORT ElError ElBPDNXSparse_s
( ElConstSparseMatrix_s A,
  ElConstMatrix_s b,
  float lambda,
  ElMatrix_s x,
  ElBPDNCtrl_s ctrl );
EL_EXPORT ElError ElBPDNXSparse_d
( ElConstSparseMatrix_d A,
  ElConstMatrix_d b,
  double lambda,
  ElMatrix_d x,
  ElBPDNCtrl_d ctrl );

EL_EXPORT ElError ElBPDNXDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b,
  float lambda,
  ElDistMultiVec_s x,
  ElBPDNCtrl_s ctrl );
EL_EXPORT ElError ElBPDNXDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b,
  double lambda,
  ElDistMultiVec_d x,
  ElBPDNCtrl_d ctrl );

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

/* Expert interface
   ---------------- */
typedef struct {
  bool useALM;
  bool usePivQR;
  bool progress;
  ElInt numPivSteps;
  ElInt maxIts;
  float tau;
  float beta;
  float rho;
  float tol;
} ElRPCACtrl_s;

typedef struct {
  bool useALM;
  bool usePivQR;
  bool progress;
  ElInt numPivSteps;
  ElInt maxIts;
  double tau;
  double beta;
  double rho;
  double tol;
} ElRPCACtrl_d;

EL_EXPORT ElError ElRPCACtrlDefault_s( ElRPCACtrl_s* ctrl );
EL_EXPORT ElError ElRPCACtrlDefault_d( ElRPCACtrl_d* ctrl );

EL_EXPORT ElError ElRPCAX_s
( ElConstMatrix_s M, ElMatrix_s L, ElMatrix_s S, ElRPCACtrl_s ctrl );
EL_EXPORT ElError ElRPCAX_d
( ElConstMatrix_d M, ElMatrix_d L, ElMatrix_d S, ElRPCACtrl_d ctrl );
EL_EXPORT ElError ElRPCAX_c
( ElConstMatrix_c M, ElMatrix_c L, ElMatrix_c S, ElRPCACtrl_s ctrl );
EL_EXPORT ElError ElRPCAX_z
( ElConstMatrix_z M, ElMatrix_z L, ElMatrix_z S, ElRPCACtrl_d ctrl );

EL_EXPORT ElError ElRPCAXDist_s
( ElConstDistMatrix_s M, ElDistMatrix_s L, ElDistMatrix_s S,
  ElRPCACtrl_s ctrl );
EL_EXPORT ElError ElRPCAXDist_d
( ElConstDistMatrix_d M, ElDistMatrix_d L, ElDistMatrix_d S,
  ElRPCACtrl_d ctrl );
EL_EXPORT ElError ElRPCAXDist_c
( ElConstDistMatrix_c M, ElDistMatrix_c L, ElDistMatrix_c S,
  ElRPCACtrl_s ctrl );
EL_EXPORT ElError ElRPCAXDist_z
( ElConstDistMatrix_z M, ElDistMatrix_z L, ElDistMatrix_z S,
  ElRPCACtrl_d ctrl );

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

/* Expert versions
   --------------- */
typedef struct {
  float rho;
  float alpha;
  ElInt maxIter;
  float absTol;
  float relTol;
  bool progress;
} ElSparseInvCovCtrl_s;

typedef struct {
  double rho;
  double alpha;
  ElInt maxIter;
  double absTol;
  double relTol;
  bool progress;
} ElSparseInvCovCtrl_d;

EL_EXPORT ElError ElSparseInvCovCtrlDefault_s( ElSparseInvCovCtrl_s* ctrl );
EL_EXPORT ElError ElSparseInvCovCtrlDefault_d( ElSparseInvCovCtrl_d* ctrl );

EL_EXPORT ElError ElSparseInvCovX_s
( ElConstMatrix_s D, float lambda, ElMatrix_s Z,
  ElSparseInvCovCtrl_s ctrl, ElInt* numIts );
EL_EXPORT ElError ElSparseInvCovX_d
( ElConstMatrix_d D, double lambda, ElMatrix_d Z,
  ElSparseInvCovCtrl_d ctrl, ElInt* numIts );
EL_EXPORT ElError ElSparseInvCovX_c
( ElConstMatrix_c D, float lambda, ElMatrix_c Z,
  ElSparseInvCovCtrl_s ctrl, ElInt* numIts );
EL_EXPORT ElError ElSparseInvCovX_z
( ElConstMatrix_z D, double lambda, ElMatrix_z Z,
  ElSparseInvCovCtrl_d ctrl, ElInt* numIts );

EL_EXPORT ElError ElSparseInvCovXDist_s
( ElConstDistMatrix_s D, float lambda, ElDistMatrix_s Z,
  ElSparseInvCovCtrl_s ctrl, ElInt* numIts );
EL_EXPORT ElError ElSparseInvCovXDist_d
( ElConstDistMatrix_d D, double lambda, ElDistMatrix_d Z,
  ElSparseInvCovCtrl_d ctrl, ElInt* numIts );
EL_EXPORT ElError ElSparseInvCovXDist_c
( ElConstDistMatrix_c D, float lambda, ElDistMatrix_c Z,
  ElSparseInvCovCtrl_s ctrl, ElInt* numIts );
EL_EXPORT ElError ElSparseInvCovXDist_z
( ElConstDistMatrix_z D, double lambda, ElDistMatrix_z Z,
  ElSparseInvCovCtrl_d ctrl, ElInt* numIts );

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
typedef struct
{
  ElQPAffineCtrl_s ipmCtrl;
} ElSVMCtrl_s;

typedef struct
{
  ElQPAffineCtrl_d ipmCtrl;
} ElSVMCtrl_d;

EL_EXPORT ElError ElSVMCtrlDefault_s( ElSVMCtrl_s* ctrl );
EL_EXPORT ElError ElSVMCtrlDefault_d( ElSVMCtrl_d* ctrl );

EL_EXPORT ElError ElSVMX_s
( ElConstMatrix_s A, ElConstMatrix_s d, float lambda,
  ElMatrix_s x, ElSVMCtrl_s ctrl );
EL_EXPORT ElError ElSVMX_d
( ElConstMatrix_d A, ElConstMatrix_d d, double lambda,
  ElMatrix_d x, ElSVMCtrl_d ctrl );

EL_EXPORT ElError ElSVMXDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s d, float lambda,
  ElDistMatrix_s x, ElSVMCtrl_s ctrl );
EL_EXPORT ElError ElSVMXDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d d, double lambda,
  ElDistMatrix_d x, ElSVMCtrl_d ctrl );

EL_EXPORT ElError ElSVMXSparse_s
( ElConstSparseMatrix_s A, ElConstMatrix_s d, float lambda,
  ElMatrix_s x, ElSVMCtrl_s ctrl );
EL_EXPORT ElError ElSVMXSparse_d
( ElConstSparseMatrix_d A, ElConstMatrix_d d, double lambda,
  ElMatrix_d x, ElSVMCtrl_d ctrl );

EL_EXPORT ElError ElSVMXDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s d, float lambda,
  ElDistMultiVec_s x, ElSVMCtrl_s ctrl );
EL_EXPORT ElError ElSVMXDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d d, double lambda,
  ElDistMultiVec_d x, ElSVMCtrl_d ctrl );

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

/* Expert versions
   --------------- */
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

/* Long-only Portfolio
   =================== */
EL_EXPORT ElError ElLongOnlyPortfolioSparse_s
( ElConstMatrix_s d,
  ElConstSparseMatrix_s F,
  ElConstMatrix_s c,
  float gamma,
  ElMatrix_s x );
EL_EXPORT ElError ElLongOnlyPortfolioSparse_d
( ElConstMatrix_d d,
  ElConstSparseMatrix_d F,
  ElConstMatrix_d c,
  double gamma,
  ElMatrix_d x );

EL_EXPORT ElError ElLongOnlyPortfolioDistSparse_s
( ElConstDistMultiVec_s d,
  ElConstDistSparseMatrix_s F,
  ElConstDistMultiVec_s c,
  float gamma,
  ElDistMultiVec_s x );
EL_EXPORT ElError ElLongOnlyPortfolioDistSparse_d
( ElConstDistMultiVec_d d,
  ElConstDistSparseMatrix_d F,
  ElConstDistMultiVec_d c,
  double gamma,
  ElDistMultiVec_d x );

/* Expert versions
   --------------- */
EL_EXPORT ElError ElLongOnlyPortfolioXSparse_s
( ElConstMatrix_s d,
  ElConstSparseMatrix_s F,
  ElConstMatrix_s c,
  float gamma,
  ElMatrix_s x,
  ElSOCPAffineCtrl_s ctrl );
EL_EXPORT ElError ElLongOnlyPortfolioXSparse_d
( ElConstMatrix_d d,
  ElConstSparseMatrix_d F,
  ElConstMatrix_d c,
  double gamma,
  ElMatrix_d x,
  ElSOCPAffineCtrl_d ctrl );

EL_EXPORT ElError ElLongOnlyPortfolioXDistSparse_s
( ElConstDistMultiVec_s d,
  ElConstDistSparseMatrix_s F,
  ElConstDistMultiVec_s c,
  float gamma,
  ElDistMultiVec_s x,
  ElSOCPAffineCtrl_s ctrl );
EL_EXPORT ElError ElLongOnlyPortfolioXDistSparse_d
( ElConstDistMultiVec_d d,
  ElConstDistSparseMatrix_d F,
  ElConstDistMultiVec_d c,
  double gamma,
  ElDistMultiVec_d x,
  ElSOCPAffineCtrl_d ctrl );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_OPTIMIZATION_MODELS_C_H */
