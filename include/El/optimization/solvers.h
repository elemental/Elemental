/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_OPTIMIZATION_SOLVERS_C_H
#define EL_OPTIMIZATION_SOLVERS_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Linear programs
   =============== */
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
  ElMatrix_s X, ElInt* numIts );
EL_EXPORT ElError ElQPBoxADMM_d
( ElConstMatrix_d Q, ElConstMatrix_d C, double lb, double ub, 
  ElMatrix_d X, ElInt* numIts );

EL_EXPORT ElError ElQPBoxADMMDist_s
( ElConstDistMatrix_s Q, ElConstDistMatrix_s C, float lb, float ub, 
  ElDistMatrix_s X, ElInt* numIts );
EL_EXPORT ElError ElQPBoxADMMDist_d
( ElConstDistMatrix_d Q, ElConstDistMatrix_d C, double lb, double ub, 
  ElDistMatrix_d X, ElInt* numIts );

/* TODO: Expert versions */

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_OPTIMIZATION_SOLVERS_C_H */
