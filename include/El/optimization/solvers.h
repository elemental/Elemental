/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_SOLVERS_C_H
#define EL_OPTIMIZATION_SOLVERS_C_H

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  EL_FULL_KKT,
  EL_AUGMENTED_KKT,
  EL_NORMAL_KKT
} ElKKTSystem;

/* Mehrotra Predictor-Corrector IPM
   ================================ */
/* See the C++ structure for documentation of the members */
typedef struct {
  bool primalInit, dualInit;
  float minTol;
  float targetTol;
  ElInt maxIts;
  float maxStepRatio;
  ElKKTSystem system;
  bool mehrotra;
  float (*centralityRule)(float,float,float,float);
  bool standardInitShift;
  float balanceTol;
  bool forceSameStep;
  ElRegSolveCtrl_s solveCtrl;
  bool resolveReg;
  bool outerEquil;
  ElInt basisSize;
  bool print;
  bool time;
  float wSafeMaxNorm;
  float wMaxLimit;
  float ruizEquilTol;
  ElInt ruizMaxIter;
  float diagEquilTol;
  bool checkResiduals;

  float reg0Tmp;
  float reg1Tmp;
  float reg2Tmp;
  float reg0Perm;
  float reg1Perm;
  float reg2Perm;
} ElMehrotraCtrl_s;

typedef struct {
  bool primalInit, dualInit;
  double minTol;
  double targetTol;
  ElInt maxIts;
  double maxStepRatio;
  ElKKTSystem system;
  bool mehrotra;
  double (*centralityRule)(double,double,double,double);
  bool standardInitShift;
  double balanceTol;
  bool forceSameStep;
  ElRegSolveCtrl_d solveCtrl;
  bool resolveReg;
  bool outerEquil;
  ElInt basisSize;
  bool print;
  bool time;
  double wSafeMaxNorm;
  double wMaxLimit;
  double ruizEquilTol;
  ElInt ruizMaxIter;
  double diagEquilTol;
  bool checkResiduals;

  double reg0Tmp;
  double reg1Tmp;
  double reg2Tmp;
  double reg0Perm;
  double reg1Perm;
  double reg2Perm;
} ElMehrotraCtrl_d;

EL_EXPORT ElError ElMehrotraCtrlDefault_s( ElMehrotraCtrl_s* ctrl );
EL_EXPORT ElError ElMehrotraCtrlDefault_d( ElMehrotraCtrl_d* ctrl );

/* Alternating Direction Method of Multipliers
   =========================================== */
typedef struct {
  float rho; 
  float alpha;
  ElInt maxIter;
  float absTol;
  float relTol;
  bool inv;
  bool print;
} ElADMMCtrl_s;

typedef struct {
  double rho; 
  double alpha;
  ElInt maxIter;
  double absTol;
  double relTol;
  bool inv;
  bool print;
} ElADMMCtrl_d;

EL_EXPORT ElError ElADMMCtrlDefault_s( ElADMMCtrl_s* ctrl );
EL_EXPORT ElError ElADMMCtrlDefault_d( ElADMMCtrl_d* ctrl );

/* Linear programs
   =============== */
typedef enum {
  EL_LP_ADMM,
  EL_LP_MEHROTRA
} ElLPApproach;

/* Direct conic form
   ----------------- */
EL_EXPORT ElError ElLPDirect_s
( ElConstMatrix_s A,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z );
EL_EXPORT ElError ElLPDirect_d
( ElConstMatrix_d A,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z );

EL_EXPORT ElError ElLPDirectDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b,
  ElConstDistMatrix_s c,
  ElDistMatrix_s x,
  ElDistMatrix_s y,
  ElDistMatrix_s z );
EL_EXPORT ElError ElLPDirectDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b,
  ElConstDistMatrix_d c,
  ElDistMatrix_d x,
  ElDistMatrix_d y,
  ElDistMatrix_d z );

EL_EXPORT ElError ElLPDirectSparse_s
( ElConstSparseMatrix_s A,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z );
EL_EXPORT ElError ElLPDirectSparse_d
( ElConstSparseMatrix_d A,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z );

EL_EXPORT ElError ElLPDirectDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b,
  ElConstDistMultiVec_s c,
  ElDistMultiVec_s x,
  ElDistMultiVec_s y,
  ElDistMultiVec_s z );
EL_EXPORT ElError ElLPDirectDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b,
  ElConstDistMultiVec_d c,
  ElDistMultiVec_d x,
  ElDistMultiVec_d y,
  ElDistMultiVec_d z );

/* Expert versions
   ^^^^^^^^^^^^^^^ */
typedef struct {
  ElLPApproach approach; 
  ElADMMCtrl_s admmCtrl;
  ElMehrotraCtrl_s mehrotraCtrl;
} ElLPDirectCtrl_s;
typedef struct {
  ElLPApproach approach; 
  ElADMMCtrl_d admmCtrl;
  ElMehrotraCtrl_d mehrotraCtrl;
} ElLPDirectCtrl_d;

EL_EXPORT ElError ElLPDirectCtrlDefault_s
( ElLPDirectCtrl_s* ctrl, bool isSparse );
EL_EXPORT ElError ElLPDirectCtrlDefault_d
( ElLPDirectCtrl_d* ctrl, bool isSparse );

EL_EXPORT ElError ElLPDirectX_s
( ElConstMatrix_s A,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElLPDirectCtrl_s ctrl );
EL_EXPORT ElError ElLPDirectX_d
( ElConstMatrix_d A,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElLPDirectCtrl_d ctrl );

EL_EXPORT ElError ElLPDirectXDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b,
  ElConstDistMatrix_s c,
  ElDistMatrix_s x,
  ElDistMatrix_s y,
  ElDistMatrix_s z,
  ElLPDirectCtrl_s ctrl );
EL_EXPORT ElError ElLPDirectXDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b,
  ElConstDistMatrix_d c,
  ElDistMatrix_d x,
  ElDistMatrix_d y,
  ElDistMatrix_d z,
  ElLPDirectCtrl_d ctrl );

EL_EXPORT ElError ElLPDirectXSparse_s
( ElConstSparseMatrix_s A,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElLPDirectCtrl_s ctrl );
EL_EXPORT ElError ElLPDirectXSparse_d
( ElConstSparseMatrix_d A,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElLPDirectCtrl_d ctrl );

EL_EXPORT ElError ElLPDirectXDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b,
  ElConstDistMultiVec_s c,
  ElDistMultiVec_s x,
  ElDistMultiVec_s y,
  ElDistMultiVec_s z,
  ElLPDirectCtrl_s ctrl );
EL_EXPORT ElError ElLPDirectXDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b,
  ElConstDistMultiVec_d c,
  ElDistMultiVec_d x,
  ElDistMultiVec_d y,
  ElDistMultiVec_d z,
  ElLPDirectCtrl_d ctrl );

/* Affine conic form
   ----------------- */
EL_EXPORT ElError ElLPAffine_s
( ElConstMatrix_s A,
  ElConstMatrix_s G,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_s h,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElMatrix_s s );
EL_EXPORT ElError ElLPAffine_d
( ElConstMatrix_d A,
  ElConstMatrix_d G,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_d h,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElMatrix_d s );

EL_EXPORT ElError ElLPAffineDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s G,
  ElConstDistMatrix_s b,
  ElConstDistMatrix_s c,
  ElConstDistMatrix_s h,
  ElDistMatrix_s x,
  ElDistMatrix_s y,
  ElDistMatrix_s z,
  ElDistMatrix_s s );
EL_EXPORT ElError ElLPAffineDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d G,
  ElConstDistMatrix_d b,
  ElConstDistMatrix_d c,
  ElConstDistMatrix_d h,
  ElDistMatrix_d x,
  ElDistMatrix_d y,
  ElDistMatrix_d z,
  ElDistMatrix_d s );

EL_EXPORT ElError ElLPAffineSparse_s
( ElConstSparseMatrix_s A,
  ElConstSparseMatrix_s G,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_s h,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElMatrix_s s );
EL_EXPORT ElError ElLPAffineSparse_d
( ElConstSparseMatrix_d A,
  ElConstSparseMatrix_d G,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_d h,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElMatrix_d s );

EL_EXPORT ElError ElLPAffineDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistSparseMatrix_s G,
  ElConstDistMultiVec_s b,
  ElConstDistMultiVec_s c,
  ElConstDistMultiVec_s h,
  ElDistMultiVec_s x,
  ElDistMultiVec_s y,
  ElDistMultiVec_s z,
  ElDistMultiVec_s s );
EL_EXPORT ElError ElLPAffineDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistSparseMatrix_d G,
  ElConstDistMultiVec_d b,
  ElConstDistMultiVec_d c,
  ElConstDistMultiVec_d h,
  ElDistMultiVec_d x,
  ElDistMultiVec_d y,
  ElDistMultiVec_d z,
  ElDistMultiVec_d s );

/* Expert versions
   ^^^^^^^^^^^^^^^ */
typedef struct {
  ElLPApproach approach; 
  ElMehrotraCtrl_s mehrotraCtrl;
} ElLPAffineCtrl_s;
typedef struct {
  ElLPApproach approach; 
  ElMehrotraCtrl_d mehrotraCtrl;
} ElLPAffineCtrl_d;

EL_EXPORT ElError ElLPAffineCtrlDefault_s( ElLPAffineCtrl_s* ctrl );
EL_EXPORT ElError ElLPAffineCtrlDefault_d( ElLPAffineCtrl_d* ctrl );

EL_EXPORT ElError ElLPAffineX_s
( ElConstMatrix_s A,
  ElConstMatrix_s G,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_s h,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElMatrix_s s,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElLPAffineX_d
( ElConstMatrix_d A,
  ElConstMatrix_d G,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_d h,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElMatrix_d s,
  ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElLPAffineXDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s G,
  ElConstDistMatrix_s b,
  ElConstDistMatrix_s c,
  ElConstDistMatrix_s h,
  ElDistMatrix_s x,
  ElDistMatrix_s y,
  ElDistMatrix_s z,
  ElDistMatrix_s s,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElLPAffineXDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d G,
  ElConstDistMatrix_d b,
  ElConstDistMatrix_d c,
  ElConstDistMatrix_d h,
  ElDistMatrix_d x,
  ElDistMatrix_d y,
  ElDistMatrix_d z,
  ElDistMatrix_d s,
  ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElLPAffineXSparse_s
( ElConstSparseMatrix_s A,
  ElConstSparseMatrix_s G,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_s h,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElMatrix_s s,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElLPAffineXSparse_d
( ElConstSparseMatrix_d A,
  ElConstSparseMatrix_d G,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_d h,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElMatrix_d s,
  ElLPAffineCtrl_d ctrl );

EL_EXPORT ElError ElLPAffineXDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistSparseMatrix_s G,
  ElConstDistMultiVec_s b,
  ElConstDistMultiVec_s c,
  ElConstDistMultiVec_s h,
  ElDistMultiVec_s x,
  ElDistMultiVec_s y,
  ElDistMultiVec_s z,
  ElDistMultiVec_s s,
  ElLPAffineCtrl_s ctrl );
EL_EXPORT ElError ElLPAffineXDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistSparseMatrix_d G,
  ElConstDistMultiVec_d b,
  ElConstDistMultiVec_d c,
  ElConstDistMultiVec_d h,
  ElDistMultiVec_d x,
  ElDistMultiVec_d y,
  ElDistMultiVec_d z,
  ElDistMultiVec_d s,
  ElLPAffineCtrl_d ctrl );

/* Quadratic program
   ================= */
typedef enum {
  EL_QP_ADMM,
  EL_QP_MEHROTRA
} ElQPApproach;

/* Direct conic form
   ----------------- */
EL_EXPORT ElError ElQPDirect_s
( ElConstMatrix_s Q,
  ElConstMatrix_s A,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z );
EL_EXPORT ElError ElQPDirect_d
( ElConstMatrix_d Q,
  ElConstMatrix_d A,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z );

EL_EXPORT ElError ElQPDirectDist_s
( ElConstDistMatrix_s Q,
  ElConstDistMatrix_s A,
  ElConstDistMatrix_s b,
  ElConstDistMatrix_s c,
  ElDistMatrix_s x,
  ElDistMatrix_s y,
  ElDistMatrix_s z );
EL_EXPORT ElError ElQPDirectDist_d
( ElConstDistMatrix_d Q,
  ElConstDistMatrix_d A,
  ElConstDistMatrix_d b,
  ElConstDistMatrix_d c,
  ElDistMatrix_d x,
  ElDistMatrix_d y,
  ElDistMatrix_d z );

EL_EXPORT ElError ElQPDirectSparse_s
( ElConstSparseMatrix_s Q,
  ElConstSparseMatrix_s A,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z );
EL_EXPORT ElError ElQPDirectSparse_d
( ElConstSparseMatrix_d Q,
  ElConstSparseMatrix_d A,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z );

EL_EXPORT ElError ElQPDirectDistSparse_s
( ElConstDistSparseMatrix_s Q,
  ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b,
  ElConstDistMultiVec_s c,
  ElDistMultiVec_s x,
  ElDistMultiVec_s y,
  ElDistMultiVec_s z );
EL_EXPORT ElError ElQPDirectDistSparse_d
( ElConstDistSparseMatrix_d Q,
  ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b,
  ElConstDistMultiVec_d c,
  ElDistMultiVec_d x,
  ElDistMultiVec_d y,
  ElDistMultiVec_d z );

/* Expert versions
   ^^^^^^^^^^^^^^^ */
typedef struct {
  ElQPApproach approach; 
  ElMehrotraCtrl_s mehrotraCtrl;
} ElQPDirectCtrl_s;
typedef struct {
  ElQPApproach approach; 
  ElMehrotraCtrl_d mehrotraCtrl;
} ElQPDirectCtrl_d;

EL_EXPORT ElError ElQPDirectCtrlDefault_s( ElQPDirectCtrl_s* ctrl );
EL_EXPORT ElError ElQPDirectCtrlDefault_d( ElQPDirectCtrl_d* ctrl );

EL_EXPORT ElError ElQPDirectX_s
( ElConstMatrix_s Q,
  ElConstMatrix_s A,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElQPDirectCtrl_s ctrl );
EL_EXPORT ElError ElQPDirectX_d
( ElConstMatrix_d Q,
  ElConstMatrix_d A,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElQPDirectCtrl_d ctrl );

EL_EXPORT ElError ElQPDirectXDist_s
( ElConstDistMatrix_s Q,
  ElConstDistMatrix_s A,
  ElConstDistMatrix_s b,
  ElConstDistMatrix_s c,
  ElDistMatrix_s x,
  ElDistMatrix_s y,
  ElDistMatrix_s z,
  ElQPDirectCtrl_s ctrl );
EL_EXPORT ElError ElQPDirectXDist_d
( ElConstDistMatrix_d Q,
  ElConstDistMatrix_d A,
  ElConstDistMatrix_d b,
  ElConstDistMatrix_d c,
  ElDistMatrix_d x,
  ElDistMatrix_d y,
  ElDistMatrix_d z,
  ElQPDirectCtrl_d ctrl );

EL_EXPORT ElError ElQPDirectXSparse_s
( ElConstSparseMatrix_s Q,
  ElConstSparseMatrix_s A,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElQPDirectCtrl_s ctrl );
EL_EXPORT ElError ElQPDirectXSparse_d
( ElConstSparseMatrix_d Q,
  ElConstSparseMatrix_d A,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElQPDirectCtrl_d ctrl );

EL_EXPORT ElError ElQPDirectXDistSparse_s
( ElConstDistSparseMatrix_s Q,
  ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b,
  ElConstDistMultiVec_s c,
  ElDistMultiVec_s x,
  ElDistMultiVec_s y,
  ElDistMultiVec_s z,
  ElQPDirectCtrl_s ctrl );
EL_EXPORT ElError ElQPDirectXDistSparse_d
( ElConstDistSparseMatrix_d Q,
  ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b,
  ElConstDistMultiVec_d c,
  ElDistMultiVec_d x,
  ElDistMultiVec_d y,
  ElDistMultiVec_d z,
  ElQPDirectCtrl_d ctrl );

/* Affine conic form
   ----------------- */
EL_EXPORT ElError ElQPAffine_s
( ElConstMatrix_s Q,
  ElConstMatrix_s A,
  ElConstMatrix_s G,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_s h,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElMatrix_s s );
EL_EXPORT ElError ElQPAffine_d
( ElConstMatrix_d Q,
  ElConstMatrix_d A,
  ElConstMatrix_d G,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_d h,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElMatrix_d s );

EL_EXPORT ElError ElQPAffineDist_s
( ElConstDistMatrix_s Q,
  ElConstDistMatrix_s A,
  ElConstDistMatrix_s G,
  ElConstDistMatrix_s b,
  ElConstDistMatrix_s c,
  ElConstDistMatrix_s h,
  ElDistMatrix_s x,
  ElDistMatrix_s y,
  ElDistMatrix_s z,
  ElDistMatrix_s s );
EL_EXPORT ElError ElQPAffineDist_d
( ElConstDistMatrix_d Q,
  ElConstDistMatrix_d A,
  ElConstDistMatrix_d G,
  ElConstDistMatrix_d b,
  ElConstDistMatrix_d c,
  ElConstDistMatrix_d h,
  ElDistMatrix_d x,
  ElDistMatrix_d y,
  ElDistMatrix_d z,
  ElDistMatrix_d s );

EL_EXPORT ElError ElQPAffineSparse_s
( ElConstSparseMatrix_s Q,
  ElConstSparseMatrix_s A,
  ElConstSparseMatrix_s G,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_s h,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElMatrix_s s );
EL_EXPORT ElError ElQPAffineSparse_d
( ElConstSparseMatrix_d Q,
  ElConstSparseMatrix_d A,
  ElConstSparseMatrix_d G,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_d h,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElMatrix_d s );

EL_EXPORT ElError ElQPAffineDistSparse_s
( ElConstDistSparseMatrix_s Q,
  ElConstDistSparseMatrix_s A,
  ElConstDistSparseMatrix_s G,
  ElConstDistMultiVec_s b,
  ElConstDistMultiVec_s c,
  ElConstDistMultiVec_s h,
  ElDistMultiVec_s x,
  ElDistMultiVec_s y,
  ElDistMultiVec_s z,
  ElDistMultiVec_s s );
EL_EXPORT ElError ElQPAffineDistSparse_d
( ElConstDistSparseMatrix_d Q,
  ElConstDistSparseMatrix_d A,
  ElConstDistSparseMatrix_d G,
  ElConstDistMultiVec_d b,
  ElConstDistMultiVec_d c,
  ElConstDistMultiVec_d h,
  ElDistMultiVec_d x,
  ElDistMultiVec_d y,
  ElDistMultiVec_d z,
  ElDistMultiVec_d s );

/* Expert versions
   ^^^^^^^^^^^^^^^ */
typedef struct {
  ElQPApproach approach; 
  ElMehrotraCtrl_s mehrotraCtrl;
} ElQPAffineCtrl_s;
typedef struct {
  ElQPApproach approach; 
  ElMehrotraCtrl_d mehrotraCtrl;
} ElQPAffineCtrl_d;

EL_EXPORT ElError ElQPAffineCtrlDefault_s( ElQPAffineCtrl_s* ctrl );
EL_EXPORT ElError ElQPAffineCtrlDefault_d( ElQPAffineCtrl_d* ctrl );

EL_EXPORT ElError ElQPAffineX_s
( ElConstMatrix_s Q,
  ElConstMatrix_s A,
  ElConstMatrix_s G,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_s h,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElMatrix_s s,
  ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElQPAffineX_d
( ElConstMatrix_d Q,
  ElConstMatrix_d A,
  ElConstMatrix_d G,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_d h,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElMatrix_d s,
  ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElQPAffineXDist_s
( ElConstDistMatrix_s Q,
  ElConstDistMatrix_s A,
  ElConstDistMatrix_s G,
  ElConstDistMatrix_s b,
  ElConstDistMatrix_s c,
  ElConstDistMatrix_s h,
  ElDistMatrix_s x,
  ElDistMatrix_s y,
  ElDistMatrix_s z,
  ElDistMatrix_s s,
  ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElQPAffineXDist_d
( ElConstDistMatrix_d Q,
  ElConstDistMatrix_d A,
  ElConstDistMatrix_d G,
  ElConstDistMatrix_d b,
  ElConstDistMatrix_d c,
  ElConstDistMatrix_d h,
  ElDistMatrix_d x,
  ElDistMatrix_d y,
  ElDistMatrix_d z,
  ElDistMatrix_d s,
  ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElQPAffineXSparse_s
( ElConstSparseMatrix_s Q,
  ElConstSparseMatrix_s A,
  ElConstSparseMatrix_s G,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_s h,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElMatrix_s s,
  ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElQPAffineXSparse_d
( ElConstSparseMatrix_d Q,
  ElConstSparseMatrix_d A,
  ElConstSparseMatrix_d G,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_d h,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElMatrix_d s,
  ElQPAffineCtrl_d ctrl );

EL_EXPORT ElError ElQPAffineXDistSparse_s
( ElConstDistSparseMatrix_s Q,
  ElConstDistSparseMatrix_s A,
  ElConstDistSparseMatrix_s G,
  ElConstDistMultiVec_s b,
  ElConstDistMultiVec_s c,
  ElConstDistMultiVec_s h,
  ElDistMultiVec_s x,
  ElDistMultiVec_s y,
  ElDistMultiVec_s z,
  ElDistMultiVec_s s,
  ElQPAffineCtrl_s ctrl );
EL_EXPORT ElError ElQPAffineXDistSparse_d
( ElConstDistSparseMatrix_d Q,
  ElConstDistSparseMatrix_d A,
  ElConstDistSparseMatrix_d G,
  ElConstDistMultiVec_d b,
  ElConstDistMultiVec_d c,
  ElConstDistMultiVec_d h,
  ElDistMultiVec_d x,
  ElDistMultiVec_d y,
  ElDistMultiVec_d z,
  ElDistMultiVec_d s,
  ElQPAffineCtrl_d ctrl );

/* Box form
   -------- */

/* Alternating Direction Method of Multipliers
   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
EL_EXPORT ElError ElQPBoxADMM_s
( ElConstMatrix_s Q,
  ElConstMatrix_s C,
  float lb, float ub,
  ElMatrix_s X,
  ElInt* numIts );
EL_EXPORT ElError ElQPBoxADMM_d
( ElConstMatrix_d Q,
  ElConstMatrix_d C,
  double lb, double ub,
  ElMatrix_d X,
  ElInt* numIts );

EL_EXPORT ElError ElQPBoxADMMDist_s
( ElConstDistMatrix_s Q,
  ElConstDistMatrix_s C,
  float lb, float ub,
  ElDistMatrix_s X,
  ElInt* numIts );
EL_EXPORT ElError ElQPBoxADMMDist_d
( ElConstDistMatrix_d Q,
  ElConstDistMatrix_d C,
  double lb, double ub,
  ElDistMatrix_d X,
  ElInt* numIts );

/* Expert versions
   """"""""""""""" */
EL_EXPORT ElError ElQPBoxADMMX_s
( ElConstMatrix_s Q,
  ElConstMatrix_s C,
  float lb, float ub,
  ElMatrix_s X,
  ElADMMCtrl_s ctrl,
  ElInt* numIts );
EL_EXPORT ElError ElQPBoxADMMX_d
( ElConstMatrix_d Q,
  ElConstMatrix_d C,
  double lb, double ub,
  ElMatrix_d X,
  ElADMMCtrl_d ctrl,
  ElInt* numIts );

EL_EXPORT ElError ElQPBoxADMMXDist_s
( ElConstDistMatrix_s Q,
  ElConstDistMatrix_s C,
  float lb, float ub,
  ElDistMatrix_s X,
  ElADMMCtrl_s ctrl,
  ElInt* numIts );
EL_EXPORT ElError ElQPBoxADMMXDist_d
( ElConstDistMatrix_d Q,
  ElConstDistMatrix_d C,
  double lb, double ub,
  ElDistMatrix_d X,
  ElADMMCtrl_d ctrl,
  ElInt* numIts );

/* Second-order cone programs
   ========================== */
typedef enum {
  EL_SOCP_ADMM,
  EL_SOCP_MEHROTRA
} ElSOCPApproach;

/* Direct conic form
   ----------------- */
EL_EXPORT ElError ElSOCPDirect_s
( ElConstMatrix_s A,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z );
EL_EXPORT ElError ElSOCPDirect_d
( ElConstMatrix_d A,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z );

EL_EXPORT ElError ElSOCPDirectDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b,
  ElConstDistMatrix_s c,
  ElConstDistMatrix_i orders,
  ElConstDistMatrix_i firstInds,
  ElDistMatrix_s x,
  ElDistMatrix_s y,
  ElDistMatrix_s z );
EL_EXPORT ElError ElSOCPDirectDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b,
  ElConstDistMatrix_d c,
  ElConstDistMatrix_i orders,
  ElConstDistMatrix_i firstInds,
  ElDistMatrix_d x,
  ElDistMatrix_d y,
  ElDistMatrix_d z );

EL_EXPORT ElError ElSOCPDirectSparse_s
( ElConstSparseMatrix_s A,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z );
EL_EXPORT ElError ElSOCPDirectSparse_d
( ElConstSparseMatrix_d A,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z );

EL_EXPORT ElError ElSOCPDirectDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b,
  ElConstDistMultiVec_s c,
  ElConstDistMultiVec_i orders,
  ElConstDistMultiVec_i firstInds,
  ElDistMultiVec_s x,
  ElDistMultiVec_s y,
  ElDistMultiVec_s z );
EL_EXPORT ElError ElSOCPDirectDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b,
  ElConstDistMultiVec_d c,
  ElConstDistMultiVec_i orders,
  ElConstDistMultiVec_i firstInds,
  ElDistMultiVec_d x,
  ElDistMultiVec_d y,
  ElDistMultiVec_d z );

/* Expert versions
   ^^^^^^^^^^^^^^^ */
typedef struct {
  ElSOCPApproach approach; 
  ElMehrotraCtrl_s mehrotraCtrl;
} ElSOCPDirectCtrl_s;
typedef struct {
  ElSOCPApproach approach; 
  ElMehrotraCtrl_d mehrotraCtrl;
} ElSOCPDirectCtrl_d;

EL_EXPORT ElError ElSOCPDirectCtrlDefault_s( ElSOCPDirectCtrl_s* ctrl );
EL_EXPORT ElError ElSOCPDirectCtrlDefault_d( ElSOCPDirectCtrl_d* ctrl );

EL_EXPORT ElError ElSOCPDirectX_s
( ElConstMatrix_s A,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElSOCPDirectCtrl_s ctrl );
EL_EXPORT ElError ElSOCPDirectX_d
( ElConstMatrix_d A,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElSOCPDirectCtrl_d ctrl );

EL_EXPORT ElError ElSOCPDirectXDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s b,
  ElConstDistMatrix_s c,
  ElConstDistMatrix_i orders,
  ElConstDistMatrix_i firstInds,
  ElDistMatrix_s x,
  ElDistMatrix_s y,
  ElDistMatrix_s z,
  ElSOCPDirectCtrl_s ctrl );
EL_EXPORT ElError ElSOCPDirectXDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d b,
  ElConstDistMatrix_d c,
  ElConstDistMatrix_i orders,
  ElConstDistMatrix_i firstInds,
  ElDistMatrix_d x,
  ElDistMatrix_d y,
  ElDistMatrix_d z,
  ElSOCPDirectCtrl_d ctrl );

EL_EXPORT ElError ElSOCPDirectXSparse_s
( ElConstSparseMatrix_s A,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElSOCPDirectCtrl_s ctrl );
EL_EXPORT ElError ElSOCPDirectXSparse_d
( ElConstSparseMatrix_d A,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElSOCPDirectCtrl_d ctrl );

EL_EXPORT ElError ElSOCPDirectXDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistMultiVec_s b,
  ElConstDistMultiVec_s c,
  ElConstDistMultiVec_i orders,
  ElConstDistMultiVec_i firstInds,
  ElDistMultiVec_s x,
  ElDistMultiVec_s y,
  ElDistMultiVec_s z,
  ElSOCPDirectCtrl_s ctrl );
EL_EXPORT ElError ElSOCPDirectXDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistMultiVec_d b,
  ElConstDistMultiVec_d c,
  ElConstDistMultiVec_i orders,
  ElConstDistMultiVec_i firstInds,
  ElDistMultiVec_d x,
  ElDistMultiVec_d y,
  ElDistMultiVec_d z,
  ElSOCPDirectCtrl_d ctrl );

/* Affine conic form
   ----------------- */
EL_EXPORT ElError ElSOCPAffine_s
( ElConstMatrix_s A,
  ElConstMatrix_s G,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_s h,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElMatrix_s s );
EL_EXPORT ElError ElSOCPAffine_d
( ElConstMatrix_d A,
  ElConstMatrix_d G,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_d h,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElMatrix_d s );

EL_EXPORT ElError ElSOCPAffineDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s G,
  ElConstDistMatrix_s b,
  ElConstDistMatrix_s c,
  ElConstDistMatrix_s h,
  ElConstDistMatrix_i orders,
  ElConstDistMatrix_i firstInds,
  ElDistMatrix_s x,
  ElDistMatrix_s y,
  ElDistMatrix_s z,
  ElDistMatrix_s s );
EL_EXPORT ElError ElSOCPAffineDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d G,
  ElConstDistMatrix_d b,
  ElConstDistMatrix_d c,
  ElConstDistMatrix_d h,
  ElConstDistMatrix_i orders,
  ElConstDistMatrix_i firstInds,
  ElDistMatrix_d x,
  ElDistMatrix_d y,
  ElDistMatrix_d z,
  ElDistMatrix_d s );

EL_EXPORT ElError ElSOCPAffineSparse_s
( ElConstSparseMatrix_s A,
  ElConstSparseMatrix_s G,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_s h,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElMatrix_s s );
EL_EXPORT ElError ElSOCPAffineSparse_d
( ElConstSparseMatrix_d A,
  ElConstSparseMatrix_d G,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_d h,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElMatrix_d s );

EL_EXPORT ElError ElSOCPAffineDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistSparseMatrix_s G,
  ElConstDistMultiVec_s b,
  ElConstDistMultiVec_s c,
  ElConstDistMultiVec_s h,
  ElConstDistMultiVec_i orders,
  ElConstDistMultiVec_i firstInds,
  ElDistMultiVec_s x,
  ElDistMultiVec_s y,
  ElDistMultiVec_s z,
  ElDistMultiVec_s s );
EL_EXPORT ElError ElSOCPAffineDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistSparseMatrix_d G,
  ElConstDistMultiVec_d b,
  ElConstDistMultiVec_d c,
  ElConstDistMultiVec_d h,
  ElConstDistMultiVec_i orders,
  ElConstDistMultiVec_i firstInds,
  ElDistMultiVec_d x,
  ElDistMultiVec_d y,
  ElDistMultiVec_d z,
  ElDistMultiVec_d s );

/* Expert versions
   ^^^^^^^^^^^^^^^ */
typedef struct {
  ElSOCPApproach approach; 
  ElMehrotraCtrl_s mehrotraCtrl;
} ElSOCPAffineCtrl_s;
typedef struct {
  ElSOCPApproach approach; 
  ElMehrotraCtrl_d mehrotraCtrl;
} ElSOCPAffineCtrl_d;

EL_EXPORT ElError ElSOCPAffineCtrlDefault_s( ElSOCPAffineCtrl_s* ctrl );
EL_EXPORT ElError ElSOCPAffineCtrlDefault_d( ElSOCPAffineCtrl_d* ctrl );

EL_EXPORT ElError ElSOCPAffineX_s
( ElConstMatrix_s A,
  ElConstMatrix_s G,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_s h,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElMatrix_s s,
  ElSOCPAffineCtrl_s ctrl );
EL_EXPORT ElError ElSOCPAffineX_d
( ElConstMatrix_d A,
  ElConstMatrix_d G,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_d h,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElMatrix_d s,
  ElSOCPAffineCtrl_d ctrl );

EL_EXPORT ElError ElSOCPAffineXDist_s
( ElConstDistMatrix_s A,
  ElConstDistMatrix_s G,
  ElConstDistMatrix_s b,
  ElConstDistMatrix_s c,
  ElConstDistMatrix_s h,
  ElConstDistMatrix_i orders,
  ElConstDistMatrix_i firstInds,
  ElDistMatrix_s x,
  ElDistMatrix_s y,
  ElDistMatrix_s z,
  ElDistMatrix_s s,
  ElSOCPAffineCtrl_s ctrl );
EL_EXPORT ElError ElSOCPAffineXDist_d
( ElConstDistMatrix_d A,
  ElConstDistMatrix_d G,
  ElConstDistMatrix_d b,
  ElConstDistMatrix_d c,
  ElConstDistMatrix_d h,
  ElConstDistMatrix_i orders,
  ElConstDistMatrix_i firstInds,
  ElDistMatrix_d x,
  ElDistMatrix_d y,
  ElDistMatrix_d z,
  ElDistMatrix_d s,
  ElSOCPAffineCtrl_d ctrl );

EL_EXPORT ElError ElSOCPAffineXSparse_s
( ElConstSparseMatrix_s A,
  ElConstSparseMatrix_s G,
  ElConstMatrix_s b,
  ElConstMatrix_s c,
  ElConstMatrix_s h,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_s x,
  ElMatrix_s y,
  ElMatrix_s z,
  ElMatrix_s s,
  ElSOCPAffineCtrl_s ctrl );
EL_EXPORT ElError ElSOCPAffineXSparse_d
( ElConstSparseMatrix_d A,
  ElConstSparseMatrix_d G,
  ElConstMatrix_d b,
  ElConstMatrix_d c,
  ElConstMatrix_d h,
  ElConstMatrix_i orders,
  ElConstMatrix_i firstInds,
  ElMatrix_d x,
  ElMatrix_d y,
  ElMatrix_d z,
  ElMatrix_d s,
  ElSOCPAffineCtrl_d ctrl );

EL_EXPORT ElError ElSOCPAffineXDistSparse_s
( ElConstDistSparseMatrix_s A,
  ElConstDistSparseMatrix_s G,
  ElConstDistMultiVec_s b,
  ElConstDistMultiVec_s c,
  ElConstDistMultiVec_s h,
  ElConstDistMultiVec_i orders,
  ElConstDistMultiVec_i firstInds,
  ElDistMultiVec_s x,
  ElDistMultiVec_s y,
  ElDistMultiVec_s z,
  ElDistMultiVec_s s,
  ElSOCPAffineCtrl_s ctrl );
EL_EXPORT ElError ElSOCPAffineXDistSparse_d
( ElConstDistSparseMatrix_d A,
  ElConstDistSparseMatrix_d G,
  ElConstDistMultiVec_d b,
  ElConstDistMultiVec_d c,
  ElConstDistMultiVec_d h,
  ElConstDistMultiVec_i orders,
  ElConstDistMultiVec_i firstInds,
  ElDistMultiVec_d x,
  ElDistMultiVec_d y,
  ElDistMultiVec_d z,
  ElDistMultiVec_d s,
  ElSOCPAffineCtrl_d ctrl );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_OPTIMIZATION_SOLVERS_C_H */
