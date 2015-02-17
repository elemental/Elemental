/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El.h"
using namespace El;

extern "C" {

/* Linear programs
   =============== */
ElError ElLPIPFLineSearchCtrlDefault_s( ElLPIPFLineSearchCtrl_s* ctrl )
{
    ctrl->gamma = 1e-3;
    ctrl->beta = 2;
    ctrl->psi = 100;
    ctrl->stepRatio = 1.5;
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElLPIPFLineSearchCtrlDefault_d( ElLPIPFLineSearchCtrl_d* ctrl )
{
    ctrl->gamma = 1e-3;
    ctrl->beta = 2;
    ctrl->psi = 100;
    ctrl->stepRatio = 1.5;
    ctrl->print = false;
    return EL_SUCCESS;
}

/* Direct conic form
   ----------------- */
ElError ElLPDirectADMMCtrlDefault_s( ElLPDirectADMMCtrl_s* ctrl )
{
    ctrl->rho = 1;
    ctrl->alpha = 1.2;
    ctrl->maxIter = 500;
    ctrl->absTol = 1e-6;
    ctrl->relTol = 1e-4;
    ctrl->inv = true;
    ctrl->print = true;
    return EL_SUCCESS;
}

ElError ElLPDirectADMMCtrlDefault_d( ElLPDirectADMMCtrl_d* ctrl )
{
    ctrl->rho = 1;
    ctrl->alpha = 1.2;
    ctrl->maxIter = 500;
    ctrl->absTol = 1e-6;
    ctrl->relTol = 1e-4;
    ctrl->inv = true;
    ctrl->print = true;
    return EL_SUCCESS;
}

ElError ElLPDirectIPFCtrlDefault_s( ElLPDirectIPFCtrl_s* ctrl, bool isSparse )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 1000;
    ctrl->centering = 0.9;
    ctrl->system = ( isSparse ? EL_LP_PRIMAL_AUGMENTED_KKT
                              : EL_LP_PRIMAL_NORMAL_KKT );
    ElLPIPFLineSearchCtrlDefault_s( &ctrl->lineSearchCtrl );
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElLPDirectIPFCtrlDefault_d( ElLPDirectIPFCtrl_d* ctrl, bool isSparse )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 1000;
    ctrl->centering = 0.9;
    ctrl->system = ( isSparse ? EL_LP_PRIMAL_AUGMENTED_KKT
                              : EL_LP_PRIMAL_NORMAL_KKT );
    ElLPIPFLineSearchCtrlDefault_d( &ctrl->lineSearchCtrl );
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElLPDirectMehrotraCtrlDefault_s
( ElLPDirectMehrotraCtrl_s* ctrl, bool isSparse )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 100;
    ctrl->maxStepRatio = 0.99;
    ctrl->print = false;
    ctrl->system = ( isSparse ? EL_LP_PRIMAL_AUGMENTED_KKT 
                              : EL_LP_PRIMAL_NORMAL_KKT );;
    return EL_SUCCESS;
}

ElError ElLPDirectMehrotraCtrlDefault_d
( ElLPDirectMehrotraCtrl_d* ctrl, bool isSparse )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 100;
    ctrl->maxStepRatio = 0.99;
    ctrl->print = false;
    ctrl->system = ( isSparse ? EL_LP_PRIMAL_AUGMENTED_KKT 
                              : EL_LP_PRIMAL_NORMAL_KKT );;
    return EL_SUCCESS;
}

ElError ElLPDirectCtrlDefault_s( ElLPDirectCtrl_s* ctrl, bool isSparse )
{
    ctrl->approach = EL_LP_MEHROTRA;
    ElLPDirectADMMCtrlDefault_s( &ctrl->admmCtrl );
    ElLPDirectIPFCtrlDefault_s( &ctrl->ipfCtrl, isSparse );
    ElLPDirectMehrotraCtrlDefault_s( &ctrl->mehrotraCtrl, isSparse );
    return EL_SUCCESS;
}

ElError ElLPDirectCtrlDefault_d
( ElLPDirectCtrl_d* ctrl, bool isSparse )
{
    ctrl->approach = EL_LP_MEHROTRA;
    ElLPDirectADMMCtrlDefault_d( &ctrl->admmCtrl );
    ElLPDirectIPFCtrlDefault_d( &ctrl->ipfCtrl, isSparse );
    ElLPDirectMehrotraCtrlDefault_d( &ctrl->mehrotraCtrl, isSparse );
    return EL_SUCCESS;
}

/* Affine conic form
   ----------------- */
ElError ElLPAffineIPFCtrlDefault_s( ElLPAffineIPFCtrl_s* ctrl )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 1000;
    ctrl->centering = 0.9;
    ElLPIPFLineSearchCtrlDefault_s( &ctrl->lineSearchCtrl );
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElLPAffineIPFCtrlDefault_d( ElLPAffineIPFCtrl_d* ctrl )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 1000;
    ctrl->centering = 0.9;
    ElLPIPFLineSearchCtrlDefault_d( &ctrl->lineSearchCtrl );
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElLPAffineMehrotraCtrlDefault_s( ElLPAffineMehrotraCtrl_s* ctrl )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 100;
    ctrl->maxStepRatio = 0.99;
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElLPAffineMehrotraCtrlDefault_d( ElLPAffineMehrotraCtrl_d* ctrl )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 100;
    ctrl->maxStepRatio = 0.99;
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElLPAffineCtrlDefault_s( ElLPAffineCtrl_s* ctrl )
{
    ctrl->approach = EL_LP_MEHROTRA;
    ElLPAffineIPFCtrlDefault_s( &ctrl->ipfCtrl );
    ElLPAffineMehrotraCtrlDefault_s( &ctrl->mehrotraCtrl );
    return EL_SUCCESS;
}

ElError ElLPAffineCtrlDefault_d( ElLPAffineCtrl_d* ctrl )
{
    ctrl->approach = EL_LP_MEHROTRA;
    ElLPAffineIPFCtrlDefault_d( &ctrl->ipfCtrl );
    ElLPAffineMehrotraCtrlDefault_d( &ctrl->mehrotraCtrl );
    return EL_SUCCESS;
}

/* Quadratic programs
   ================== */
ElError ElQPIPFLineSearchCtrlDefault_s( ElQPIPFLineSearchCtrl_s* ctrl )
{
    ctrl->gamma = 1e-3;
    ctrl->beta = 2;
    ctrl->psi = 100;
    ctrl->stepRatio = 1.5;
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElQPIPFLineSearchCtrlDefault_d( ElQPIPFLineSearchCtrl_d* ctrl )
{
    ctrl->gamma = 1e-3;
    ctrl->beta = 2;
    ctrl->psi = 100;
    ctrl->stepRatio = 1.5;
    ctrl->print = false;
    return EL_SUCCESS;
}

/* Direct conic form
   ----------------- */
ElError ElQPDirectIPFCtrlDefault_s( ElQPDirectIPFCtrl_s* ctrl )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 1000;
    ctrl->centering = 0.9;
    ctrl->system = EL_QP_PRIMAL_AUGMENTED_KKT;
    ElQPIPFLineSearchCtrlDefault_s( &ctrl->lineSearchCtrl );
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElQPDirectIPFCtrlDefault_d( ElQPDirectIPFCtrl_d* ctrl )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 1000;
    ctrl->centering = 0.9;
    ctrl->system = EL_QP_PRIMAL_AUGMENTED_KKT;
    ElQPIPFLineSearchCtrlDefault_d( &ctrl->lineSearchCtrl );
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElQPDirectMehrotraCtrlDefault_s( ElQPDirectMehrotraCtrl_s* ctrl )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 100;
    ctrl->maxStepRatio = 0.99;
    ctrl->print = false;
    ctrl->system = EL_QP_PRIMAL_AUGMENTED_KKT;
    return EL_SUCCESS;
}

ElError ElQPDirectMehrotraCtrlDefault_d( ElQPDirectMehrotraCtrl_d* ctrl )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 100;
    ctrl->maxStepRatio = 0.99;
    ctrl->print = false;
    ctrl->system = EL_QP_PRIMAL_AUGMENTED_KKT;
    return EL_SUCCESS;
}

ElError ElQPDirectCtrlDefault_s( ElQPDirectCtrl_s* ctrl )
{
    ctrl->approach = EL_QP_MEHROTRA;
    ElQPDirectIPFCtrlDefault_s( &ctrl->ipfCtrl );
    ElQPDirectMehrotraCtrlDefault_s( &ctrl->mehrotraCtrl );
    return EL_SUCCESS;
}

ElError ElQPDirectCtrlDefault_d( ElQPDirectCtrl_d* ctrl )
{
    ctrl->approach = EL_QP_MEHROTRA;
    ElQPDirectIPFCtrlDefault_d( &ctrl->ipfCtrl );
    ElQPDirectMehrotraCtrlDefault_d( &ctrl->mehrotraCtrl );
    return EL_SUCCESS;
}

/* Affine conic form
   ----------------- */
ElError ElQPAffineIPFCtrlDefault_s( ElQPAffineIPFCtrl_s* ctrl )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 1000;
    ctrl->centering = 0.9;
    ElQPIPFLineSearchCtrlDefault_s( &ctrl->lineSearchCtrl );
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElQPAffineIPFCtrlDefault_d( ElQPAffineIPFCtrl_d* ctrl )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 1000;
    ctrl->centering = 0.9;
    ElQPIPFLineSearchCtrlDefault_d( &ctrl->lineSearchCtrl );
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElQPAffineMehrotraCtrlDefault_s( ElQPAffineMehrotraCtrl_s* ctrl )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 100;
    ctrl->maxStepRatio = 0.99;
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElQPAffineMehrotraCtrlDefault_d( ElQPAffineMehrotraCtrl_d* ctrl )
{
    ctrl->primalInitialized = false;
    ctrl->dualInitialized = false;
    ctrl->tol = 1e-8;
    ctrl->maxIts = 100;
    ctrl->maxStepRatio = 0.99;
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElQPAffineCtrlDefault_s( ElQPAffineCtrl_s* ctrl )
{
    ctrl->approach = EL_QP_MEHROTRA;
    ElQPAffineIPFCtrlDefault_s( &ctrl->ipfCtrl );
    ElQPAffineMehrotraCtrlDefault_s( &ctrl->mehrotraCtrl );
    return EL_SUCCESS;
}

ElError ElQPAffineCtrlDefault_d( ElQPAffineCtrl_d* ctrl )
{
    ctrl->approach = EL_QP_MEHROTRA;
    ElQPAffineIPFCtrlDefault_d( &ctrl->ipfCtrl );
    ElQPAffineMehrotraCtrlDefault_d( &ctrl->mehrotraCtrl );
    return EL_SUCCESS;
}

#define C_PROTO_REAL(SIG,Real) \
  /* Linear program
     ============== */ \
  /* Direct conic form
     ----------------- */ \
  ElError ElLPDirect_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x,      ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z ) \
  { EL_TRY( LP( *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  ElError ElLPDirectDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElDistMatrix_ ## SIG x,      ElDistMatrix_ ## SIG y, \
    ElDistMatrix_ ## SIG z ) \
  { EL_TRY( LP( *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  ElError ElLPDirectSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b,       ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x,            ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z ) \
  { EL_TRY( LP( *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  ElError ElLPDirectDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b,     ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG x,          ElDistMultiVec_ ## SIG y, \
    ElDistMultiVec_ ## SIG z ) \
  { EL_TRY( LP( *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  /* Expert version
     ^^^^^^^^^^^^^^ */ \
  ElError ElLPDirectX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x,      ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z, \
    ElLPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( LP( *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), \
       CReflect(ctrl) ) ) } \
  ElError ElLPDirectXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElDistMatrix_ ## SIG x,      ElDistMatrix_ ## SIG y, \
    ElDistMatrix_ ## SIG z, \
    ElLPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( LP( *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), \
       CReflect(ctrl) ) ) } \
  ElError ElLPDirectXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b,       ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x,            ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z, \
    ElLPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( LP( *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), \
       CReflect(ctrl) ) ) } \
  ElError ElLPDirectXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b,     ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG x,          ElDistMultiVec_ ## SIG y, \
    ElDistMultiVec_ ## SIG z, \
    ElLPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( LP( *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), \
       CReflect(ctrl) ) ) } \
  /* Affine conic form
     ----------------- */ \
  ElError ElLPAffine_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG G, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElConstMatrix_ ## SIG h, \
    ElMatrix_ ## SIG x,      ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z,      ElMatrix_ ## SIG s ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s) ) ) } \
  ElError ElLPAffineDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG G, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElConstDistMatrix_ ## SIG h, \
    ElDistMatrix_ ## SIG x,      ElDistMatrix_ ## SIG y, \
    ElDistMatrix_ ## SIG z,      ElDistMatrix_ ## SIG s ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s) ) ) } \
  ElError ElLPAffineSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstSparseMatrix_ ## SIG G, \
    ElConstMatrix_ ## SIG b,       ElConstMatrix_ ## SIG c, \
    ElConstMatrix_ ## SIG h, \
    ElMatrix_ ## SIG x,            ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z,            ElMatrix_ ## SIG s ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s) ) ) } \
  ElError ElLPAffineDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistSparseMatrix_ ## SIG G, \
    ElConstDistMultiVec_ ## SIG b,     ElConstDistMultiVec_ ## SIG c, \
    ElConstDistMultiVec_ ## SIG h, \
    ElDistMultiVec_ ## SIG x,          ElDistMultiVec_ ## SIG y, \
    ElDistMultiVec_ ## SIG z,          ElDistMultiVec_ ## SIG s ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s) ) ) } \
  /* Expert version
     ^^^^^^^^^^^^^^ */ \
  ElError ElLPAffineX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG G, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElConstMatrix_ ## SIG h, \
    ElMatrix_ ## SIG x,      ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z,      ElMatrix_ ## SIG s, \
    ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s), \
       CReflect(ctrl) ) ) } \
  ElError ElLPAffineXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG G, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElConstDistMatrix_ ## SIG h, \
    ElDistMatrix_ ## SIG x,      ElDistMatrix_ ## SIG y, \
    ElDistMatrix_ ## SIG z,      ElDistMatrix_ ## SIG s, \
    ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s), \
       CReflect(ctrl) ) ) } \
  ElError ElLPAffineXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstSparseMatrix_ ## SIG G, \
    ElConstMatrix_ ## SIG b,       ElConstMatrix_ ## SIG c, \
    ElConstMatrix_ ## SIG h, \
    ElMatrix_ ## SIG x,            ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z,            ElMatrix_ ## SIG s, \
    ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s), \
       CReflect(ctrl) ) ) } \
  ElError ElLPAffineXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistSparseMatrix_ ## SIG G, \
    ElConstDistMultiVec_ ## SIG b,     ElConstDistMultiVec_ ## SIG c, \
    ElConstDistMultiVec_ ## SIG h, \
    ElDistMultiVec_ ## SIG x,          ElDistMultiVec_ ## SIG y, \
    ElDistMultiVec_ ## SIG z,          ElDistMultiVec_ ## SIG s, \
    ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s), \
       CReflect(ctrl) ) ) } \
  /* Self-dual conic form
     -------------------- */ \
  /* TODO */ \
  /* Quadratic program
     ================= */ \
  /* Direct conic form
     ----------------- */ \
  ElError ElQPDirect_ ## SIG \
  ( ElConstMatrix_ ## SIG Q, ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x,      ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  ElError ElQPDirectDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG Q, ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElDistMatrix_ ## SIG x,      ElDistMatrix_ ## SIG y, \
    ElDistMatrix_ ## SIG z ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  ElError ElQPDirectSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG Q, ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b,       ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x,            ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  ElError ElQPDirectDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG Q, ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b,     ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG x,          ElDistMultiVec_ ## SIG y, \
    ElDistMultiVec_ ## SIG z ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  /* Expert version
     ^^^^^^^^^^^^^^ */ \
  ElError ElQPDirectX_ ## SIG \
  ( ElConstMatrix_ ## SIG Q, ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x,      ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z, \
    ElQPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), \
       CReflect(ctrl) ) ) } \
  ElError ElQPDirectXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG Q, ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElDistMatrix_ ## SIG x,      ElDistMatrix_ ## SIG y, \
    ElDistMatrix_ ## SIG z, \
    ElQPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), \
       CReflect(ctrl) ) ) } \
  ElError ElQPDirectXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG Q, ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b,       ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x,            ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z, \
    ElQPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), \
       CReflect(ctrl) ) ) } \
  ElError ElQPDirectXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG Q, ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b,     ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG x,          ElDistMultiVec_ ## SIG y, \
    ElDistMultiVec_ ## SIG z, \
    ElQPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), \
       CReflect(ctrl) ) ) } \
  /* Affine conic form
     ----------------- */ \
  ElError ElQPAffine_ ## SIG \
  ( ElConstMatrix_ ## SIG Q, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG G, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElConstMatrix_ ## SIG h, \
    ElMatrix_ ## SIG x,      ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z,      ElMatrix_ ## SIG s ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s) ) ) } \
  ElError ElQPAffineDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG Q, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG G, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElConstDistMatrix_ ## SIG h, \
    ElDistMatrix_ ## SIG x,      ElDistMatrix_ ## SIG y, \
    ElDistMatrix_ ## SIG z,      ElDistMatrix_ ## SIG s ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s) ) ) } \
  ElError ElQPAffineSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG Q, \
    ElConstSparseMatrix_ ## SIG A, ElConstSparseMatrix_ ## SIG G, \
    ElConstMatrix_ ## SIG b,       ElConstMatrix_ ## SIG c, \
    ElConstMatrix_ ## SIG h, \
    ElMatrix_ ## SIG x,            ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z,            ElMatrix_ ## SIG s ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s) ) ) } \
  ElError ElQPAffineDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG Q, \
    ElConstDistSparseMatrix_ ## SIG A, ElConstDistSparseMatrix_ ## SIG G, \
    ElConstDistMultiVec_ ## SIG b,     ElConstDistMultiVec_ ## SIG c, \
    ElConstDistMultiVec_ ## SIG h, \
    ElDistMultiVec_ ## SIG x,          ElDistMultiVec_ ## SIG y, \
    ElDistMultiVec_ ## SIG z,          ElDistMultiVec_ ## SIG s ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s) ) ) } \
  /* Expert version
     ^^^^^^^^^^^^^^ */ \
  ElError ElQPAffineX_ ## SIG \
  ( ElConstMatrix_ ## SIG Q, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG G, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElConstMatrix_ ## SIG h, \
    ElMatrix_ ## SIG x,      ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z,      ElMatrix_ ## SIG s, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s), \
       CReflect(ctrl) ) ) } \
  ElError ElQPAffineXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG Q, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG G, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElConstDistMatrix_ ## SIG h, \
    ElDistMatrix_ ## SIG x,      ElDistMatrix_ ## SIG y, \
    ElDistMatrix_ ## SIG z,      ElDistMatrix_ ## SIG s, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s), \
       CReflect(ctrl) ) ) } \
  ElError ElQPAffineXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG Q, \
    ElConstSparseMatrix_ ## SIG A, ElConstSparseMatrix_ ## SIG G, \
    ElConstMatrix_ ## SIG b,       ElConstMatrix_ ## SIG c, \
    ElConstMatrix_ ## SIG h, \
    ElMatrix_ ## SIG x,            ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z,            ElMatrix_ ## SIG s, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s), \
       CReflect(ctrl) ) ) } \
  ElError ElQPAffineXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG Q, \
    ElConstDistSparseMatrix_ ## SIG A, ElConstDistSparseMatrix_ ## SIG G, \
    ElConstDistMultiVec_ ## SIG b,     ElConstDistMultiVec_ ## SIG c, \
    ElConstDistMultiVec_ ## SIG h, \
    ElDistMultiVec_ ## SIG x,          ElDistMultiVec_ ## SIG y, \
    ElDistMultiVec_ ## SIG z,          ElDistMultiVec_ ## SIG s, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( QP( *CReflect(Q), *CReflect(A), *CReflect(G), \
      *CReflect(b), *CReflect(c), *CReflect(h), \
      *CReflect(x), *CReflect(y), *CReflect(z), *CReflect(s), \
       CReflect(ctrl) ) ) } \
  /* Self-dual conic form
     -------------------- */ \
  /* TODO */ \
  /* Box form (no linear equalities)
     ------------------------------- */ \
  ElError ElQPBoxADMM_ ## SIG \
  ( ElConstMatrix_ ## SIG Q, ElConstMatrix_ ## SIG C, \
    Real lb, Real ub, ElMatrix_ ## SIG Z, ElInt* numIts ) \
  { EL_TRY( *numIts = qp::box::ADMM( *CReflect(Q), *CReflect(C), lb, ub, \
      *CReflect(Z) ) ) } \
  ElError ElQPBoxADMMDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG Q, ElConstDistMatrix_ ## SIG C, \
    Real lb, Real ub, ElDistMatrix_ ## SIG Z, ElInt* numIts ) \
  { EL_TRY( *numIts = qp::box::ADMM( *CReflect(Q), *CReflect(C), lb, ub, \
      *CReflect(Z) ) ) } \
  /* Affine conic form
     ----------------- */ \
  /* TODO */ \
  /* Self-dual conic form
     -------------------- */ \
  /* TODO */

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
