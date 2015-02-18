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

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* Basis pursuit
     ============= */ \
  /* ADMM
     ---- */ \
  ElError ElBPADMM_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG z, ElInt* numIts ) \
  { EL_TRY( *numIts = bp::ADMM \
      ( *CReflect(A), *CReflect(b), *CReflect(z) ) ) } \
  ElError ElBPADMMDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG z, ElInt* numIts ) \
  { EL_TRY( *numIts = bp::ADMM \
      ( *CReflect(A), *CReflect(b), *CReflect(z) ) ) } \
  /* BPDN
     ==== */ \
  /* ADMM
     ---- */ \
  ElError ElBPDNADMM_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, Base<F> lambda, \
    ElMatrix_ ## SIG z, ElInt* numIts ) \
  { EL_TRY( *numIts = bpdn::ADMM( *CReflect(A), *CReflect(b), lambda, \
      *CReflect(z) ) ) } \
  ElError ElBPDNADMMDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, Base<F> lambda, \
    ElDistMatrix_ ## SIG z, ElInt* numIts ) \
  { EL_TRY( *numIts = bpdn::ADMM( *CReflect(A), *CReflect(b), lambda, \
      *CReflect(z) ) ) } \
  /* Robust Principal Component Analysis
     =================================== */ \
  ElError ElRPCA_ ## SIG \
  ( ElConstMatrix_ ## SIG M, ElMatrix_ ## SIG L, ElMatrix_ ## SIG S ) \
  { EL_TRY( RPCA( *CReflect(M), *CReflect(L), *CReflect(S) ) ) } \
  ElError ElRPCADist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG M, \
    ElDistMatrix_ ## SIG L, ElDistMatrix_ ## SIG S ) \
  { EL_TRY( RPCA( *CReflect(M), *CReflect(L), *CReflect(S) ) ) } \
  /* Sparse inverse covariance selection
     =================================== */ \
  ElError ElSparseInvCov_ ## SIG \
  ( ElConstMatrix_ ## SIG D, Base<F> lambda, ElMatrix_ ## SIG Z, \
    ElInt* numIts ) \
  { EL_TRY( *numIts = SparseInvCov( *CReflect(D), lambda, *CReflect(Z) ) ) } \
  ElError ElSparseInvCovDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG D, Base<F> lambda, ElDistMatrix_ ## SIG Z, \
    ElInt* numIts ) \
  { EL_TRY( *numIts = SparseInvCov( *CReflect(D), lambda, *CReflect(Z) ) ) }

#define C_PROTO_REAL(SIG,Real) \
  C_PROTO_FIELD(SIG,SIG,Real) \
  /* Basis Pursuit
     ============= */ \
  ElError ElBP_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( BP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElBPDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG x ) \
  { EL_TRY( BP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElBPSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( BP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElBPDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( BP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElBPX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x, ElLPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( BP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElBPXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG x, ElLPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( BP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElBPXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x, ElLPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( BP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElBPXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElDistMultiVec_ ## SIG x, ElLPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( BP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  /* Chebyshev point
     =============== */ \
  ElError ElCP_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( CP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElCPDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG x ) \
  { EL_TRY( CP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElCPSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( CP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElCPDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( CP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElCPX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( CP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElCPXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( CP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElCPXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( CP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElCPXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElDistMultiVec_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( CP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  /* Dantzig Selector
     ================ */ \
  ElError ElDS_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x ) \
  { EL_TRY( DS( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElDSDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real lambda, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( DS( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElDSSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x ) \
  { EL_TRY( DS( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElDSDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real lambda, ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( DS( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElDSX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( DS \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElDSXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real lambda, ElDistMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( DS \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElDSXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( DS \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElDSXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real lambda, ElDistMultiVec_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( DS \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  /* Least Absolute Value regression
     =============================== */ \
  ElError ElLAV_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( LAV( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElLAVDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG x ) \
  { EL_TRY( LAV( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElLAVSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( LAV( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElLAVDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( LAV( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElLAVX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( LAV \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElLAVXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( LAV \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElLAVXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( LAV \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElLAVXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElDistMultiVec_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( LAV \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  /* Basis Pursuit Denoising
     ======================= */ \
  ElError ElBPDN_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x ) \
  { EL_TRY( BPDN( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElBPDNDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real lambda, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( BPDN( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElBPDNSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x ) \
  { EL_TRY( BPDN( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElBPDNDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real lambda, ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( BPDN( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElBPDNX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x, ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( BPDN \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElBPDNXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real lambda, ElDistMatrix_ ## SIG x, ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( BPDN \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElBPDNXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x, ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( BPDN \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElBPDNXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real lambda, ElDistMultiVec_ ## SIG x, ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( BPDN \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  /* Elastic net
     =========== */ \
  ElError ElEN_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda1, Real lambda2, ElMatrix_ ## SIG x ) \
  { EL_TRY( EN( *CReflect(A), *CReflect(b), \
      lambda1, lambda2, *CReflect(x) ) ) } \
  ElError ElENDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real lambda1, Real lambda2, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( EN( *CReflect(A), *CReflect(b), \
      lambda1, lambda2, *CReflect(x) ) ) } \
  ElError ElENSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda1, Real lambda2, ElMatrix_ ## SIG x ) \
  { EL_TRY( EN( *CReflect(A), *CReflect(b), \
      lambda1, lambda2, *CReflect(x) ) ) } \
  ElError ElENDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real lambda1, Real lambda2, ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( EN( *CReflect(A), *CReflect(b), \
      lambda1, lambda2, *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElENX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda1, Real lambda2, ElMatrix_ ## SIG x, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( EN \
      ( *CReflect(A), *CReflect(b), lambda1, lambda2, \
        *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElENXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real lambda1, Real lambda2, ElDistMatrix_ ## SIG x, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( EN \
      ( *CReflect(A), *CReflect(b), lambda1, lambda2, \
        *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElENXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda1, Real lambda2, ElMatrix_ ## SIG x, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( EN \
      ( *CReflect(A), *CReflect(b), lambda1, lambda2, \
        *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElENXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real lambda1, Real lambda2, ElDistMultiVec_ ## SIG x, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( EN \
      ( *CReflect(A), *CReflect(b), lambda1, lambda2, \
        *CReflect(x), CReflect(ctrl) ) ) } \
  /* Support Vector Machine (soft-margin) 
     ==================================== */ \
  ElError ElSVM_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG d, \
    Real lambda, ElMatrix_ ## SIG x ) \
  { EL_TRY( SVM( *CReflect(A), *CReflect(d), lambda, *CReflect(x) ) ) } \
  ElError ElSVMDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG d, \
    Real lambda, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( SVM( *CReflect(A), *CReflect(d), lambda, *CReflect(x) ) ) } \
  ElError ElSVMSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG d, \
    Real lambda, ElMatrix_ ## SIG x ) \
  { EL_TRY( SVM( *CReflect(A), *CReflect(d), lambda, *CReflect(x) ) ) } \
  ElError ElSVMDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG d, \
    Real lambda, ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( SVM( *CReflect(A), *CReflect(d), lambda, *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElSVMX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG d, \
    Real lambda, ElMatrix_ ## SIG x, ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( SVM \
      ( *CReflect(A), *CReflect(d), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElSVMXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG d, \
    Real lambda, ElDistMatrix_ ## SIG x, ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( SVM \
      ( *CReflect(A), *CReflect(d), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElSVMXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG d, \
    Real lambda, ElMatrix_ ## SIG x, ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( SVM \
      ( *CReflect(A), *CReflect(d), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElSVMXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG d, \
    Real lambda, ElDistMultiVec_ ## SIG x, ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( SVM \
      ( *CReflect(A), *CReflect(d), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  /* Total variation denoising 
     ========================= */ \
  ElError ElTV_ ## SIG \
  ( ElConstMatrix_ ## SIG b, Real lambda, ElMatrix_ ## SIG x ) \
  { EL_TRY( TV( *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElTVDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG b, Real lambda, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( TV( *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElTVDistSparse_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG b, Real lambda, ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( TV( *CReflect(b), lambda, *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElTVX_ ## SIG \
  ( ElConstMatrix_ ## SIG b, Real lambda, ElMatrix_ ## SIG x, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( TV \
      ( *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElTVXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG b, Real lambda, ElDistMatrix_ ## SIG x, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( TV \
      ( *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElTVXDistSparse_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG b, Real lambda, ElDistMultiVec_ ## SIG x, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( TV \
      ( *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  /* Logistic regression
     =================== */ \
  ElError ElLogisticRegression_ ## SIG \
  ( ElConstMatrix_ ## SIG G, ElConstMatrix_ ## SIG q, \
    ElMatrix_ ## SIG z, Real gamma, ElRegularization penalty, \
    ElInt* numIts ) \
  { EL_TRY( *numIts = LogisticRegression( *CReflect(G), *CReflect(q), \
      *CReflect(z), gamma, CReflect(penalty) ) ) } \
  ElError ElLogisticRegressionDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG G, ElConstDistMatrix_ ## SIG q, \
    ElDistMatrix_ ## SIG z, Real gamma, ElRegularization penalty, \
    ElInt* numIts ) \
  { EL_TRY( *numIts = LogisticRegression( *CReflect(G), *CReflect(q), \
      *CReflect(z), gamma, CReflect(penalty) ) ) } \
  /* ModelFit 
     ======== */ \
  ElError ElModelFit_ ## SIG \
  ( void (*lossProx)(ElMatrix_ ## SIG,Real), \
    void (*regProx)(ElMatrix_ ## SIG,Real), \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG w, Real rho, ElInt* numIts ) \
  { try { \
      auto lossLambda = \
        [&]( Matrix<Real>& B, Real tau ) { lossProx(CReflect(&B),tau); }; \
      auto regLambda = \
        [&]( Matrix<Real>& B, Real tau ) { regProx(CReflect(&B),tau); }; \
      *numIts = ModelFit \
        ( function<void(Matrix<Real>&,Real)>(lossLambda), \
          function<void(Matrix<Real>&,Real)>(regLambda), \
          *CReflect(A), *CReflect(b), *CReflect(w), rho ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElModelFitDist_ ## SIG \
  ( void (*lossProx)(ElDistMatrix_ ## SIG,Real), \
    void (*regProx)(ElDistMatrix_ ## SIG,Real), \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG w, Real rho, ElInt* numIts ) \
  { try { \
      auto lossLambda = \
        [&]( DistMatrix<Real>& B, Real tau ) { lossProx(CReflect(&B),tau); }; \
      auto regLambda = \
        [&]( DistMatrix<Real>& B, Real tau ) { regProx(CReflect(&B),tau); }; \
      *numIts = ModelFit \
        ( function<void(DistMatrix<Real>&,Real)>(lossLambda), \
          function<void(DistMatrix<Real>&,Real)>(regLambda), \
          *CReflect(A), *CReflect(b), *CReflect(w), rho ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Non-negative matrix factorization
     ================================= */ \
  ElError ElNMF_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG X, ElMatrix_ ## SIG Y ) \
  { EL_TRY( NMF( *CReflect(A), *CReflect(X), *CReflect(Y) ) ) } \
  ElError ElNMFDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG X, ElDistMatrix_ ## SIG Y ) \
  { EL_TRY( NMF( *CReflect(A), *CReflect(X), *CReflect(Y) ) ) } \
  /* Expert versions */ \
  ElError ElNMFX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG X, ElMatrix_ ## SIG Y, \
    ElQPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( NMF( *CReflect(A), *CReflect(X), *CReflect(Y), \
      CReflect(ctrl) ) ) } \
  ElError ElNMFXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG X, ElDistMatrix_ ## SIG Y, \
    ElQPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( NMF( *CReflect(A), *CReflect(X), *CReflect(Y), \
      CReflect(ctrl) ) ) } \
  /* Non-negative least squares
     ========================== */ \
  ElError ElNNLS_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElMatrix_ ## SIG X ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X) ) ) } \
  ElError ElNNLSDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG X ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X) ) ) } \
  ElError ElNNLSSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElMatrix_ ## SIG X ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X) ) ) } \
  ElError ElNNLSDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG B, \
    ElDistMultiVec_ ## SIG X ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X) ) ) } \
  /* Expert version */ \
  ElError ElNNLSX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElMatrix_ ## SIG X, ElQPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X), \
      CReflect(ctrl) ) ) } \
  ElError ElNNLSXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG X, ElQPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X), \
      CReflect(ctrl) ) ) } \
  ElError ElNNLSXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElMatrix_ ## SIG X, ElQPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X), \
      CReflect(ctrl) ) ) } \
  ElError ElNNLSXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG B, \
    ElDistMultiVec_ ## SIG X, ElQPDirectCtrl_ ## SIG ctrl ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X), \
      CReflect(ctrl) ) ) } \
  /* ADMM
     ---- */ \
  ElError ElNNLSADMM_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElMatrix_ ## SIG X, ElInt* numIts ) \
  { EL_TRY( *numIts = nnls::ADMM( *CReflect(A), *CReflect(B), \
      *CReflect(X) ) ) } \
  ElError ElNNLSADMMDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG X, ElInt* numIts ) \
  { EL_TRY( *numIts = nnls::ADMM( *CReflect(A), *CReflect(B), \
      *CReflect(X) ) ) } \
  /* Support Vector Machine
     ====================== */ \
  /* TODO */ \
  /* ADMM
     ---- */ \
  ElError ElSVMADMM_ ## SIG \
  ( ElConstMatrix_ ## SIG G, ElConstMatrix_ ## SIG q, \
    Real gamma, ElMatrix_ ## SIG z, ElInt* numIts ) \
  { EL_TRY( *numIts = \
      svm::ADMM( *CReflect(G), *CReflect(q), gamma, *CReflect(z) ) ) } \
  ElError ElSVMADMMDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG G, ElConstDistMatrix_ ## SIG q, \
    Real gamma, ElDistMatrix_ ## SIG z, ElInt* numIts ) \
  { EL_TRY( *numIts = \
      svm::ADMM( *CReflect(G), *CReflect(q), gamma, *CReflect(z) ) ) }

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
