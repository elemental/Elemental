/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/lapack_like/factor.hpp>
#include <El-lite.h>
#include <El/lapack_like/factor.h>
using namespace El;

extern "C" {

ElError ElRegSolveCtrlDefault_s( ElRegSolveCtrl_s* ctrl )
{
    const float eps = limits::Epsilon<float>();
    ctrl->alg = EL_REG_SOLVE_FGMRES;
    ctrl->relTol = Pow(eps,float(0.5));
    ctrl->relTolRefine = Pow(eps,float(0.8));
    ctrl->maxIts = 4;
    ctrl->maxRefineIts = 2;
    ctrl->restart = 4;
    ctrl->progress = false;
    ctrl->time = false;
    return EL_SUCCESS;
}

ElError ElRegSolveCtrlDefault_d( ElRegSolveCtrl_d* ctrl )
{
    const double eps = limits::Epsilon<double>();
    ctrl->alg = EL_REG_SOLVE_FGMRES;
    ctrl->relTol = Pow(eps,0.5);
    ctrl->relTolRefine = Pow(eps,0.8);
    ctrl->maxIts = 4;
    ctrl->maxRefineIts = 2;
    ctrl->restart = 4;
    ctrl->progress = false;
    ctrl->time = false;
    return EL_SUCCESS;
}

ElError ElQRCtrlDefault_s( ElQRCtrl_s* ctrl )
{
    ctrl->colPiv = false;
    ctrl->boundRank = false;
    ctrl->maxRank = 0;
    ctrl->adaptive = false;
    ctrl->tol = 0;
    ctrl->alwaysRecomputeNorms = false;
    ctrl->smallestFirst = false;
    return EL_SUCCESS;
}

ElError ElQRCtrlDefault_d( ElQRCtrl_d* ctrl )
{
    ctrl->colPiv = false;
    ctrl->boundRank = false;
    ctrl->maxRank = 0;
    ctrl->adaptive = false;
    ctrl->tol = 0;
    ctrl->alwaysRecomputeNorms = false;
    ctrl->smallestFirst = false;
    return EL_SUCCESS;
}

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* Cholesky
     ======== */ \
  /* Cholesky (no pivoting) */ \
  ElError ElCholesky_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( Cholesky( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElCholeskyDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Cholesky( CReflect(uplo), *CReflect(A) ) ) } \
  /* Reverse Cholesky (no pivoting) */ \
  ElError ElReverseCholesky_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( ReverseCholesky( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElReverseCholeskyDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( ReverseCholesky( CReflect(uplo), *CReflect(A) ) ) } \
  /* Cholesky (full pivoting) */ \
  ElError ElCholeskyPiv_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElPermutation p ) \
  { EL_TRY( Cholesky( CReflect(uplo), *CReflect(A), *CReflect(p) ) ) } \
  ElError ElCholeskyPivDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistPermutation p ) \
  { EL_TRY( Cholesky( CReflect(uplo), *CReflect(A), *CReflect(p) ) ) } \
  /* Cholesky low-rank modification */ \
  ElError ElCholeskyMod_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG T, \
    Base<F> alpha, ElMatrix_ ## SIG V ) \
  { EL_TRY( \
      CholeskyMod( CReflect(uplo), *CReflect(T), alpha, *CReflect(V) ) ) } \
  ElError ElCholeskyModDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG T, \
    Base<F> alpha, ElDistMatrix_ ## SIG V ) \
  { EL_TRY( \
      CholeskyMod( CReflect(uplo), *CReflect(T), alpha, *CReflect(V) ) ) } \
  /* Hermitian Positive Semi-Definite Cholesky */ \
  ElError ElHPSDCholesky_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HPSDCholesky( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHPSDCholeskyDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HPSDCholesky( CReflect(uplo), *CReflect(A) ) ) } \
  /* Solve after a Cholesky factorization (without pivoting) */ \
  ElError ElSolveAfterCholesky_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( \
      cholesky::SolveAfter( \
        CReflect(uplo), CReflect(orientation), \
        *CReflect(A), *CReflect(B) ) ) } \
  ElError ElSolveAfterCholeskyDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( \
      cholesky::SolveAfter( \
        CReflect(uplo), CReflect(orientation), \
        *CReflect(A), *CReflect(B) ) ) } \
  /* Solve after a Cholesky factorization (full pivoting) */ \
  ElError ElSolveAfterCholeskyPiv_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstPermutation p, ElMatrix_ ## SIG B ) \
  { EL_TRY( \
      cholesky::SolveAfter( \
        CReflect(uplo), CReflect(orientation), \
        *CReflect(A), *CReflect(p), *CReflect(B) ) ) } \
  ElError ElSolveAfterCholeskyPivDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistPermutation p, \
    ElDistMatrix_ ## SIG B ) \
  { EL_TRY( \
      cholesky::SolveAfter( \
        CReflect(uplo), CReflect(orientation), \
        *CReflect(A), *CReflect(p), *CReflect(B) ) ) } \
  /* Generalized QR
     ============== */ \
  ElError ElGQR_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG tA, ElMatrix_ ## SIGBASE dA, \
    ElMatrix_ ## SIG B, ElMatrix_ ## SIG tB, ElMatrix_ ## SIGBASE dB ) \
  { EL_TRY( GQR( *CReflect(A), *CReflect(tA), *CReflect(dA), \
                 *CReflect(B), *CReflect(tB), *CReflect(dB) ) ) } \
  ElError ElGQRDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG tA, ElDistMatrix_ ## SIGBASE dA, \
    ElDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG tB, ElDistMatrix_ ## SIGBASE dB ) \
  { EL_TRY( GQR( *CReflect(A), *CReflect(tA), *CReflect(dA), \
                 *CReflect(B), *CReflect(tB), *CReflect(dB) ) ) } \
  ElError ElGQRExplicitTriang_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( gqr::ExplicitTriang( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElGQRExplicitTriangDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( gqr::ExplicitTriang( *CReflect(A), *CReflect(B) ) ) } \
  /* Generalized RQ
     ============== */ \
  ElError ElGRQ_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG tA, ElMatrix_ ## SIGBASE dA, \
    ElMatrix_ ## SIG B, ElMatrix_ ## SIG tB, ElMatrix_ ## SIGBASE dB ) \
  { EL_TRY( GRQ( *CReflect(A), *CReflect(tA), *CReflect(dA), \
                 *CReflect(B), *CReflect(tB), *CReflect(dB) ) ) } \
  ElError ElGRQDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG tA, ElDistMatrix_ ## SIGBASE dA, \
    ElDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG tB, ElDistMatrix_ ## SIGBASE dB ) \
  { EL_TRY( GRQ( *CReflect(A), *CReflect(tA), *CReflect(dA), \
                 *CReflect(B), *CReflect(tB), *CReflect(dB) ) ) } \
  ElError ElGRQExplicitTriang_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( grq::ExplicitTriang( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElGRQExplicitTriangDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( grq::ExplicitTriang( *CReflect(A), *CReflect(B) ) ) } \
  /* Interpolative Decomposition
     =========================== */ \
  ElError ElID_ ## SIG \
  ( ElMatrix_ ## SIG A, ElPermutation Omega, ElMatrix_ ## SIG Z, \
    ElQRCtrl_ ## SIGBASE ctrl, bool canOverwrite ) \
  { EL_TRY( \
      ID( *CReflect(A), *CReflect(Omega), *CReflect(Z), \
          CReflect(ctrl), canOverwrite ) ) } \
  ElError ElIDDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistPermutation Omega, ElDistMatrix_ ## SIG Z, \
    ElQRCtrl_ ## SIGBASE ctrl, bool canOverwrite ) \
  { EL_TRY( \
      ID( *CReflect(A), *CReflect(Omega), *CReflect(Z), \
          CReflect(ctrl), canOverwrite ) ) } \
  /* LDL factorization
     ================= */ \
  /* Return the inertia given diagonal and subdiagonal from an LDL^H fact */ \
  ElError ElInertiaAfterLDL_ ## SIG \
  ( ElConstMatrix_ ## SIGBASE d, ElConstMatrix_ ## SIG dSub, \
    ElInertiaType* inertia ) \
  { EL_TRY( *inertia = \
      CReflect(ldl::Inertia(*CReflect(d),*CReflect(dSub))) ) } \
  ElError ElInertiaAfterLDLDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIGBASE d, ElConstDistMatrix_ ## SIG dSub, \
    ElInertiaType* inertia ) \
  { EL_TRY( *inertia = CReflect(ldl::Inertia( \
      *CReflect(d), *CReflect(dSub))) ) } \
  /* LQ factorization
     ================ */ \
  /* Return the packed LQ factorization */ \
  ElError ElLQ_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG t, ElMatrix_ ## SIGBASE d ) \
  { EL_TRY( LQ( *CReflect(A), *CReflect(t), *CReflect(d) ) ) } \
  ElError ElLQDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG t, \
    ElDistMatrix_ ## SIGBASE d ) \
  { EL_TRY( LQ( *CReflect(A), *CReflect(t), *CReflect(d) ) ) } \
  /* Explicitly return both factors */ \
  ElError ElLQExplicit_ ## SIG \
  ( ElMatrix_ ## SIG L, ElMatrix_ ## SIG A ) \
  { EL_TRY( lq::Explicit( *CReflect(L), *CReflect(A) ) ) } \
  ElError ElLQExplicitDist_ ## SIG \
  ( ElDistMatrix_ ## SIG L, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( lq::Explicit( *CReflect(L), *CReflect(A) ) ) } \
  /* Only return the triangular factor */ \
  ElError ElLQExplicitTriang_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( lq::ExplicitTriang( *CReflect(A) ) ) } \
  ElError ElLQExplicitTriangDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( lq::ExplicitTriang( *CReflect(A) ) ) } \
  /* Only return the unitary factor */ \
  ElError ElLQExplicitUnitary_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( lq::ExplicitUnitary( *CReflect(A) ) ) } \
  ElError ElLQExplicitUnitaryDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( lq::ExplicitUnitary( *CReflect(A) ) ) } \
  /* Apply Q after an LQ factorization */ \
  ElError ElApplyQAfterLQ_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG t, \
    ElConstMatrix_ ## SIGBASE d, ElMatrix_ ## SIG B ) \
  { EL_TRY( lq::ApplyQ( \
      CReflect(side), CReflect(orientation), \
      *CReflect(A), *CReflect(t), \
      *CReflect(d), *CReflect(B) ) ) } \
  ElError ElApplyQAfterLQDist_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG t, \
    ElConstDistMatrix_ ## SIGBASE d, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( lq::ApplyQ( \
      CReflect(side), CReflect(orientation), \
      *CReflect(A), *CReflect(t), *CReflect(d), *CReflect(B) ) ) } \
  /* Solve against vectors after LQ factorization */ \
  ElError ElSolveAfterLQ_ ## SIG \
  ( ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG t, \
    ElConstMatrix_ ## SIGBASE d, ElConstMatrix_ ## SIG B, ElMatrix_ ## SIG X ) \
  { EL_TRY( lq::SolveAfter( \
      CReflect(orientation), *CReflect(A), *CReflect(t), \
      *CReflect(d), *CReflect(B), *CReflect(X) ) ) } \
  ElError ElSolveAfterLQDist_ ## SIG \
  ( ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG t, \
    ElConstDistMatrix_ ## SIGBASE d, ElConstDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG X ) \
  { EL_TRY( lq::SolveAfter( \
      CReflect(orientation), *CReflect(A), *CReflect(t), \
      *CReflect(d), *CReflect(B), *CReflect(X) ) ) } \
  /* LU factorization
     ================ */ \
  /* LU without pivoting */ \
  ElError ElLU_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( LU( *CReflect(A) ) ) } \
  ElError ElLUDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( LU( *CReflect(A) ) ) } \
  /* LU with partial pivoting */ \
  ElError ElLUPartialPiv_ ## SIG ( ElMatrix_ ## SIG A, ElPermutation P ) \
  { EL_TRY( LU( *CReflect(A), *CReflect(P) ) ) } \
  ElError ElLUPartialPivDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistPermutation P ) \
  { EL_TRY( LU( *CReflect(A), *CReflect(P) ) ) } \
  /* LU with full pivoting */ \
  ElError ElLUFullPiv_ ## SIG \
  ( ElMatrix_ ## SIG A, ElPermutation P, ElPermutation Q ) \
  { EL_TRY( LU( *CReflect(A), *CReflect(P), *CReflect(Q) ) ) } \
  ElError ElLUFullPivDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistPermutation P, ElDistPermutation Q ) \
  { EL_TRY( LU( *CReflect(A), *CReflect(P), *CReflect(Q) ) ) } \
  /* Solve against vectors after LU with no pivoting */ \
  ElError ElSolveAfterLU_ ## SIG \
  ( ElOrientation orientation, ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( lu::SolveAfter( \
      CReflect(orientation), *CReflect(A), *CReflect(B) ) ) } \
  ElError ElSolveAfterLUDist_ ## SIG \
  ( ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( lu::SolveAfter( \
      CReflect(orientation), *CReflect(A), *CReflect(B) ) ) } \
  /* Solve against vectors after LU with partial pivoting */ \
  ElError ElSolveAfterLUPartialPiv_ ## SIG \
  ( ElOrientation orientation, ElConstMatrix_ ## SIG A, ElConstPermutation P, \
    ElMatrix_ ## SIG B ) \
  { EL_TRY( lu::SolveAfter( \
      CReflect(orientation), \
      *CReflect(A), *CReflect(P), *CReflect(B) ) ) } \
  ElError ElSolveAfterLUPartialPivDist_ ## SIG \
  ( ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistPermutation P, \
    ElDistMatrix_ ## SIG B ) \
  { EL_TRY( lu::SolveAfter( \
      CReflect(orientation), *CReflect(A), *CReflect(P), *CReflect(B) ) ) } \
  /* Solve against vectors after LU with full pivoting */ \
  ElError ElSolveAfterLUFullPiv_ ## SIG \
  ( ElOrientation orientation, ElConstMatrix_ ## SIG A, \
    ElConstPermutation P, ElConstPermutation Q, ElMatrix_ ## SIG B ) \
  { EL_TRY( lu::SolveAfter( \
      CReflect(orientation), *CReflect(A), \
      *CReflect(P), *CReflect(Q), *CReflect(B) ) ) } \
  ElError ElSolveAfterLUFullPivDist_ ## SIG \
  ( ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, \
    ElConstDistPermutation P, ElConstDistPermutation Q, \
    ElDistMatrix_ ## SIG B ) \
  { EL_TRY( lu::SolveAfter( \
      CReflect(orientation), *CReflect(A), \
      *CReflect(P), *CReflect(Q), *CReflect(B) ) ) } \
  /* QR factorization
     ================ */ \
  /* Return the packed QR factorization (with no pivoting) */ \
  ElError ElQR_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG t, ElMatrix_ ## SIGBASE d ) \
  { EL_TRY( QR( *CReflect(A), *CReflect(t), *CReflect(d) ) ) } \
  ElError ElQRDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG t, ElDistMatrix_ ## SIGBASE d ) \
  { EL_TRY( QR( *CReflect(A), *CReflect(t), *CReflect(d) ) ) } \
  /* Return the packed QR factorization (with column pivoting) */ \
  ElError ElQRColPiv_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG t, ElMatrix_ ## SIGBASE d, \
    ElPermutation Omega ) \
  { EL_TRY( QR( \
      *CReflect(A), *CReflect(t), *CReflect(d), *CReflect(Omega) ) ) } \
  ElError ElQRColPivDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG t, ElDistMatrix_ ## SIGBASE d, \
    ElDistPermutation Omega ) \
  { EL_TRY( QR( *CReflect(A), *CReflect(t), *CReflect(d), \
                *CReflect(Omega) ) ) } \
  /* Return the packed QR factorization (with column pivoting, eXpert) */ \
  ElError ElQRColPivX_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG t, ElMatrix_ ## SIGBASE d, \
    ElPermutation Omega, ElQRCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( QR( \
      *CReflect(A), *CReflect(t), *CReflect(d), *CReflect(Omega),\
      CReflect(ctrl) ) ) } \
  ElError ElQRColPivXDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG t, ElDistMatrix_ ## SIGBASE d, \
    ElDistPermutation Omega, \
    ElQRCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( QR( \
      *CReflect(A), *CReflect(t), *CReflect(d), \
      *CReflect(Omega), CReflect(ctrl) ) ) } \
  /* Explicitly return Q and R (with no pivoting) */ \
  ElError ElQRExplicit_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG R ) \
  { EL_TRY( qr::Explicit( *CReflect(A), *CReflect(R) ) ) } \
  ElError ElQRExplicitDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG R ) \
  { EL_TRY( qr::Explicit( *CReflect(A), *CReflect(R) ) ) } \
  /* Explicitly return Q, R, and Omega (with column pivoting) */ \
  ElError ElQRColPivExplicit_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG R, ElMatrix_i Omega ) \
  { EL_TRY( qr::Explicit( *CReflect(A), *CReflect(R), *CReflect(Omega) ) ) } \
  ElError ElQRColPivExplicitDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG R, ElDistMatrix_i Omega ) \
  { EL_TRY( qr::Explicit( *CReflect(A), *CReflect(R), *CReflect(Omega) ) ) } \
  /* Return the triangular factor from QR */ \
  ElError ElQRExplicitTriang_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( qr::ExplicitTriang( *CReflect(A) ) ) } \
  ElError ElQRExplicitTriangDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( qr::ExplicitTriang( *CReflect(A) ) ) } \
  /* Return the unitary factor from QR */ \
  ElError ElQRExplicitUnitary_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( qr::ExplicitUnitary( *CReflect(A) ) ) } \
  ElError ElQRExplicitUnitaryDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( qr::ExplicitUnitary( *CReflect(A) ) ) } \
  /* Cholesky-based QR factorization */ \
  ElError ElCholeskyQR_ ## SIG ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG R ) \
  { EL_TRY( qr::Cholesky( *CReflect(A), *CReflect(R) ) ) } \
  ElError ElCholeskyQRDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG R ) \
  { EL_TRY( qr::Cholesky( *CReflect(A), *CReflect(R) ) ) } \
  /* Apply Q after a QR factorization */ \
  ElError ElApplyQAfterQR_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG t, \
    ElConstMatrix_ ## SIGBASE d, ElMatrix_ ## SIG B ) \
  { EL_TRY( qr::ApplyQ( \
      CReflect(side), CReflect(orientation), \
      *CReflect(A), *CReflect(t), \
      *CReflect(d), *CReflect(B) ) ) } \
  ElError ElApplyQAfterQRDist_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG t, \
    ElConstDistMatrix_ ## SIGBASE d, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( qr::ApplyQ( \
      CReflect(side), CReflect(orientation), \
      *CReflect(A), *CReflect(t), *CReflect(d), *CReflect(B) ) ) } \
  /* Solve against vectors after QR factorization */ \
  ElError ElSolveAfterQR_ ## SIG \
  ( ElOrientation orientation, ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG t, ElConstMatrix_ ## SIGBASE d, \
    ElConstMatrix_ ## SIG B, ElMatrix_ ## SIG X ) \
  { EL_TRY( qr::SolveAfter( \
      CReflect(orientation), *CReflect(A), \
      *CReflect(t), *CReflect(d), *CReflect(B), *CReflect(X) ) ) } \
  ElError ElSolveAfterQRDist_ ## SIG \
  ( ElOrientation orientation, ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG t, ElConstDistMatrix_ ## SIGBASE d, \
    ElConstDistMatrix_ ## SIG B, ElDistMatrix_ ## SIG X ) \
  { EL_TRY( qr::SolveAfter( \
      CReflect(orientation), *CReflect(A), \
      *CReflect(t), *CReflect(d), *CReflect(B), *CReflect(X) ) ) } \
  /* RQ factorization
     ================ */ \
  /* Return the packed RQ factorization */ \
  ElError ElRQ_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG t, ElMatrix_ ## SIGBASE d ) \
  { EL_TRY( RQ( *CReflect(A), *CReflect(t), *CReflect(d) ) ) } \
  ElError ElRQDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG t, \
    ElDistMatrix_ ## SIGBASE d ) \
  { EL_TRY( RQ( *CReflect(A), *CReflect(t), *CReflect(d) ) ) } \
  /* TODO: Explicitly return both factors */ \
  /* Only return the triangular factor */ \
  ElError ElRQExplicitTriang_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( rq::ExplicitTriang( *CReflect(A) ) ) } \
  ElError ElRQExplicitTriangDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( rq::ExplicitTriang( *CReflect(A) ) ) } \
  /* TODO: Only return the unitary factor */ \
  /* Apply Q after an RQ factorization */ \
  ElError ElApplyQAfterRQ_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG t, \
    ElConstMatrix_ ## SIGBASE d, ElMatrix_ ## SIG B ) \
  { EL_TRY( rq::ApplyQ( \
      CReflect(side), CReflect(orientation), \
      *CReflect(A), *CReflect(t), *CReflect(d), *CReflect(B) ) ) } \
  ElError ElApplyQAfterRQDist_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG t, \
    ElConstDistMatrix_ ## SIGBASE d, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( rq::ApplyQ( \
      CReflect(side), CReflect(orientation), \
      *CReflect(A), *CReflect(t), *CReflect(d), *CReflect(B) ) ) } \
  /* Solve against vectors after RQ factorization */ \
  ElError ElSolveAfterRQ_ ## SIG \
  ( ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG t, \
    ElConstMatrix_ ## SIGBASE d, ElConstMatrix_ ## SIG B, ElMatrix_ ## SIG X ) \
  { EL_TRY( rq::SolveAfter( \
      CReflect(orientation), *CReflect(A), *CReflect(t), \
      *CReflect(d), *CReflect(B), *CReflect(X) ) ) } \
  ElError ElSolveAfterRQDist_ ## SIG \
  ( ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG t, \
    ElConstDistMatrix_ ## SIGBASE d, ElConstDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG X ) \
  { EL_TRY( rq::SolveAfter( \
      CReflect(orientation), *CReflect(A), *CReflect(t), \
      *CReflect(d), *CReflect(B), *CReflect(X) ) ) } \
  /* Skeleton factorization
     ====================== */ \
  ElError ElSkeleton_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElPermutation PR, ElPermutation PC, \
    ElMatrix_ ## SIG Z, ElQRCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( Skeleton( *CReflect(A), *CReflect(PR), *CReflect(PC), \
                      *CReflect(Z), CReflect(ctrl) ) ) } \
  ElError ElSkeletonDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistPermutation PR, ElDistPermutation PC, \
    ElDistMatrix_ ## SIG Z, ElQRCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( Skeleton( *CReflect(A), *CReflect(PR), *CReflect(PC), \
                      *CReflect(Z), CReflect(ctrl) ) ) }

#define C_PROTO_REAL(SIG,Real) \
  C_PROTO_FIELD(SIG,SIG,Real) \
  /* LDL factorization
     ================= */ \
  ElError ElLDLPivotConstant_ ## SIG ( ElLDLPivotType pivotType, Real* gamma ) \
  { EL_TRY( *gamma = LDLPivotConstant<Real>(CReflect(pivotType)) ) } \
  /* Return the packed LDL factorization (without pivoting) */ \
  ElError ElLDL_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( LDL( *CReflect(A), false ) ) } \
  ElError ElLDLDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( LDL( *CReflect(A), false ) ) } \
  /* Return the packed LDL factorization with pivoting */ \
  ElError ElLDLPiv_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG dSub, ElPermutation p ) \
  { EL_TRY( LDL( *CReflect(A), *CReflect(dSub), *CReflect(p), false ) ) } \
  ElError ElLDLPivDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG dSub, ElDistPermutation p ) \
  { EL_TRY( LDL( *CReflect(A), *CReflect(dSub), *CReflect(p), false ) ) } \
  /* Expert versions */ \
  ElError ElLDLPivX_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG dSub, ElPermutation p, \
    ElLDLPivotCtrl_ ## SIG ctrl ) \
  { EL_TRY( LDL( *CReflect(A), *CReflect(dSub), *CReflect(p), false, \
                 CReflect(ctrl) ) ) } \
  ElError ElLDLPivXDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG dSub, ElDistPermutation p, \
    ElLDLPivotCtrl_ ## SIG ctrl ) \
  { EL_TRY( LDL( *CReflect(A), *CReflect(dSub), *CReflect(p), false, \
                 CReflect(ctrl) ) ) } \
  /* Multiply vectors after an unpivoted LDL factorization */ \
  ElError ElMultiplyAfterLDL_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( ldl::MultiplyAfter( *CReflect(A), *CReflect(B), false ) ) } \
  ElError ElMultiplyAfterLDLDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( ldl::MultiplyAfter( *CReflect(A), *CReflect(B), false ) ) } \
  /* Multiply vectors after a pivoted LDL factorization */ \
  ElError ElMultiplyAfterLDLPiv_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG dSub, ElConstPermutation p, \
    ElMatrix_ ## SIG B ) \
  { EL_TRY( ldl::MultiplyAfter( \
      *CReflect(A), *CReflect(dSub), *CReflect(p), *CReflect(B), \
      false ) ) } \
  ElError ElMultiplyAfterLDLPivDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG dSub, \
    ElConstDistPermutation p, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( ldl::MultiplyAfter( \
      *CReflect(A), *CReflect(dSub), *CReflect(p), *CReflect(B), false ) ) } \
  /* Solve against vectors after an unpivoted LDL factorization */ \
  ElError ElSolveAfterLDL_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( ldl::SolveAfter( *CReflect(A), *CReflect(B), false ) ) } \
  ElError ElSolveAfterLDLDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( ldl::SolveAfter( *CReflect(A), *CReflect(B), false ) ) } \
  /* Solve against vectors after a pivoted LDL factorization */ \
  ElError ElSolveAfterLDLPiv_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG dSub, ElConstPermutation p, \
    ElMatrix_ ## SIG B ) \
  { EL_TRY( ldl::SolveAfter( \
      *CReflect(A), *CReflect(dSub), *CReflect(p), *CReflect(B), \
      false ) ) } \
  ElError ElSolveAfterLDLPivDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG dSub, \
    ElConstDistPermutation p, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( ldl::SolveAfter( \
      *CReflect(A), *CReflect(dSub), *CReflect(p), *CReflect(B), false ) ) } \
  /* LU factorization
     ================ */ \
  /* Rank-one LU factorization modification */ \
  ElError ElLUMod_ ## SIG \
  ( ElMatrix_ ## SIG A, ElPermutation P, \
    ElConstMatrix_ ## SIG u, ElConstMatrix_ ## SIG v, Real tau ) \
  { EL_TRY( LUMod( \
      *CReflect(A), *CReflect(P), \
      *CReflect(u), *CReflect(v), false, tau ) ) } \
  ElError ElLUModDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistPermutation P, \
    ElConstDistMatrix_ ## SIG u, ElConstDistMatrix_ ## SIG v, Real tau ) \
  { EL_TRY( LUMod( \
      *CReflect(A), *CReflect(P), \
      *CReflect(u), *CReflect(v), false, tau ) ) }

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* LDL factorization
     ================= */ \
  /* Return the packed LDL factorization (without pivoting) */ \
  ElError ElLDL_ ## SIG ( ElMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( LDL( *CReflect(A), conjugate ) ) } \
  ElError ElLDLDist_ ## SIG ( ElDistMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( LDL( *CReflect(A), conjugate ) ) } \
  /* Return the packed LDL factorization with pivoting */ \
  ElError ElLDLPiv_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElMatrix_ ## SIG dSub, \
    ElPermutation p, \
    bool conjugate ) \
  { EL_TRY( LDL( *CReflect(A), *CReflect(dSub), *CReflect(p), conjugate ) ) } \
  ElError ElLDLPivDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG dSub, \
    ElDistPermutation p, \
    bool conjugate ) \
  { EL_TRY( LDL( *CReflect(A), *CReflect(dSub), *CReflect(p), conjugate ) ) } \
  /* Expert versions */ \
  ElError ElLDLPivX_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElMatrix_ ## SIG dSub, \
    ElPermutation p, \
    bool conjugate, \
    ElLDLPivotCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( LDL( *CReflect(A), *CReflect(dSub), *CReflect(p), conjugate, \
                 CReflect(ctrl) ) ) } \
  ElError ElLDLPivXDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG dSub, \
    ElDistPermutation p, \
    bool conjugate, ElLDLPivotCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( LDL( *CReflect(A), *CReflect(dSub), *CReflect(p), conjugate, \
                 CReflect(ctrl) ) ) } \
  /* Multiply vectors after an unpivoted LDL factorization */ \
  ElError ElMultiplyAfterLDL_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B, bool conjugate ) \
  { EL_TRY( ldl::MultiplyAfter( *CReflect(A), *CReflect(B), conjugate ) ) } \
  ElError ElMultiplyAfterLDLDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B, bool conjugate ) \
  { EL_TRY( ldl::MultiplyAfter( *CReflect(A), *CReflect(B), conjugate ) ) } \
  /* Multiply vectors after a pivoted LDL factorization */ \
  ElError ElMultiplyAfterLDLPiv_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG dSub, \
    ElConstPermutation p, \
    ElMatrix_ ## SIG B, bool conjugate ) \
  { EL_TRY( ldl::MultiplyAfter( \
      *CReflect(A), *CReflect(dSub), *CReflect(p), *CReflect(B), \
      conjugate ) ) } \
  ElError ElMultiplyAfterLDLPivDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG dSub, \
    ElConstDistPermutation p, \
    ElDistMatrix_ ## SIG B, \
    bool conjugate ) \
  { EL_TRY( ldl::MultiplyAfter( \
      *CReflect(A), *CReflect(dSub), *CReflect(p), *CReflect(B), \
      conjugate ) ) } \
  /* Solve against vectors after an unpivoted LDL factorization */ \
  ElError ElSolveAfterLDL_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B, bool conjugate ) \
  { EL_TRY( ldl::SolveAfter( *CReflect(A), *CReflect(B), conjugate ) ) } \
  ElError ElSolveAfterLDLDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B, bool conjugate ) \
  { EL_TRY( ldl::SolveAfter( *CReflect(A), *CReflect(B), conjugate ) ) } \
  /* Solve against vectors after a pivoted LDL factorization */ \
  ElError ElSolveAfterLDLPiv_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG dSub, \
    ElConstPermutation p, \
    ElMatrix_ ## SIG B, \
    bool conjugate ) \
  { EL_TRY( ldl::SolveAfter( \
      *CReflect(A), *CReflect(dSub), *CReflect(p), *CReflect(B), \
      conjugate ) ) } \
  ElError ElSolveAfterLDLPivDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG dSub, \
    ElConstDistPermutation p, \
    ElDistMatrix_ ## SIG B, \
    bool conjugate ) \
  { EL_TRY( ldl::SolveAfter( \
      *CReflect(A), *CReflect(dSub), *CReflect(p), *CReflect(B), \
      conjugate ) ) } \
  /* LU factorization
     ================ */ \
  /* Rank-one LU factorization modification */ \
  ElError ElLUMod_ ## SIG \
  ( ElMatrix_ ## SIG A, ElPermutation P, \
    ElConstMatrix_ ## SIG u, ElConstMatrix_ ## SIG v, \
    bool conjugate, Base<F> tau ) \
  { EL_TRY( LUMod( \
      *CReflect(A), *CReflect(P), \
      *CReflect(u), *CReflect(v), conjugate, tau ) ) } \
  ElError ElLUModDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistPermutation P, \
    ElConstDistMatrix_ ## SIG u, ElConstDistMatrix_ ## SIG v, \
    bool conjugate, Base<F> tau ) \
  { EL_TRY( LUMod( \
      *CReflect(A), *CReflect(P), \
      *CReflect(u), *CReflect(v), conjugate, tau ) ) }

#define EL_NO_INT_PROTO
#include <El/macros/CInstantiate.h>

} // extern "C"
