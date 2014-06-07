/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FACTOR_HPP
#define EL_FACTOR_HPP

namespace El {

template<typename Real>
struct QRCtrl
{
    bool boundRank;  
    Int maxRank;

    bool adaptive;
    Real tol;

    bool alwaysRecomputeNorms;

    // TODO: Add Chan ratio

    QRCtrl()
    : boundRank(false), maxRank(0), adaptive(false), tol(0),
      alwaysRecomputeNorms(false)
    { }
};

// Cholesky
// ========
template<typename F>
void Cholesky( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void Cholesky( UpperOrLower uplo, DistMatrix<F>& A );

template<typename F>
void ReverseCholesky( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void ReverseCholesky( UpperOrLower uplo, DistMatrix<F>& A );

template<typename F>
void Cholesky( UpperOrLower uplo, Matrix<F>& A, Matrix<Int>& pPerm );
template<typename F,Dist UPerm>
void Cholesky
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& pPerm );

template<typename F>
void CholeskyMod
( UpperOrLower uplo, Matrix<F>& T, Base<F> alpha, Matrix<F>& V );
template<typename F>
void CholeskyMod
( UpperOrLower uplo, DistMatrix<F>& T, Base<F> alpha, DistMatrix<F>& V );

template<typename F>
void HPSDCholesky( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void HPSDCholesky( UpperOrLower uplo, DistMatrix<F>& A );

namespace cholesky {

template<typename F>
void SolveAfter
( UpperOrLower uplo, Orientation orientation, 
  const Matrix<F>& A, Matrix<F>& B );
template<typename F>
void SolveAfter
( UpperOrLower uplo, Orientation orientation,
  const DistMatrix<F>& A, DistMatrix<F>& B );

template<typename F>
void SolveAfter
( UpperOrLower uplo, Orientation orientation,
  const Matrix<F>& A, const Matrix<Int>& pPerm, 
        Matrix<F>& B );
template<typename F,Dist UPerm>
void SolveAfter
( UpperOrLower uplo, Orientation orientation,
  const DistMatrix<F>& A, const DistMatrix<Int,UPerm,STAR>& pPerm,
        DistMatrix<F>& B );

} // namespace cholesky

// Generalized QR
// ==============
template<typename F>
void GQR( Matrix<F>& A, Matrix<F>& B );
template<typename F>
void GQR( DistMatrix<F>& A, DistMatrix<F>& B );

template<typename F>
void GQR
( Matrix<F>& A, Matrix<F>& tA, Matrix<Base<F>>& dA,
  Matrix<F>& B, Matrix<F>& tB, Matrix<Base<F>>& dB );
template<typename F>
void GQR
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& tA, DistMatrix<Base<F>,MD,STAR>& dA,
  DistMatrix<F>& B, DistMatrix<F,MD,STAR>& tB, DistMatrix<Base<F>,MD,STAR>& dB 
);

// Generalized RQ
// ==============
template<typename F>
void GRQ( Matrix<F>& A, Matrix<F>& B );
template<typename F>
void GRQ( DistMatrix<F>& A, DistMatrix<F>& B );

template<typename F>
void GRQ
( Matrix<F>& A, Matrix<F>& tA, Matrix<Base<F>>& dA,
  Matrix<F>& B, Matrix<F>& tB, Matrix<Base<F>>& dB );
template<typename F>
void GRQ
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& tA, DistMatrix<Base<F>,MD,STAR>& dA,
  DistMatrix<F>& B, DistMatrix<F,MD,STAR>& tB, DistMatrix<Base<F>,MD,STAR>& dB 
);

// Interpolative Decomposition
// ===========================
template<typename F>
void ID
( const Matrix<F>& A, Matrix<Int>& pPerm, 
        Matrix<F>& Z, 
  const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>() );
// NOTE: This is only instantiated for UPerm=VR
template<typename F,Dist UPerm>
void ID
( const DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& pPerm,
        DistMatrix<F,STAR,VR>& Z, 
  const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>() );

template<typename F>
void ID
( Matrix<F>& A, Matrix<Int>& pPerm, 
  Matrix<F>& Z, const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>(), 
  bool canOverwrite=false );
// NOTE: This is only instantiated for UPerm=VR
template<typename F,Dist UPerm>
void ID
( DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& pPerm,
  DistMatrix<F,STAR,VR>& Z, const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>(), 
  bool canOverwrite=false );

// LDL
// ===
namespace LDLPivotTypeNS {
enum LDLPivotType
{
    BUNCH_KAUFMAN_A=0,
    BUNCH_KAUFMAN_C=1,
    BUNCH_KAUFMAN_D=2,
    BUNCH_KAUFMAN_BOUNDED=3,
    BUNCH_PARLETT=4
};
}
using namespace LDLPivotTypeNS;

struct LDLPivot
{
    Int nb;
    Int from[2];
};

// Return the L from an LDL factorization of A
// -------------------------------------------
template<typename F>
void LDL( Matrix<F>& A, bool conjugate );
template<typename F>
void LDL( DistMatrix<F>& A, bool conjugate );

// Return an implicit representation of a pivoted LDL factorization of A
// ---------------------------------------------------------------------
template<typename F>
void LDL
( Matrix<F>& A, Matrix<F>& dSub,
  Matrix<Int>& pPerm, bool conjugate,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A );
// NOTE: Only instantiated for UPerm=VC
template<typename F,Dist UPerm>
void LDL
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& dSub,
  DistMatrix<Int,UPerm,STAR>& pPerm, bool conjugate,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A );

namespace ldl {

// Compute the inertia triplet of a Hermitian matrix's LDL^H factorization
// -----------------------------------------------------------------------
template<typename F>
InertiaType Inertia
( const Matrix<Base<F>>& d, const Matrix<F>& dSub );
// NOTE: Only instantiated for (U,V) = (MD,STAR)
template<typename F,Dist U,Dist V>
InertiaType Inertia
( const DistMatrix<Base<F>,U,V>& d, const DistMatrix<F,U,V>& dSub );

// Multiply vectors using an implicit representation of an LDL factorization
// -------------------------------------------------------------------------
template<typename F>
void MultiplyAfter( const Matrix<F>& A, Matrix<F>& B, bool conjugated );
template<typename F>
void MultiplyAfter( const DistMatrix<F>& A, DistMatrix<F>& B, bool conjugated );

// Multiply vectors using an implicit representation of a pivoted LDL fact.
// ------------------------------------------------------------------------
template<typename F>
void MultiplyAfter
( const Matrix<F>& A, const Matrix<F>& dSub,
  const Matrix<Int>& pPerm, Matrix<F>& B, bool conjugated );
// NOTE: Only instantiated for UPerm=VR
template<typename F,Dist UPerm>
void MultiplyAfter
( const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& dSub,
  const DistMatrix<Int,UPerm,STAR>& pPerm, DistMatrix<F>& B, bool conjugated );

// Solve linear systems using an implicit LDL factorization
// --------------------------------------------------------
template<typename F>
void SolveAfter( const Matrix<F>& A, Matrix<F>& B, bool conjugated );
template<typename F>
void SolveAfter( const DistMatrix<F>& A, DistMatrix<F>& B, bool conjugated );

// Solve linear system with the implicit representations of L, D, and P
// --------------------------------------------------------------------
template<typename F>
void SolveAfter
( const Matrix<F>& A, const Matrix<F>& dSub,
  const Matrix<Int>& pPerm, Matrix<F>& B, bool conjugated );
// NOTE: Only instantiated for UPerm=VR
template<typename F,Dist UPerm>
void SolveAfter
( const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& dSub,
  const DistMatrix<Int,UPerm,STAR>& pPerm, DistMatrix<F>& B, bool conjugated );

} // namespace ldl

// LQ
// ==

// Overwrite A with L
// ------------------
template<typename F>
void LQ( Matrix<F>& A );
template<typename F>
void LQ( DistMatrix<F>& A );

// Overwrite A with both L and the scaled Householder vectors
// ----------------------------------------------------------
template<typename F>
void LQ( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d );
template<typename F>
void LQ
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d );

namespace lq {

// Apply Q using its implicit representation
// -----------------------------------------
template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& t,
  const Matrix<Base<F>>& d, Matrix<F>& B );
// NOTE: Only instantiated for (Ut,Vt)=(Ud,Vd)=(MD,STAR) 
template<typename F,Dist Ut,Dist Vt,Dist Ud,Dist Vd>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const DistMatrix<F>& A, const DistMatrix<F,Ut,Vt>& t,
  const DistMatrix<Base<F>,Ud,Vd>& d, DistMatrix<F>& B );

// Overwrite A with Q
// ------------------
template<typename F>
void Explicit( Matrix<F>& A );
template<typename F>
void Explicit( DistMatrix<F>& A );

// Return both L and Q such that A = L Q
// -------------------------------------
template<typename F>
void Explicit( Matrix<F>& L, Matrix<F>& A );
template<typename F>
void Explicit( DistMatrix<F>& L, DistMatrix<F>& A );

// Solve a linear system with the implicit representations of L and Q 
// ------------------------------------------------------------------
template<typename F>
void SolveAfter
( Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& t,
  const Matrix<Base<F>>& d, const Matrix<F>& B,
        Matrix<F>& X );
template<typename F>
void SolveAfter
( Orientation orientation,
  const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t,
  const DistMatrix<Base<F>,MD,STAR>& d, const DistMatrix<F>& B,
        DistMatrix<F>& X );

} // namespace lq

// LU
// ==

// LU without pivoting
// -------------------
template<typename F>
void LU( Matrix<F>& A );
template<typename F>
void LU( DistMatrix<F>& A );

// LU with partial pivoting
// ------------------------
template<typename F>
void LU( Matrix<F>& A, Matrix<Int>& pPerm );
// NOTE: Only instantiated for UPerm=VC
template<typename F,Dist UPerm>
void LU( DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& pPerm );

// LU with full pivoting
// ---------------------
template<typename F>
void LU( Matrix<F>& A, Matrix<Int>& pPerm, Matrix<Int>& qPerm );
// NOTE: Only instantiated for UPerm=VC
template<typename F,Dist UPerm>
void LU
( DistMatrix<F>& A, 
  DistMatrix<Int,UPerm,STAR>& pPerm, DistMatrix<Int,UPerm,STAR>& qPerm );

// Rank-one modification of a partially-pivoted LU factorization
// -------------------------------------------------------------
template<typename F>
void LUMod
( Matrix<F>& A, Matrix<Int>& perm,
  const Matrix<F>& u, const Matrix<F>& v, bool conjugate=true,
  Base<F> tau=0.1 );
// NOTE: Only instantiated for UPerm=VC
template<typename F,Dist UPerm>
void LUMod
( DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& perm,
  const DistMatrix<F>& u, const DistMatrix<F>& v, bool conjugate=true,
  Base<F> tau=0.1 );

namespace lu {

// Perform a panel factorization
// -----------------------------
// NOTE: This is currently only needed within GaussianElimination and
//       could ideally be removed
template<typename F>
void Panel( Matrix<F>& APan, Matrix<Int>& p1 );
template<typename F>
void Panel
( DistMatrix<F,  STAR,STAR>& A11, 
  DistMatrix<F,  MC,  STAR>& A21, 
  DistMatrix<Int,STAR,STAR>& p1 );

// Solve linear systems using an implicit unpivoted LU factorization
// -----------------------------------------------------------------
template<typename F>
void SolveAfter( Orientation orientation, const Matrix<F>& A, Matrix<F>& B );
template<typename F>
void SolveAfter
( Orientation orientation, const DistMatrix<F>& A, DistMatrix<F>& B );

// Solve linear systems using an implicit partially-pivoted LU factorization
// -------------------------------------------------------------------------
template<typename F>
void SolveAfter
( Orientation orientation, const Matrix<F>& A,
  const Matrix<Int>& pPerm, Matrix<F>& B );
// NOTE: Only instantiated for UPerm=VC
template<typename F,Dist UPerm>
void SolveAfter
( Orientation orientation, const DistMatrix<F>& A,
  const DistMatrix<Int,UPerm,STAR>& pPerm, DistMatrix<F>& B );

// Solve linear systems using an implicit fully-pivoted LU factorization
// ---------------------------------------------------------------------
template<typename F>
void SolveAfter
( Orientation orientation, const Matrix<F>& A,
  const Matrix<Int>& pPerm,
  const Matrix<Int>& qPerm,
        Matrix<F>& B );
// NOTE: Only instantiated for UPerm=VC
template<typename F,Dist UPerm>
void SolveAfter
( Orientation orientation, const DistMatrix<F>& A,
  const DistMatrix<Int,UPerm,STAR>& pPerm,
  const DistMatrix<Int,UPerm,STAR>& qPerm,
        DistMatrix<F>& B );

} // namespace lu

// QR
// ==

// Return R (with non-negative diagonal) such that A = Q R
// -------------------------------------------------------
template<typename F>
void QR( Matrix<F>& A );
template<typename F>
void QR( DistMatrix<F>& A );

// Return an implicit representation of Q and R such that A = Q R
// --------------------------------------------------------------
template<typename F>
void QR
( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d );
template<typename F>
void QR
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d );

// Return R and (implicit) P such that A P ~= Q R
// ----------------------------------------------
template<typename F>
Int QR
( Matrix<F>& A, Matrix<Int>& pPerm, 
  const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>() );
// NOTE: Only instantiated for UPerm=VR
template<typename F,Dist UPerm>
Int QR
( DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& pPerm, 
  const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>() );

// Return an implicit representation of (Q,R,P) such that A P ~= Q R
// -----------------------------------------------------------------
template<typename F>
Int QR
( Matrix<F>& A, Matrix<F>& t, 
  Matrix<Base<F>>& d, Matrix<Int>& pPerm,
  const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>() );
// NOTE: Only instantiated for UPerm=VR
template<typename F,Dist UPerm>
Int QR
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, 
  DistMatrix<Base<F>,MD,STAR>& d, DistMatrix<Int,UPerm,STAR>& pPerm,
  const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>() );

namespace qr {

// Apply Q using its implicit representation
// -----------------------------------------
template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& t,
  const Matrix<Base<F>>& d, Matrix<F>& B );
// NOTE: Only instantiated for (Ut,Vt)=(Ud,Vd)=(MD,STAR) 
template<typename F,Dist Ut,Dist Vt,Dist Ud,Dist Vd>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const DistMatrix<F>& A, const DistMatrix<F,Ut,Vt>& t,
  const DistMatrix<Base<F>,Ud,Vd>& d, DistMatrix<F>& B );

// Cholesky-based QR
// -----------------
template<typename F>
void Cholesky( Matrix<F>& A, Matrix<F>& R );
template<typename F>
void Cholesky( DistMatrix<F,VC,STAR>& A, DistMatrix<F,STAR,STAR>& R );

// Overwrite A with Q
// ------------------
template<typename F>
void Explicit( Matrix<F>& A, bool colPiv=false );
template<typename F>
void Explicit( DistMatrix<F>& A, bool colPiv=false );
// TODO: Version which accepts QRCtrl

// Return both Q and R such that A = Q R or A P = Q R
// --------------------------------------------------
template<typename F>
void Explicit( Matrix<F>& A, Matrix<F>& R, bool colPiv=false );
template<typename F>
void Explicit( DistMatrix<F>& A, DistMatrix<F>& R, bool colPiv=false );
// TODO: Version which accepts QRCtrl

// Return (Q,R,P) such that A P = Q R
// ----------------------------------
// TODO: Return the explicit permutation matrix instead of the representative
//       vector
template<typename F>
void Explicit( Matrix<F>& A, Matrix<F>& R, Matrix<Int>& pPerm );
// NOTE: Only instantiated for UPerm=VR
template<typename F,Dist UPerm>
void Explicit
( DistMatrix<F>& A, DistMatrix<F>& R, DistMatrix<Int,UPerm,STAR>& pPerm );
// TODO: Version which accepts QRCtrl

// Solve a linear system with the implicit QR factorization
// --------------------------------------------------------
template<typename F>
void SolveAfter
( Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& t,
  const Matrix<Base<F>>& d, const Matrix<F>& B,
        Matrix<F>& X );
template<typename F>
void SolveAfter
( Orientation orientation,
  const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t,
  const DistMatrix<Base<F>,MD,STAR>& d, const DistMatrix<F>& B,
        DistMatrix<F>& X );
// TODO: Version which involves permutation matrix

// TODO: Add support for solving with an implicit pivoted QR factorizatino

template<typename F>
struct TreeData
{
    Matrix<F> QR0, t0;
    Matrix<Base<F>> d0;
    std::vector<Matrix<F>> QRList;
    std::vector<Matrix<F>> tList;
    std::vector<Matrix<Base<F>>> dList;

    TreeData( Int numStages=0 )
    : QRList(numStages), tList(numStages), dList(numStages)
    { }

    TreeData( TreeData<F>&& treeData )
    : QR0(std::move(treeData.QR0)),
      t0(std::move(treeData.t0)),
      d0(std::move(treeData.d0)),
      QRList(std::move(treeData.QRList)),
      tList(std::move(treeData.tList)),
      dList(std::move(treeData.dList))
    { }

    TreeData<F>& operator=( TreeData<F>&& treeData )
    {
        QR0 = std::move(treeData.QR0);
        t0 = std::move(treeData.t0);
        d0 = std::move(treeData.d0);
        QRList = std::move(treeData.QRList);
        tList = std::move(treeData.tList);
        dList = std::move(treeData.dList);
        return *this;
    }
};

// Return an implicit tall-skinny QR factorization
// NOTE: Only instantiated for U=VC
template<typename F,Dist U>
TreeData<F> TS( const DistMatrix<F,U,STAR>& A );

// Return an explicit tall-skinny QR factorization
// NOTE: Only instantiated for U=VC
template<typename F,Dist U>
void ExplicitTS( DistMatrix<F,U,STAR>& A, DistMatrix<F,STAR,STAR>& R );

namespace ts {

// NOTE: Only instantiated for U=VC
template<typename F,Dist U>
Matrix<F>& RootQR
( const DistMatrix<F,U,STAR>& A, TreeData<F>& treeData );

// NOTE: Only instantiated for U=VC
template<typename F,Dist U>
const Matrix<F>& RootQR
( const DistMatrix<F,U,STAR>& A, const TreeData<F>& treeData );

// NOTE: Only instantiated for U=VC
template<typename F,Dist U>
void Reduce( const DistMatrix<F,U,STAR>& A, TreeData<F>& treeData );

// NOTE: Only instantiated for U=VC
template<typename F,Dist U>
void Scatter( DistMatrix<F,U,STAR>& A, const TreeData<F>& treeData );

} // namespace ts

} // namespace qr

// RQ
// ==
template<typename F>
void RQ( Matrix<F>& A );
template<typename F>
void RQ( DistMatrix<F>& A );

template<typename F>
void RQ( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d );
template<typename F>
void RQ
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d );

namespace rq {

template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& t,
  const Matrix<Base<F>>& d, Matrix<F>& B );
// NOTE: Only instantiated for (Ut,Vt)=(Ud,Vd)=(MD,STAR) 
template<typename F,Dist Ut,Dist Vt,Dist Ud,Dist Vd>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const DistMatrix<F>& A, const DistMatrix<F,Ut,Vt>& t,
  const DistMatrix<Base<F>,Ud,Vd>& d, DistMatrix<F>& B );

// TODO: Think about ensuring this ordering is consistent with lq::Explicit
template<typename F>
void Cholesky( Matrix<F>& A, Matrix<F>& R );
template<typename F>
void Cholesky( DistMatrix<F,STAR,VR>& A, DistMatrix<F,STAR,STAR>& R );

template<typename F>
void SolveAfter
( Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& t,
  const Matrix<Base<F>>& d, const Matrix<F>& B,
        Matrix<F>& X );
template<typename F>
void SolveAfter
( Orientation orientation,
  const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t,
  const DistMatrix<Base<F>,MD,STAR>& d, const DistMatrix<F>& B,
        DistMatrix<F>& X );

} // namespace rq

// Skeleton
// ========
template<typename F>
void Skeleton
( const Matrix<F>& A,
  Matrix<Int>& permR, Matrix<Int>& permC,
  Matrix<F>& Z, const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>() );
// NOTE: Only instantiated for UPerm=VR
template<typename F,Dist UPerm>
void Skeleton
( const DistMatrix<F>& A,
  DistMatrix<Int,UPerm,STAR>& permR, DistMatrix<Int,UPerm,STAR>& permC,
  DistMatrix<F>& Z, const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>() );

// Inline convenience functions
// ############################
// TODO: Statically instantiate routines with modest numbers of instances

// Cholesky
// ========
template<typename F>
inline void
LocalCholesky( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalCholesky"))
    Cholesky( uplo, A.Matrix() );
}

template<typename F>
inline void
LocalReverseCholesky( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalReverseCholesky"))
    ReverseCholesky( uplo, A.Matrix() );
}

// LDL
// ===
template<typename F>
inline void
LocalLDL( DistMatrix<F,STAR,STAR>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("LocalLDL"))
    LDL( A.Matrix(), conjugate );
}

// LU
// ==
template<typename F>
inline void
LocalLU( DistMatrix<F,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalLU"))
    LU( A.Matrix() );
}

} // namespace El

#endif // ifndef EL_FACTOR_HPP
