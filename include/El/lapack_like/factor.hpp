/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FACTOR_HPP
#define EL_FACTOR_HPP

#include "El/lapack_like/factor/symbolic.hpp"
#include "El/lapack_like/factor/numeric.hpp"

namespace El {

// Cholesky
// ========
template<typename F>
void Cholesky( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void Cholesky( UpperOrLower uplo, AbstractDistMatrix<F>& A );

template<typename F>
void ReverseCholesky( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void ReverseCholesky( UpperOrLower uplo, AbstractDistMatrix<F>& A );

template<typename F>
void Cholesky( UpperOrLower uplo, Matrix<F>& A, Matrix<Int>& p );
template<typename F>
void Cholesky
( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& p );

template<typename F>
void CholeskyMod
( UpperOrLower uplo, Matrix<F>& T, Base<F> alpha, Matrix<F>& V );
template<typename F>
void CholeskyMod
( UpperOrLower uplo, AbstractDistMatrix<F>& T, 
  Base<F> alpha, AbstractDistMatrix<F>& V );

template<typename F>
void HPSDCholesky( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void HPSDCholesky( UpperOrLower uplo, AbstractDistMatrix<F>& A );

namespace cholesky {

template<typename F>
void SolveAfter
( UpperOrLower uplo, Orientation orientation, 
  const Matrix<F>& A, Matrix<F>& B );
template<typename F>
void SolveAfter
( UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B );

template<typename F>
void SolveAfter
( UpperOrLower uplo, Orientation orientation,
  const Matrix<F>& A, const Matrix<Int>& p, 
        Matrix<F>& B );
template<typename F>
void SolveAfter
( UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<F>& A, const AbstractDistMatrix<Int>& p,
        AbstractDistMatrix<F>& B );

} // namespace cholesky

// LDL
// ===
namespace LDLPivotTypeNS {
enum LDLPivotType
{
    BUNCH_KAUFMAN_A,
    BUNCH_KAUFMAN_C,
    BUNCH_KAUFMAN_D,
    BUNCH_KAUFMAN_BOUNDED,
    BUNCH_PARLETT,
    LDL_WITHOUT_PIVOTING
    /* TODO: Diagonal pivoting? */
};
}
using namespace LDLPivotTypeNS;

template<typename Real>
Real LDLPivotConstant( LDLPivotType pivType )
{
    // TODO: Check that the Bunch-Parlett choice is optimal
    switch( pivType )
    {
    case BUNCH_KAUFMAN_A:
    case BUNCH_PARLETT:   return (1+Sqrt(Real(17)))/8;
    case BUNCH_KAUFMAN_D: return Real(0.525);
    default: 
        LogicError("No default constant exists for this pivot type");
        return 0;
    }
}

struct LDLPivot
{
    Int nb;
    Int from[2];
};

template<typename Real>
struct LDLPivotCtrl {
  LDLPivotType pivotType;
  Real gamma;

  LDLPivotCtrl( LDLPivotType piv=BUNCH_KAUFMAN_A ) 
  : pivotType(piv), gamma(LDLPivotConstant<Real>(piv)) { }
};

// Return the L (and D) from an LDL factorization of A (without pivoting)
// ----------------------------------------------------------------------
template<typename F>
void LDL( Matrix<F>& A, bool conjugate );
template<typename F>
void LDL( AbstractDistMatrix<F>& A, bool conjugate );

// Return an implicit representation of a pivoted LDL factorization of A
// ---------------------------------------------------------------------
template<typename F>
void LDL
( Matrix<F>& A, Matrix<F>& dSub,
  Matrix<Int>& p, bool conjugate,
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );
template<typename F>
void LDL
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& dSub,
  AbstractDistMatrix<Int>& p, bool conjugate,
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );

// All fronts of L are required to be initialized to the expansions of the 
// original sparse matrix before calling LDL.
template<typename F>
void LDL
( const SymmNodeInfo& info, SymmFront<F>& L, 
  SymmFrontType newType=LDL_2D );
template<typename F>
void LDL
( const DistSymmNodeInfo& info, DistSymmFront<F>& L, 
  SymmFrontType newType=LDL_2D );

namespace ldl {

// Compute the inertia triplet of a Hermitian matrix's LDL^H factorization
// -----------------------------------------------------------------------
template<typename F>
InertiaType Inertia
( const Matrix<Base<F>>& d, const Matrix<F>& dSub );
template<typename F>
InertiaType Inertia
( const AbstractDistMatrix<Base<F>>& d, const AbstractDistMatrix<F>& dSub );

// Multiply vectors using an implicit representation of an LDL factorization
// -------------------------------------------------------------------------
template<typename F>
void MultiplyAfter( const Matrix<F>& A, Matrix<F>& B, bool conjugated );
template<typename F>
void MultiplyAfter
( const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B, bool conjugated );

// Multiply vectors using an implicit representation of a pivoted LDL fact.
// ------------------------------------------------------------------------
template<typename F>
void MultiplyAfter
( const Matrix<F>& A, const Matrix<F>& dSub,
  const Matrix<Int>& p, Matrix<F>& B, bool conjugated );
template<typename F>
void MultiplyAfter
( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& dSub,
  const AbstractDistMatrix<Int>& p, AbstractDistMatrix<F>& B, bool conjugated );

// Solve linear systems using an implicit LDL factorization
// --------------------------------------------------------
template<typename F>
void SolveAfter( const Matrix<F>& A, Matrix<F>& B, bool conjugated );
template<typename F>
void SolveAfter
( const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B, bool conjugated );

template<typename F>
void SolveAfter
( const vector<Int>& invMap, const SymmNodeInfo& info,
  const SymmFront<F>& front, Matrix<F>& X ); 
template<typename F>
void SolveAfter
( const DistMap& invMap, const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front, DistMultiVec<F>& X ); 

template<typename F>
void SolveAfter
( const SymmNodeInfo& info,
  const SymmFront<F>& front, MatrixNode<F>& X );
template<typename F>
void SolveAfter
( const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front, DistMultiVecNode<F>& X );
template<typename F>
void SolveAfter
( const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front, DistMatrixNode<F>& X );

template<typename F>
Int SolveWithIterativeRefinement
( const SparseMatrix<F>& A,
  const vector<Int>& invMap, const SymmNodeInfo& info,
  const SymmFront<F>& front, Matrix<F>& y,
  Base<F> minReductionFactor=2, Int maxRefineIts=10 );
template<typename F>
Int SolveWithIterativeRefinement
( const DistSparseMatrix<F>& A,
  const DistMap& invMap, const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front, DistMultiVec<F>& y,
  Base<F> minReductionFactor=2, Int maxRefineIts=10 );

// Solve linear system with the implicit representations of L, D, and P
// --------------------------------------------------------------------
template<typename F>
void SolveAfter
( const Matrix<F>& A, const Matrix<F>& dSub,
  const Matrix<Int>& p, Matrix<F>& B, bool conjugated );
template<typename F>
void SolveAfter
( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& dSub,
  const AbstractDistMatrix<Int>& p, AbstractDistMatrix<F>& B, 
  bool conjugated );

} // namespace ldl

// Regularized Quasi-semidefinite LDL^H factorization
// ==================================================
// NOTE: If the pivot candidate is not at least as large as the pivot tolerance
//       and with the implied sign, then it is increased by the specified value.
//       If aPriori is set to true, then this regularization is *always* 
//       applied.
template<typename F>
void RegularizedQSDLDL
( Matrix<F>& A, Base<F> pivTol,
  const Matrix<Base<F>>& regCand,
        Matrix<Base<F>>& reg, 
  bool aPriori=false );
template<typename F>
void RegularizedQSDLDL
( AbstractDistMatrix<F>& A, Base<F> pivTol,
  const AbstractDistMatrix<Base<F>>& regCand,
        AbstractDistMatrix<Base<F>>& reg,
  bool aPriori=false );

template<typename F>
void RegularizedQSDLDL
( const SymmNodeInfo& info, SymmFront<F>& front,
  Base<F> pivTol,
  const MatrixNode<Base<F>>& regCand,
        MatrixNode<Base<F>>& reg,
  bool aPriori=false,
  SymmFrontType newType=LDL_2D );
template<typename F>
void RegularizedQSDLDL
( const DistSymmNodeInfo& info, DistSymmFront<F>& front,
  Base<F> pivTol,
  const DistMultiVecNode<Base<F>>& regCand,
        DistMultiVecNode<Base<F>>& reg,
  bool aPriori=false,
  SymmFrontType newType=LDL_2D );

enum RegQSDRefineAlg
{
  REG_REFINE_FGMRES,
  REG_REFINE_LGMRES,
  REG_REFINE_IR,
  REG_REFINE_IR_MOD
};

namespace reg_qsd_ldl {

template<typename F>
Int RegularizedSolveAfter
( const SparseMatrix<F>& A,   const Matrix<Base<F>>& reg,
  const vector<Int>& invMap,  const SymmNodeInfo& info,
  const SymmFront<F>& front,        Matrix<F>& y,
  Base<F> minReductionFactor,       Int maxRefineIts,
  bool progress );
template<typename F>
Int RegularizedSolveAfter
( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap,             const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,           DistMultiVec<F>& y,
  Base<F> minReductionFactor,              Int maxRefineIts,
  bool progress );

template<typename F>
Int SolveAfter
( const SparseMatrix<F>& A,   const Matrix<Base<F>>& reg,
  const vector<Int>& invMap,  const SymmNodeInfo& info,
  const SymmFront<F>& front,        Matrix<F>& y,
  RegQSDRefineAlg refineAlg,
  Base<F> minReductionFactor,       Int maxRefineIts,
  bool progress );
template<typename F>
Int SolveAfter
( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap,             const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,           DistMultiVec<F>& y,
  RegQSDRefineAlg refineAlg,
  Base<F> minReductionFactor,              Int maxRefineIts,
  bool progress );

} // namespace reg_qsd_ldl

// LU
// ==

// NOTE: This is not yet made use of, but the fully-pivoted version of LU
//       should (soon?) accept it as an argument and potentially return one or
//       more of the permutation matrices as the identity
namespace LUPivotTypeNS {
enum LUPivotType
{
    LU_PARTIAL, 
    LU_FULL,
    LU_ROOK, /* not yet supported */
    LU_WITHOUT_PIVOTING
};
}
using namespace LUPivotTypeNS;

// LU without pivoting
// -------------------
template<typename F>
void LU( Matrix<F>& A );
template<typename F>
void LU( AbstractDistMatrix<F>& A );

// LU with partial pivoting
// ------------------------
template<typename F>
void LU( Matrix<F>& A, Matrix<Int>& p );
template<typename F>
void LU( AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& p );

// LU with full pivoting
// ---------------------
template<typename F>
void LU( Matrix<F>& A, Matrix<Int>& p, Matrix<Int>& q );
template<typename F>
void LU
( AbstractDistMatrix<F>& A, 
  AbstractDistMatrix<Int>& p, AbstractDistMatrix<Int>& q );

// Rank-one modification of a partially-pivoted LU factorization
// -------------------------------------------------------------
template<typename F>
void LUMod
( Matrix<F>& A, Matrix<Int>& p,
  const Matrix<F>& u, const Matrix<F>& v, bool conjugate=true,
  Base<F> tau=0.1 );
template<typename F>
void LUMod
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& p,
  const AbstractDistMatrix<F>& u, const AbstractDistMatrix<F>& v, 
  bool conjugate=true, Base<F> tau=0.1 );

namespace lu {

// Solve linear systems using an implicit unpivoted LU factorization
// -----------------------------------------------------------------
template<typename F>
void SolveAfter( Orientation orientation, const Matrix<F>& A, Matrix<F>& B );
template<typename F>
void SolveAfter
( Orientation orientation, 
  const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B );

// Solve linear systems using an implicit partially-pivoted LU factorization
// -------------------------------------------------------------------------
template<typename F>
void SolveAfter
( Orientation orientation, const Matrix<F>& A,
  const Matrix<Int>& p, Matrix<F>& B );
template<typename F>
void SolveAfter
( Orientation orientation, const AbstractDistMatrix<F>& A,
  const AbstractDistMatrix<Int>& p, AbstractDistMatrix<F>& B );

// Solve linear systems using an implicit fully-pivoted LU factorization
// ---------------------------------------------------------------------
template<typename F>
void SolveAfter
( Orientation orientation, const Matrix<F>& A,
  const Matrix<Int>& p, const Matrix<Int>& q,
        Matrix<F>& B );
template<typename F>
void SolveAfter
( Orientation orientation, const AbstractDistMatrix<F>& A,
  const AbstractDistMatrix<Int>& p, const AbstractDistMatrix<Int>& q,
        AbstractDistMatrix<F>& B );

} // namespace lu

// LQ
// ==

// Overwrite A with both L and the scaled Householder vectors
// ----------------------------------------------------------
template<typename F>
void LQ( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d );
template<typename F>
void LQ
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& t, 
  AbstractDistMatrix<Base<F>>& d );

namespace lq {

// Apply Q using its implicit representation
// -----------------------------------------
template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& t,
  const Matrix<Base<F>>& d, Matrix<F>& B );
template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& t,
  const AbstractDistMatrix<Base<F>>& d, AbstractDistMatrix<F>& B );

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
  const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& t,
  const AbstractDistMatrix<Base<F>>& d, const AbstractDistMatrix<F>& B,
        AbstractDistMatrix<F>& X );

// Overwrite A with L
// ------------------
template<typename F>
void ExplicitTriang( Matrix<F>& A );
template<typename F>
void ExplicitTriang( AbstractDistMatrix<F>& A );

// Overwrite A with Q
// ------------------
template<typename F>
void ExplicitUnitary( Matrix<F>& A );
template<typename F>
void ExplicitUnitary( AbstractDistMatrix<F>& A );

// Return both L and Q such that A = L Q
// -------------------------------------
template<typename F>
void Explicit( Matrix<F>& L, Matrix<F>& A );
template<typename F>
void Explicit( AbstractDistMatrix<F>& L, AbstractDistMatrix<F>& A );

} // namespace lq

// QR factorization
// ================

template<typename Real>
struct QRCtrl
{
    bool colPiv;

    bool boundRank;  
    Int maxRank;

    bool adaptive;
    Real tol;

    bool alwaysRecomputeNorms;

    QRCtrl()
    : colPiv(false), 
      boundRank(false), maxRank(0), adaptive(false), tol(0),
      alwaysRecomputeNorms(false)
    { }
};

// Return an implicit representation of Q and R such that A = Q R
// --------------------------------------------------------------
template<typename F>
void QR
( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d );
template<typename F>
void QR
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& t, 
  AbstractDistMatrix<Base<F>>& d );

// Return an implicit representation of (Q,R,P) such that A P ~= Q R
// -----------------------------------------------------------------
template<typename F>
void QR
( Matrix<F>& A, Matrix<F>& t, 
  Matrix<Base<F>>& d, Matrix<Int>& p,
  const QRCtrl<Base<F>>& ctrl=QRCtrl<Base<F>>() );
template<typename F>
void QR
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& t, 
  AbstractDistMatrix<Base<F>>& d, AbstractDistMatrix<Int>& p,
  const QRCtrl<Base<F>>& ctrl=QRCtrl<Base<F>>() );

namespace qr {

// Apply Q using its implicit representation
// -----------------------------------------
template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& t,
  const Matrix<Base<F>>& d, Matrix<F>& B );
template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& t,
  const AbstractDistMatrix<Base<F>>& d, AbstractDistMatrix<F>& B );

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
  const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& t,
  const AbstractDistMatrix<Base<F>>& d, const AbstractDistMatrix<F>& B,
        AbstractDistMatrix<F>& X );
// TODO: Version which involves permutation matrix

// Cholesky-based QR
// -----------------
template<typename F>
void Cholesky( Matrix<F>& A, Matrix<F>& R );
template<typename F>
void Cholesky( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& R );

// Return R (with non-negative diagonal) such that A = Q R or A P = Q R
// --------------------------------------------------------------------
template<typename F>
void ExplicitTriang
( Matrix<F>& A, const QRCtrl<Base<F>>& ctrl=QRCtrl<Base<F>>() );
template<typename F>
void ExplicitTriang
( AbstractDistMatrix<F>& A, const QRCtrl<Base<F>>& ctrl=QRCtrl<Base<F>>() );

// Return Q such that either A = Q R or A P = Q R
// ----------------------------------------------
template<typename F>
void ExplicitUnitary
( Matrix<F>& A, const QRCtrl<Base<F>>& ctrl=QRCtrl<Base<F>>() );
template<typename F>
void ExplicitUnitary
( AbstractDistMatrix<F>& A, const QRCtrl<Base<F>>& ctrl=QRCtrl<Base<F>>() );

// Return both Q and R such that A = Q R or A P = Q R
// --------------------------------------------------
template<typename F>
void Explicit
( Matrix<F>& A, Matrix<F>& R, 
  const QRCtrl<Base<F>>& ctrl=QRCtrl<Base<F>>() );
template<typename F>
void Explicit
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& R, 
  const QRCtrl<Base<F>>& ctrl=QRCtrl<Base<F>>() );

// Return (Q,R,P) such that A P = Q R
// ----------------------------------
// NOTE: Column-pivoting is performed regardless of the value of ctrl.colPiv 
template<typename F>
void Explicit
( Matrix<F>& A, Matrix<F>& R, 
  Matrix<Int>& P, const QRCtrl<Base<F>>& ctrl=QRCtrl<Base<F>>() );
template<typename F>
void Explicit
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& R, 
  AbstractDistMatrix<Int>& P, const QRCtrl<Base<F>>& ctrl=QRCtrl<Base<F>>() );

template<typename F>
struct TreeData
{
    Matrix<F> QR0, t0;
    Matrix<Base<F>> d0;
    vector<Matrix<F>> QRList;
    vector<Matrix<F>> tList;
    vector<Matrix<Base<F>>> dList;

    TreeData( Int numStages=0 )
    : QRList(numStages), tList(numStages), dList(numStages)
    { }

    TreeData( TreeData<F>&& treeData )
    : QR0(move(treeData.QR0)),
      t0(move(treeData.t0)),
      d0(move(treeData.d0)),
      QRList(move(treeData.QRList)),
      tList(move(treeData.tList)),
      dList(move(treeData.dList))
    { }

    TreeData<F>& operator=( TreeData<F>&& treeData )
    {
        QR0 = move(treeData.QR0);
        t0 = move(treeData.t0);
        d0 = move(treeData.d0);
        QRList = move(treeData.QRList);
        tList = move(treeData.tList);
        dList = move(treeData.dList);
        return *this;
    }
};

// Return an implicit tall-skinny QR factorization
template<typename F>
TreeData<F> TS( const AbstractDistMatrix<F>& A );

// Return an explicit tall-skinny QR factorization
template<typename F>
void ExplicitTS( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& R );

namespace ts {

template<typename F>
Matrix<F>& RootQR
( const AbstractDistMatrix<F>& A, TreeData<F>& treeData );

template<typename F>
const Matrix<F>& RootQR
( const AbstractDistMatrix<F>& A, const TreeData<F>& treeData );

template<typename F>
void Reduce( const AbstractDistMatrix<F>& A, TreeData<F>& treeData );

template<typename F>
void Scatter( AbstractDistMatrix<F>& A, const TreeData<F>& treeData );

} // namespace ts

} // namespace qr

// RQ
// ==
template<typename F>
void RQ( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d );
template<typename F>
void RQ
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& t, 
  AbstractDistMatrix<Base<F>>& d );

namespace rq {

template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& t,
  const Matrix<Base<F>>& d, Matrix<F>& B );
template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& t,
  const AbstractDistMatrix<Base<F>>& d, AbstractDistMatrix<F>& B );

template<typename F>
void SolveAfter
( Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& t,
  const Matrix<Base<F>>& d, const Matrix<F>& B,
        Matrix<F>& X );
template<typename F>
void SolveAfter
( Orientation orientation,
  const AbstractDistMatrix<F      >& A, const AbstractDistMatrix<F>& t,
  const AbstractDistMatrix<Base<F>>& d, const AbstractDistMatrix<F>& B,
        AbstractDistMatrix<F      >& X );

// TODO: Think about ensuring this ordering is consistent with lq::Explicit
template<typename F>
void Cholesky( Matrix<F>& A, Matrix<F>& R );
template<typename F>
void Cholesky( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& R );

template<typename F>
void ExplicitTriang( Matrix<F>& A );
template<typename F>
void ExplicitTriang( AbstractDistMatrix<F>& A );

} // namespace rq

// Generalized QR
// ==============
template<typename F>
void GQR
( Matrix<F>& A, Matrix<F>& tA, Matrix<Base<F>>& dA,
  Matrix<F>& B, Matrix<F>& tB, Matrix<Base<F>>& dB );
template<typename F>
void GQR
( AbstractDistMatrix<F>& A, 
  AbstractDistMatrix<F>& tA, AbstractDistMatrix<Base<F>>& dA,
  AbstractDistMatrix<F>& B, 
  AbstractDistMatrix<F>& tB, AbstractDistMatrix<Base<F>>& dB );

namespace gqr {

template<typename F>
void ExplicitTriang( Matrix<F>& A, Matrix<F>& B );
template<typename F>
void ExplicitTriang( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B );

} // namespace gqr

// Generalized RQ
// ==============
template<typename F>
void GRQ
( Matrix<F>& A, Matrix<F>& tA, Matrix<Base<F>>& dA,
  Matrix<F>& B, Matrix<F>& tB, Matrix<Base<F>>& dB );
template<typename F>
void GRQ
( AbstractDistMatrix<F>& A, 
  AbstractDistMatrix<F>& tA, AbstractDistMatrix<Base<F>>& dA,
  AbstractDistMatrix<F>& B, 
  AbstractDistMatrix<F>& tB, AbstractDistMatrix<Base<F>>& dB );

namespace grq {

template<typename F>
void ExplicitTriang( Matrix<F>& A, Matrix<F>& B );
template<typename F>
void ExplicitTriang( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B );

} // namespace grq

// Interpolative Decomposition
// ===========================
template<typename F>
void ID
( const Matrix<F>& A, Matrix<Int>& p, 
        Matrix<F>& Z, 
  const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>() );
template<typename F>
void ID
( const AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& p,
        AbstractDistMatrix<F>& Z, 
  const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>() );

template<typename F>
void ID
( Matrix<F>& A, Matrix<Int>& p, 
  Matrix<F>& Z, const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>(), 
  bool canOverwrite=false );
template<typename F>
void ID
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& p,
  AbstractDistMatrix<F>& Z, const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>(), 
  bool canOverwrite=false );

// Skeleton
// ========
template<typename F>
void Skeleton
( const Matrix<F>& A,
  Matrix<Int>& pR, Matrix<Int>& pC,
  Matrix<F>& Z, const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>() );
template<typename F>
void Skeleton
( const AbstractDistMatrix<F>& A,
  AbstractDistMatrix<Int>& pR, AbstractDistMatrix<Int>& pC,
  AbstractDistMatrix<F>& Z, const QRCtrl<Base<F>> ctrl=QRCtrl<Base<F>>() );

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
