/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_FACTOR_HPP
#define EL_FACTOR_HPP

#include <El/lapack_like/perm.hpp>
#include <El/lapack_like/util.hpp>
#include <El/lapack_like/factor/ldl/sparse/symbolic.hpp>
#include <El/lapack_like/factor/ldl/sparse/numeric.hpp>

namespace El {

// Cholesky
// ========
template<typename Field>
void Cholesky( UpperOrLower uplo, Matrix<Field>& A );
template<typename Field>
void Cholesky
( UpperOrLower uplo, AbstractDistMatrix<Field>& A, bool scalapack=false );
template<typename Field>
void Cholesky( UpperOrLower uplo, DistMatrix<Field,STAR,STAR>& A );

template<typename Field>
void ReverseCholesky( UpperOrLower uplo, Matrix<Field>& A );
template<typename Field>
void ReverseCholesky( UpperOrLower uplo, AbstractDistMatrix<Field>& A );
template<typename Field>
void ReverseCholesky( UpperOrLower uplo, DistMatrix<Field,STAR,STAR>& A );

template<typename Field>
void Cholesky( UpperOrLower uplo, Matrix<Field>& A, Permutation& P );
template<typename Field>
void Cholesky
( UpperOrLower uplo, AbstractDistMatrix<Field>& A, DistPermutation& P );

template<typename Field>
void CholeskyMod
( UpperOrLower uplo,
  Matrix<Field>& T,
  Base<Field> alpha,
  Matrix<Field>& V );
template<typename Field>
void CholeskyMod
( UpperOrLower uplo,
  AbstractDistMatrix<Field>& T,
  Base<Field> alpha,
  AbstractDistMatrix<Field>& V );

template<typename Field>
void HPSDCholesky( UpperOrLower uplo, Matrix<Field>& A );
template<typename Field>
void HPSDCholesky( UpperOrLower uplo, AbstractDistMatrix<Field>& A );

namespace cholesky {

template<typename Field>
void SolveAfter
( UpperOrLower uplo,
  Orientation orientation,
  const Matrix<Field>& A,
        Matrix<Field>& B );
template<typename Field>
void SolveAfter
( UpperOrLower uplo,
  Orientation orientation,
  const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B );

template<typename Field>
void SolveAfter
( UpperOrLower uplo,
  Orientation orientation,
  const Matrix<Field>& A,
  const Permutation& P,
        Matrix<Field>& B );
template<typename Field>
void SolveAfter
( UpperOrLower uplo,
  Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const DistPermutation& P,
        AbstractDistMatrix<Field>& B );

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
    /* TODO(poulson): Diagonal pivoting? */
};
}
using namespace LDLPivotTypeNS;

template<typename Real>
Real LDLPivotConstant( LDLPivotType pivType )
{
    // TODO(poulson): Check that the Bunch-Parlett choice is optimal
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
template<typename Field>
void LDL( Matrix<Field>& A, bool conjugate );
template<typename Field>
void LDL( AbstractDistMatrix<Field>& A, bool conjugate );
template<typename Field>
void LDL( DistMatrix<Field,STAR,STAR>& A, bool conjugate );

// Return an implicit representation of a pivoted LDL factorization of A
// ---------------------------------------------------------------------
template<typename Field>
void LDL
( Matrix<Field>& A,
  Matrix<Field>& dSub,
  Permutation& P,
  bool conjugate,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );
template<typename Field>
void LDL
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& dSub,
  DistPermutation& P,
  bool conjugate,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );

namespace ldl {

// Compute the inertia triplet of a Hermitian matrix's LDL^H factorization
// -----------------------------------------------------------------------
template<typename Field>
InertiaType Inertia
( const Matrix<Base<Field>>& d,
  const Matrix<Field>& dSub );
template<typename Field>
InertiaType Inertia
( const AbstractDistMatrix<Base<Field>>& d,
  const AbstractDistMatrix<Field>& dSub );

// Multiply vectors using an implicit representation of an LDL factorization
// -------------------------------------------------------------------------
template<typename Field>
void MultiplyAfter
( const Matrix<Field>& A,
        Matrix<Field>& B,
  bool conjugated );
template<typename Field>
void MultiplyAfter
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B,
  bool conjugated );

// Multiply vectors using an implicit representation of a pivoted LDL fact.
// ------------------------------------------------------------------------
template<typename Field>
void MultiplyAfter
( const Matrix<Field>& A,
  const Matrix<Field>& dSub,
  const Permutation& P,
        Matrix<Field>& B,
  bool conjugated );
template<typename Field>
void MultiplyAfter
( const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& dSub,
  const DistPermutation& P,
        AbstractDistMatrix<Field>& B,
  bool conjugated );

// Solve linear systems using an implicit LDL factorization
// --------------------------------------------------------
template<typename Field>
void SolveAfter
( const Matrix<Field>& A,
        Matrix<Field>& B,
  bool conjugated );
template<typename Field>
void SolveAfter
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B,
  bool conjugated );

// Solve linear system with the implicit representations of L, D, and P
// --------------------------------------------------------------------
template<typename Field>
void SolveAfter
( const Matrix<Field>& A,
  const Matrix<Field>& dSub,
  const Permutation& P,
        Matrix<Field>& B,
  bool conjugated );
template<typename Field>
void SolveAfter
( const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& dSub,
  const DistPermutation& P,
        AbstractDistMatrix<Field>& B,
  bool conjugated );

} // namespace ldl

// Solve a linear system with a regularized factorization
// ======================================================
enum RegSolveAlg
{
  REG_SOLVE_FGMRES,
  REG_SOLVE_LGMRES
};

template<typename Real>
struct RegSolveCtrl
{
    RegSolveAlg alg=REG_SOLVE_FGMRES;
    Real relTol;
    Real relTolRefine;
    Int maxIts=4;
    Int maxRefineIts=2;
    Int restart=4;
    bool progress=false;
    bool time=false;

    RegSolveCtrl()
    {
        const Real eps = limits::Epsilon<Real>();
        relTol = Pow(eps,Real(0.5));
        relTolRefine = Pow(eps,Real(0.8));
    }
};

namespace reg_ldl {

template<typename Field>
Int RegularizedSolveAfter
( const SparseMatrix<Field>& A,
  const Matrix<Base<Field>>& reg,
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
        Base<Field> relTolRefine,
        Int maxRefineIts,
        bool progress=false,
        bool time=false );
template<typename Field>
Int RegularizedSolveAfter
( const DistSparseMatrix<Field>& A,
  const DistMultiVec<Base<Field>>& reg,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
        Base<Field> relTolRefine,
        Int maxRefineIts,
        bool progress=false,
        bool time=false );

template<typename Field>
Int RegularizedSolveAfter
( const SparseMatrix<Field>& A,
  const Matrix<Base<Field>>& reg,
  const Matrix<Base<Field>>& d,
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
        Base<Field> relTolRefine,
        Int maxRefineIts,
        bool progress=false,
        bool time=false );
template<typename Field>
Int RegularizedSolveAfter
( const DistSparseMatrix<Field>& A,
  const DistMultiVec<Base<Field>>& reg,
  const DistMultiVec<Base<Field>>& d,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
        Base<Field> relTolRefine,
        Int maxRefineIts,
        bool progress=false,
        bool time=false );

template<typename Field>
Int SolveAfter
( const SparseMatrix<Field>& A,
  const Matrix<Base<Field>>& reg,
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
  const RegSolveCtrl<Base<Field>>& ctrl );
template<typename Field>
Int SolveAfter
( const DistSparseMatrix<Field>& A,
  const DistMultiVec<Base<Field>>& reg,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
  const RegSolveCtrl<Base<Field>>& ctrl );

template<typename Field>
Int SolveAfter
( const SparseMatrix<Field>& A,
  const Matrix<Base<Field>>& reg,
  const Matrix<Base<Field>>& d,
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
  const RegSolveCtrl<Base<Field>>& ctrl );
template<typename Field>
Int SolveAfter
( const DistSparseMatrix<Field>& A,
  const DistMultiVec<Base<Field>>& reg,
  const DistMultiVec<Base<Field>>& d,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
  const RegSolveCtrl<Base<Field>>& ctrl );

} // namespace reg_ldl

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
template<typename Field>
void LU( Matrix<Field>& A );
template<typename Field>
void LU( AbstractDistMatrix<Field>& A );
template<typename Field>
void LU( DistMatrix<Field,STAR,STAR>& A );

// LU with partial pivoting
// ------------------------
template<typename Field>
void LU( Matrix<Field>& A, Permutation& P );
template<typename Field>
void LU( AbstractDistMatrix<Field>& A, DistPermutation& P );

// LU with full pivoting
// ---------------------
// P A Q^T = L U
template<typename Field>
void LU
( Matrix<Field>& A,
  Permutation& P,
  Permutation& Q );
template<typename Field>
void LU
( AbstractDistMatrix<Field>& A,
  DistPermutation& P,
  DistPermutation& Q );

// Low-rank modification of a partially-pivoted LU factorization
// -------------------------------------------------------------
// NOTE: This routine currently performs a sequence of rank-one updates
// and will eventually be generalized to a (much faster) single-pass
// algorithm.
template<typename Field>
void LUMod
(       Matrix<Field>& A,
        Permutation& P,
  const Matrix<Field>& U,
  const Matrix<Field>& V,
  bool conjugate=true,
  Base<Field> tau=Base<Field>(1)/Base<Field>(10) );
template<typename Field>
void LUMod
(       AbstractDistMatrix<Field>& A,
        DistPermutation& P,
  const AbstractDistMatrix<Field>& U,
  const AbstractDistMatrix<Field>& V,
  bool conjugate=true,
  Base<Field> tau=Base<Field>(1)/Base<Field>(10) );

namespace lu {

// Solve linear systems using an implicit unpivoted LU factorization
// -----------------------------------------------------------------
template<typename Field>
void SolveAfter
( Orientation orientation,
  const Matrix<Field>& A,
        Matrix<Field>& B );
template<typename Field>
void SolveAfter
( Orientation orientation,
  const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B );

// Solve linear systems using an implicit partially-pivoted LU factorization
// -------------------------------------------------------------------------
template<typename Field>
void SolveAfter
( Orientation orientation,
  const Matrix<Field>& A,
  const Permutation& P,
        Matrix<Field>& B );
template<typename Field>
void SolveAfter
( Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const DistPermutation& P,
        AbstractDistMatrix<Field>& B );

// Solve linear systems using an implicit fully-pivoted LU factorization
// ---------------------------------------------------------------------
template<typename Field>
void SolveAfter
( Orientation orientation,
  const Matrix<Field>& A,
  const Permutation& P,
  const Permutation& Q,
        Matrix<Field>& B );
template<typename Field>
void SolveAfter
( Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const DistPermutation& P,
  const DistPermutation& Q,
        AbstractDistMatrix<Field>& B );

} // namespace lu

// LQ
// ==

// Overwrite A with both L and the scaled Householder vectors
// ----------------------------------------------------------
template<typename Field>
void LQ
( Matrix<Field>& A,
  Matrix<Field>& householderScalars,
  Matrix<Base<Field>>& signature );
template<typename Field>
void LQ
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& householderScalars,
  AbstractDistMatrix<Base<Field>>& signature );

namespace lq {

// Apply Q using its implicit representation
// -----------------------------------------
template<typename Field>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const Matrix<Field>& A,
  const Matrix<Field>& householderScalars,
  const Matrix<Base<Field>>& signature,
        Matrix<Field>& B );
template<typename Field>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& householderScalars,
  const AbstractDistMatrix<Base<Field>>& signature,
        AbstractDistMatrix<Field>& B );

// Solve a linear system with the implicit representations of L and Q
// ------------------------------------------------------------------
template<typename Field>
void SolveAfter
( Orientation orientation,
  const Matrix<Field>& A,
  const Matrix<Field>& householderScalars,
  const Matrix<Base<Field>>& signature,
  const Matrix<Field>& B,
        Matrix<Field>& X );
template<typename Field>
void SolveAfter
( Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& householderScalars,
  const AbstractDistMatrix<Base<Field>>& signature,
  const AbstractDistMatrix<Field>& B,
        AbstractDistMatrix<Field>& X );

// Overwrite A with L
// ------------------
template<typename Field>
void ExplicitTriang( Matrix<Field>& A );
template<typename Field>
void ExplicitTriang( AbstractDistMatrix<Field>& A );

// Overwrite A with Q
// ------------------
template<typename Field>
void ExplicitUnitary( Matrix<Field>& A );
template<typename Field>
void ExplicitUnitary( AbstractDistMatrix<Field>& A );

// Return both L and Q such that A = L Q
// -------------------------------------
template<typename Field>
void Explicit( Matrix<Field>& L, Matrix<Field>& A );
template<typename Field>
void Explicit( AbstractDistMatrix<Field>& L, AbstractDistMatrix<Field>& A );

} // namespace lq

// QR factorization
// ================

template<typename Real>
struct QRCtrl
{
    bool colPiv=false;

    bool boundRank=false;
    Int maxRank=0;

    bool adaptive=false;
    Real tol=Real(0);

    bool alwaysRecomputeNorms=false;

    // Selecting for the smallest norm first is an important preprocessing
    // step for LLL suggested by Wubben et al.
    //
    // Ideally a black-box reduction operation could be provided by the user
    // instead, as it is often the case that one may desire a custom pivoting
    // rule.
    bool smallestFirst=false;
};

// Return an implicit representation of Q and R such that A = Q R
// --------------------------------------------------------------
template<typename Field>
void QR
( Matrix<Field>& A,
  Matrix<Field>& householderScalars,
  Matrix<Base<Field>>& signature );
template<typename Field>
void QR
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& householderScalars,
  AbstractDistMatrix<Base<Field>>& signature );

// Return an implicit representation of (Q,R,Omega) such that A Omega^T ~= Q R
// ---------------------------------------------------------------------------
template<typename Field>
void QR
( Matrix<Field>& A,
  Matrix<Field>& householderScalars,
  Matrix<Base<Field>>& signature,
  Permutation& Omega,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>() );
template<typename Field>
void QR
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& householderScalars,
  AbstractDistMatrix<Base<Field>>& signature,
  DistPermutation& Omega,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>() );

namespace qr {

// Apply Q using its implicit representation
// -----------------------------------------
template<typename Field>
void ApplyQ
( LeftOrRight side,
  Orientation orientation,
  const Matrix<Field>& A,
  const Matrix<Field>& householderScalars,
  const Matrix<Base<Field>>& signature,
        Matrix<Field>& B );
template<typename Field>
void ApplyQ
( LeftOrRight side,
  Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& householderScalars,
  const AbstractDistMatrix<Base<Field>>& signature,
        AbstractDistMatrix<Field>& B );

// Solve a linear system with the implicit QR factorization
// --------------------------------------------------------
template<typename Field>
void SolveAfter
( Orientation orientation,
  const Matrix<Field>& A,
  const Matrix<Field>& householderScalars,
  const Matrix<Base<Field>>& signature,
  const Matrix<Field>& B,
        Matrix<Field>& X );
template<typename Field>
void SolveAfter
( Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& householderScalars,
  const AbstractDistMatrix<Base<Field>>& signature,
  const AbstractDistMatrix<Field>& B,
        AbstractDistMatrix<Field>& X );
// TODO(poulson): Version which involves permutation matrix

// Cholesky-based QR
// -----------------
template<typename Field>
void Cholesky( Matrix<Field>& A, Matrix<Field>& R );
template<typename Field>
void Cholesky( AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& R );

// Return R (with non-negative diagonal) such that A = Q R or A Omega^T = Q R
// --------------------------------------------------------------------------
template<typename Field>
void ExplicitTriang
( Matrix<Field>& A,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>() );
template<typename Field>
void ExplicitTriang
( AbstractDistMatrix<Field>& A,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>() );

// Return Q such that either A = Q R or A Omega^T = Q R
// ----------------------------------------------------
template<typename Field>
void ExplicitUnitary
( Matrix<Field>& A,
  bool thinQ=true,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>() );
template<typename Field>
void ExplicitUnitary
( AbstractDistMatrix<Field>& A,
  bool thinQ=true,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>() );

// Return both Q and R such that A = Q R or A Omega^T = Q R
// --------------------------------------------------------
template<typename Field>
void Explicit
( Matrix<Field>& A,
  Matrix<Field>& R,
  bool thinQ=true,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>() );
template<typename Field>
void Explicit
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& R,
  bool thinQ=true,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>() );

// Return (Q,R,Omega) such that A Omega^T = Q R
// --------------------------------------------
// NOTE: Column-pivoting is performed regardless of the value of ctrl.colPiv
template<typename Field>
void Explicit
( Matrix<Field>& A,
  Matrix<Field>& R,
  Matrix<Int>& Omega,
  bool thinQ=true,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>() );
template<typename Field>
void Explicit
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& R,
  AbstractDistMatrix<Int>& Omega,
  bool thinQ=true,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>() );

// Swap neighboring columns (j,j+1) and update the QR factorization
// ----------------------------------------------------------------
template<typename Field>
void NeighborColSwap
( Matrix<Field>& Q,
  Matrix<Field>& R,
  Int j );

// Swap disjoint sets of neighboring columns and update the QR factorization
// -------------------------------------------------------------------------
template<typename Field>
void DisjointNeighborColSwaps
(       Matrix<Field>& Q,
        Matrix<Field>& R,
  const Matrix<Int>& colSwaps );

template<typename Field>
struct TreeData
{
    Matrix<Field> QR0, householderScalars0;
    Matrix<Base<Field>> signature0;
    vector<Matrix<Field>> QRList;
    vector<Matrix<Field>> householderScalarsList;
    vector<Matrix<Base<Field>>> signatureList;

    TreeData( Int numStages=0 )
    : QRList(numStages),
      householderScalarsList(numStages),
      signatureList(numStages)
    { }

    TreeData( TreeData<Field>&& treeData )
    : QR0(move(treeData.QR0)),
      householderScalars0(move(treeData.householderScalars0)),
      signature0(move(treeData.signature0)),
      QRList(move(treeData.QRList)),
      householderScalarsList(move(treeData.householderScalarsList)),
      signatureList(move(treeData.signatureList))
    { }

    TreeData<Field>& operator=( TreeData<Field>&& treeData )
    {
        QR0 = move(treeData.QR0);
        householderScalars0 = move(treeData.householderScalars0);
        signature0 = move(treeData.signature0);
        QRList = move(treeData.QRList);
        householderScalarsList = move(treeData.householderScalarsList);
        signatureList = move(treeData.signatureList);
        return *this;
    }
};

// Return an implicit tall-skinny QR factorization
template<typename Field>
TreeData<Field> TS( const AbstractDistMatrix<Field>& A );

// Return an explicit tall-skinny QR factorization
template<typename Field>
void ExplicitTS( AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& R );

namespace ts {

template<typename Field>
Matrix<Field>& RootQR
( const AbstractDistMatrix<Field>& A, TreeData<Field>& treeData );

template<typename Field>
const Matrix<Field>& RootQR
( const AbstractDistMatrix<Field>& A, const TreeData<Field>& treeData );

template<typename Field>
void Reduce( const AbstractDistMatrix<Field>& A, TreeData<Field>& treeData );

template<typename Field>
void Scatter( AbstractDistMatrix<Field>& A, const TreeData<Field>& treeData );

} // namespace ts

} // namespace qr

// RQ
// ==
template<typename Field>
void RQ
( Matrix<Field>& A,
  Matrix<Field>& householderScalars,
  Matrix<Base<Field>>& signature );
template<typename Field>
void RQ
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& householderScalars,
  AbstractDistMatrix<Base<Field>>& signature );

namespace rq {

template<typename Field>
void ApplyQ
( LeftOrRight side,
  Orientation orientation,
  const Matrix<Field>& A,
  const Matrix<Field>& householderScalars,
  const Matrix<Base<Field>>& signature,
        Matrix<Field>& B );
template<typename Field>
void ApplyQ
( LeftOrRight side,
  Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& householderScalars,
  const AbstractDistMatrix<Base<Field>>& signature,
        AbstractDistMatrix<Field>& B );

template<typename Field>
void SolveAfter
( Orientation orientation,
  const Matrix<Field>& A,
  const Matrix<Field>& householderScalars,
  const Matrix<Base<Field>>& signature,
  const Matrix<Field>& B,
        Matrix<Field>& X );
template<typename Field>
void SolveAfter
( Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& householderScalars,
  const AbstractDistMatrix<Base<Field>>& signature,
  const AbstractDistMatrix<Field>& B,
        AbstractDistMatrix<Field>& X );

// TODO(poulson): Think about ensuring this ordering is consistent with
// lq::Explicit
template<typename Field>
void Cholesky( Matrix<Field>& A, Matrix<Field>& R );
template<typename Field>
void Cholesky( AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& R );

template<typename Field>
void ExplicitTriang( Matrix<Field>& A );
template<typename Field>
void ExplicitTriang( AbstractDistMatrix<Field>& A );

} // namespace rq

// Generalized QR
// ==============
template<typename Field>
void GQR
( Matrix<Field>& A,
  Matrix<Field>& householderScalarsA,
  Matrix<Base<Field>>& signatureA,
  Matrix<Field>& B,
  Matrix<Field>& householderScalarsB,
  Matrix<Base<Field>>& signatureB );
template<typename Field>
void GQR
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& householderScalarsA,
  AbstractDistMatrix<Base<Field>>& signatureA,
  AbstractDistMatrix<Field>& B,
  AbstractDistMatrix<Field>& householderScalarsB,
  AbstractDistMatrix<Base<Field>>& signatureB );

namespace gqr {

template<typename Field>
void ExplicitTriang( Matrix<Field>& A, Matrix<Field>& B );
template<typename Field>
void ExplicitTriang
( AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& B );

} // namespace gqr

// Generalized RQ
// ==============
template<typename Field>
void GRQ
( Matrix<Field>& A,
  Matrix<Field>& householderScalarsA,
  Matrix<Base<Field>>& signatureA,
  Matrix<Field>& B,
  Matrix<Field>& householderScalarsB,
  Matrix<Base<Field>>& signatureB );
template<typename Field>
void GRQ
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& householderScalarsA,
  AbstractDistMatrix<Base<Field>>& signatureA,
  AbstractDistMatrix<Field>& B,
  AbstractDistMatrix<Field>& householderScalarsB,
  AbstractDistMatrix<Base<Field>>& signatureB );

namespace grq {

template<typename Field>
void ExplicitTriang( Matrix<Field>& A, Matrix<Field>& B );
template<typename Field>
void ExplicitTriang
( AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& B );

} // namespace grq

// Interpolative Decomposition
// ===========================
template<typename Field>
void ID
( const Matrix<Field>& A,
        Permutation& P,
        Matrix<Field>& Z,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>() );
template<typename Field>
void ID
( const AbstractDistMatrix<Field>& A,
        DistPermutation& P,
        AbstractDistMatrix<Field>& Z,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>() );

template<typename Field>
void ID
( Matrix<Field>& A,
  Permutation& P,
  Matrix<Field>& Z,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>(),
  bool canOverwrite=false );
template<typename Field>
void ID
( AbstractDistMatrix<Field>& A,
  DistPermutation& P,
  AbstractDistMatrix<Field>& Z,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>(),
  bool canOverwrite=false );

// Skeleton
// ========
template<typename Field>
void Skeleton
( const Matrix<Field>& A,
        Permutation& PR,
        Permutation& PC,
        Matrix<Field>& Z,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>() );
template<typename Field>
void Skeleton
( const AbstractDistMatrix<Field>& A,
        DistPermutation& PR,
        DistPermutation& PC,
        AbstractDistMatrix<Field>& Z,
  const QRCtrl<Base<Field>>& ctrl=QRCtrl<Base<Field>>() );

} // namespace El

#include <El/lapack_like/factor/qr/ProxyHouseholder.hpp>

#endif // ifndef EL_FACTOR_HPP
