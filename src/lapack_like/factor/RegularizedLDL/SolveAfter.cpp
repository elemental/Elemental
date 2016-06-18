/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

namespace reg_ldl {

template<typename F>
inline Int RegularizedSolveAfterNoPromote
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& reg,
  const vector<Int>& invMap, 
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& B,
        Base<F> relTol,
        Int maxRefineIts, 
        bool progress,
        bool time )
{
    DEBUG_CSE

    // TODO: Use time in these lambdas
    auto applyA =
      [&]( const Matrix<F>& X, Matrix<F>& Y )
      {
        Y = X;
        DiagonalScale( LEFT, NORMAL, reg, Y ); 
        Multiply( NORMAL, F(1), A, X, F(1), Y );
      };
    auto applyAInv = 
      [&]( Matrix<F>& Y )
      {
        ldl::MatrixNode<F> YNodal( invMap, info, Y );
        ldl::SolveAfter( info, front, YNodal );
        YNodal.Push( invMap, info, Y );
      };

    return RefinedSolve( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename F>
inline Int RegularizedSolveAfterNoPromote
( const SparseMatrix<F>& A,
  const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d, 
  const vector<Int>& invMap,
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& B,
  Base<F> relTol,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    DEBUG_CSE

    // TODO: Use time in these lambdas
    auto applyA =
      [&]( const Matrix<F>& X, Matrix<F>& Y )
      {
        Y = X;
        DiagonalScale( LEFT, NORMAL, reg, Y ); 
        Multiply( NORMAL, F(1), A, X, F(1), Y );
      };
    auto applyAInv = 
      [&]( Matrix<F>& Y )
      {
        DiagonalSolve( LEFT, NORMAL, d, Y );
        ldl::MatrixNode<F> YNodal( invMap, info, Y );
        ldl::SolveAfter( info, front, YNodal );
        YNodal.Push( invMap, info, Y );
        DiagonalSolve( LEFT, NORMAL, d, Y );
      };

    return RefinedSolve( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename F>
inline DisableIf<IsSame<F,Promote<F>>,Int>
RegularizedSolveAfterPromote
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& reg,
  const vector<Int>& invMap, 
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& B,
  Base<F> relTol,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    DEBUG_CSE
    
    typedef Base<F> Real; 
    typedef Promote<Real> PReal;
    typedef Promote<F> PF;

    // TODO: Avoid reforming these each call
    SparseMatrix<PF> AProm;
    Copy( A, AProm );
    Matrix<PReal> regProm;
    Copy( reg, regProm );

    // TODO: Use time in these lambdas
    auto applyA =
      [&]( const Matrix<PF>& XProm, Matrix<PF>& YProm )
      {
        YProm = XProm; 
        DiagonalScale( LEFT, NORMAL, regProm, YProm ); 
        Multiply( NORMAL, PF(1), AProm, XProm, PF(1), YProm );
      }; 
    auto applyAInv =  
      [&]( Matrix<F>& Y )
      {
        ldl::MatrixNode<F> YNodal( invMap, info, Y );
        ldl::SolveAfter( info, front, YNodal );
        YNodal.Push( invMap, info, Y );
      };

    return PromotedRefinedSolve
           ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename F>
inline EnableIf<IsSame<F,Promote<F>>,Int>
RegularizedSolveAfterPromote
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& reg,
  const vector<Int>& invMap, 
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& B,
  Base<F> relTol,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    DEBUG_CSE
    return RegularizedSolveAfterNoPromote
      ( A, reg, invMap, info, front, B, relTol, maxRefineIts, progress, time );
}

template<typename F>
inline DisableIf<IsSame<F,Promote<F>>,Int>
RegularizedSolveAfterPromote
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d, 
  const vector<Int>& invMap, 
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& B,
  Base<F> relTol,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    DEBUG_CSE
    typedef Base<F> Real;
    typedef Promote<Real> PReal;
    typedef Promote<F> PF;

    // TODO: Avoid reforming these each call
    SparseMatrix<PF> AProm;
    Copy( A, AProm );
    Matrix<PReal> regProm;
    Copy( reg, regProm );

    // TODO: Use time in these lambdas
    auto applyA =
      [&]( const Matrix<PF>& XProm, Matrix<PF>& YProm )
      {
        YProm = XProm; 
        DiagonalScale( LEFT, NORMAL, regProm, YProm ); 
        Multiply( NORMAL, PF(1), AProm, XProm, PF(1), YProm );
      }; 
    auto applyAInv =  
      [&]( Matrix<F>& Y )
      {
        DiagonalSolve( LEFT, NORMAL, d, Y );
        ldl::MatrixNode<F> YNodal( invMap, info, Y );
        ldl::SolveAfter( info, front, YNodal );
        YNodal.Push( invMap, info, Y );
        DiagonalSolve( LEFT, NORMAL, d, Y );
      };

    return PromotedRefinedSolve
           ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename F>
inline EnableIf<IsSame<F,Promote<F>>,Int>
RegularizedSolveAfterPromote
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d, 
  const vector<Int>& invMap, 
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& B,
  Base<F> relTol,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    DEBUG_CSE
    return RegularizedSolveAfter
      ( A, reg, d, invMap, info, front, B,
        relTol, maxRefineIts, progress, time );
}

template<typename F>
Int RegularizedSolveAfter
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& reg,
  const vector<Int>& invMap, 
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& B,
  Base<F> relTol,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    DEBUG_CSE
    return RegularizedSolveAfterPromote
           ( A, reg, invMap, info, front, B, relTol, maxRefineIts, 
             progress, time );
}

template<typename F>
Int RegularizedSolveAfter
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d, 
  const vector<Int>& invMap,
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& B,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE
    return RegularizedSolveAfterPromote
           ( A, reg, d, invMap, info, front, 
             B, relTol, maxRefineIts, progress, time );
}

template<typename F>
inline Int RegularizedSolveAfterNoPromote
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE

    // TODO: Use time in these lambdas
    auto applyA =
      [&]( const DistMultiVec<F>& X, DistMultiVec<F>& Y )
      {
        Y = X;
        DiagonalScale( LEFT, NORMAL, reg, Y ); 
        Multiply( NORMAL, F(1), A, X, F(1), Y );
      };
    auto applyAInv = 
      [&]( DistMultiVec<F>& Y )
      {
        // TODO: Switch to DistMatrixNode with large numbers of RHS
        ldl::DistMultiVecNode<F> YNodal;
        YNodal.Pull( invMap, info, Y, meta );
        ldl::SolveAfter( info, front, YNodal );
        YNodal.Push( invMap, info, Y, meta );
      };

    return RefinedSolve( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename F>
inline Int RegularizedSolveAfterNoPromote
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE
    ldl::DistMultiVecNodeMeta meta;
    return RegularizedSolveAfterNoPromote
           ( A, reg, invMap, info, front, B, meta,
             relTol, maxRefineIts, progress, time );
}

template<typename F>
inline Int RegularizedSolveAfterNoPromote
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE

    // TODO: Use time in these lambdas
    auto applyA =
      [&]( const DistMultiVec<F>& X, DistMultiVec<F>& Y )
      {
        Y = X;
        DiagonalScale( LEFT, NORMAL, reg, Y ); 
        Multiply( NORMAL, F(1), A, X, F(1), Y );
      };
    auto applyAInv = 
      [&]( DistMultiVec<F>& Y )
      {
        // TODO: Switch to DistMatrixNode with large numbers of RHS
        DiagonalSolve( LEFT, NORMAL, d, Y );
        ldl::DistMultiVecNode<F> YNodal;
        YNodal.Pull( invMap, info, Y, meta );
        ldl::SolveAfter( info, front, YNodal );
        YNodal.Push( invMap, info, Y, meta );
        DiagonalSolve( LEFT, NORMAL, d, Y );
      };

    return RefinedSolve( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename F>
inline Int RegularizedSolveAfterNoPromote
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE
    ldl::DistMultiVecNodeMeta meta;
    return RegularizedSolveAfterNoPromote
           ( A, reg, d, invMap, info, front, B, meta,
             relTol, maxRefineIts, progress, time );
}

template<typename F>
inline DisableIf<IsSame<F,Promote<F>>,Int>
RegularizedSolveAfterPromote
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE
    typedef Base<F> Real;
    typedef Promote<Real> PReal;
    typedef Promote<F> PF;

    // TODO: Perform these conversions less frequently at a higher level
    DistSparseMatrix<PF> AProm(A.Comm());
    Copy( A, AProm );
    DistMultiVec<PReal> regProm(reg.Comm());
    Copy( reg, regProm );

    // TODO: Use time in these lambdas
    auto applyA =
      [&]( const DistMultiVec<PF>& XProm, DistMultiVec<PF>& YProm )
      {
        YProm = XProm;
        DiagonalScale( LEFT, NORMAL, regProm, YProm ); 
        Multiply( NORMAL, PF(1), AProm, XProm, PF(1), YProm );
      };
    auto applyAInv = 
      [&]( DistMultiVec<F>& Y )
      {
        // TODO: Switch to DistMatrixNode for large numbers of RHS
        ldl::DistMultiVecNode<F> YNodal;
        YNodal.Pull( invMap, info, Y, meta );
        ldl::SolveAfter( info, front, YNodal );
        YNodal.Push( invMap, info, Y, meta );
      };

    return PromotedRefinedSolve
           ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename F>
inline EnableIf<IsSame<F,Promote<F>>,Int>
RegularizedSolveAfterPromote
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE
    return RegularizedSolveAfterNoPromote
      ( A, reg, invMap, info, front, B, meta,
        relTol, maxRefineIts, progress, time );
}

template<typename F>
inline Int RegularizedSolveAfterPromote
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE
    ldl::DistMultiVecNodeMeta meta;
    return RegularizedSolveAfterPromote
           ( A, reg, invMap, info, front, B, meta,
             relTol, maxRefineIts, progress, time );
}

template<typename F>
inline DisableIf<IsSame<F,Promote<F>>,Int>
RegularizedSolveAfterPromote
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE
    typedef Base<F> Real;
    typedef Promote<Real> PReal;
    typedef Promote<F> PF;

    // TODO: Perform these conversions less frequently at a higher level
    DistSparseMatrix<PF> AProm(A.Comm());
    Copy( A, AProm );
    DistMultiVec<PReal> regProm(reg.Comm());
    Copy( reg, regProm );

    // TODO: Use time in these lambdas
    auto applyA =
      [&]( const DistMultiVec<PF>& XProm, DistMultiVec<PF>& YProm )
      {
        YProm = XProm;
        DiagonalScale( LEFT, NORMAL, regProm, YProm ); 
        Multiply( NORMAL, PF(1), AProm, XProm, PF(1), YProm );
      };
    auto applyAInv = 
      [&]( DistMultiVec<F>& Y )
      {
        DiagonalSolve( LEFT, NORMAL, d, Y );
        ldl::DistMultiVecNode<F> YNodal;
        YNodal.Pull( invMap, info, Y, meta );
        ldl::SolveAfter( info, front, YNodal );
        YNodal.Push( invMap, info, Y, meta );
        DiagonalSolve( LEFT, NORMAL, d, Y );
      };

    return PromotedRefinedSolve
           ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename F>
inline EnableIf<IsSame<F,Promote<F>>,Int>
RegularizedSolveAfterPromote
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE
    return RegularizedSolveAfterNoPromote
      ( A, reg, d, invMap, info, front, B, meta,
        relTol, maxRefineIts, progress, time );
}

template<typename F>
inline Int RegularizedSolveAfterPromote
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d, 
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE
    ldl::DistMultiVecNodeMeta meta;
    return RegularizedSolveAfterPromote
           ( A, reg, d, invMap, info, front, B, meta,
             relTol, maxRefineIts, progress, time );
}

template<typename F>
Int RegularizedSolveAfter
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE
    return RegularizedSolveAfterPromote
    ( A, reg, invMap, info, front, B, meta,
      relTol, maxRefineIts, progress, time );
}

template<typename F>
Int RegularizedSolveAfter
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE
    ldl::DistMultiVecNodeMeta meta;
    return RegularizedSolveAfter
           ( A, reg, invMap, info, front, B, meta,
             relTol, maxRefineIts, progress, time );
}

template<typename F>
Int RegularizedSolveAfter
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d, 
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE
    return RegularizedSolveAfterPromote
    ( A, reg, d, invMap, info, front, B, meta,
      relTol, maxRefineIts, progress, time );
}

template<typename F>
Int RegularizedSolveAfter
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d, 
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_CSE
    ldl::DistMultiVecNodeMeta meta;
    return RegularizedSolveAfter
    ( A, reg, d, invMap, info, front, B, meta,
      relTol, maxRefineIts, progress, time );
}

template<typename F>
Int LGMRESSolveAfter
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& reg,
  const vector<Int>& invMap, 
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& B,
  Base<F> relTol,
  Int restart,
  Int maxIts,
  Base<F> relTolRefine,
  Int maxRefineIts,
  bool progress )
{
    DEBUG_CSE

    auto applyA =
      [&]( F alpha, const Matrix<F>& X, F beta, Matrix<F>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( Matrix<F>& W )
      {
        RegularizedSolveAfter
        ( A, reg, invMap, info, front, W, 
          relTolRefine, maxRefineIts, progress );
      };

    return LGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename F>
Int LGMRESSolveAfter
( const SparseMatrix<F>& A,
  const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d,
  const vector<Int>& invMap,
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& B,
  Base<F> relTol,
  Int restart,
  Int maxIts,
  Base<F> relTolRefine,
  Int maxRefineIts,
  bool progress )
{
    DEBUG_CSE

    auto applyA =
      [&]( F alpha, const Matrix<F>& X, F beta, Matrix<F>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( Matrix<F>& W )
      {
        RegularizedSolveAfter
        ( A, reg, d, invMap, info, front, W, 
          relTolRefine, maxRefineIts, progress );
      };

    return LGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename F>
Int LGMRESSolveAfter
( const DistSparseMatrix<F>& A,
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int restart,
  Int maxIts,
  Base<F> relTolRefine,
  Int maxRefineIts,
  bool progress )
{
    DEBUG_CSE

    auto applyA =
      [&]( F alpha, const DistMultiVec<F>& X, F beta, DistMultiVec<F>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( DistMultiVec<F>& W )
      {
        RegularizedSolveAfter
        ( A, reg, invMap, info, front, W, meta,
          relTolRefine, maxRefineIts, progress );
      };

    return LGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename F>
Int LGMRESSolveAfter
( const DistSparseMatrix<F>& A,
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
  Base<F> relTol,
  Int restart,
  Int maxIts,
  Base<F> relTolRefine,
  Int maxRefineIts,
  bool progress )
{
    DEBUG_CSE
    ldl::DistMultiVecNodeMeta meta;
    return LGMRESSolveAfter
           ( A, reg, invMap, info, front, B, meta,
             relTol, restart, maxIts, relTolRefine, maxRefineIts, progress );
}

template<typename F>
Int LGMRESSolveAfter
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        ldl::DistMultiVecNodeMeta& meta,
        Base<F> relTol,
        Int restart,
        Int maxIts,
        Base<F> relTolRefine,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_CSE

    auto applyA =
      [&]( F alpha, const DistMultiVec<F>& X, F beta, DistMultiVec<F>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( DistMultiVec<F>& W )
      {
        RegularizedSolveAfter
        ( A, reg, d, invMap, info, front, W, meta,
          relTolRefine, maxRefineIts, progress );
      };

    return LGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename F>
Int LGMRESSolveAfter
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
  Base<F> relTol,
  Int restart,
  Int maxIts,
  Base<F> relTolRefine,
  Int maxRefineIts,
  bool progress )
{
    DEBUG_CSE
    ldl::DistMultiVecNodeMeta meta;
    return LGMRESSolveAfter
           ( A, reg, d, invMap, info, front, B, meta,
             relTol, restart, maxIts, relTolRefine, maxRefineIts, progress );
}

template<typename F>
Int FGMRESSolveAfter
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& reg,
  const vector<Int>& invMap, 
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& B,
        Base<F> relTol,
        Int restart,
        Int maxIts,
        Base<F> relTolRefine,
        Int maxRefineIts, 
        bool progress,
        bool time )
{
    DEBUG_CSE

    auto applyA =
      [&]( F alpha, const Matrix<F>& X, F beta, Matrix<F>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( Matrix<F>& W )
      {
        RegularizedSolveAfter
        ( A, reg, invMap, info, front, W, 
          relTolRefine, maxRefineIts, progress );
      };

    return FGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename F>
Int FGMRESSolveAfter
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d,
  const vector<Int>& invMap, 
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& B,
        Base<F> relTol,
        Int restart,
        Int maxIts,
        Base<F> relTolRefine,
        Int maxRefineIts, 
        bool progress,
        bool time )
{
    DEBUG_CSE

    auto applyA =
      [&]( F alpha, const Matrix<F>& X, F beta, Matrix<F>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( Matrix<F>& W )
      {
        RegularizedSolveAfter
        ( A, reg, d, invMap, info, front, W,
          relTolRefine, maxRefineIts, progress );
      };

    return FGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename F>
Int FGMRESSolveAfter
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        ldl::DistMultiVecNodeMeta& meta,
        Base<F> relTol,
        Int restart,
        Int maxIts,
        Base<F> relTolRefine,
        Int maxRefineIts, 
        bool progress,
        bool time )
{
    DEBUG_CSE

    auto applyA =
      [&]( F alpha, const DistMultiVec<F>& X, F beta, DistMultiVec<F>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( DistMultiVec<F>& W )
      {
        RegularizedSolveAfter
        ( A, reg, invMap, info, front, W, meta,
          relTolRefine, maxRefineIts, progress );
      };

    return FGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename F>
Int FGMRESSolveAfter
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        Base<F> relTol,
        Int restart,
        Int maxIts,
        Base<F> relTolRefine,
        Int maxRefineIts, 
        bool progress,
        bool time )
{
    DEBUG_CSE
    ldl::DistMultiVecNodeMeta meta;
    return FGMRESSolveAfter
           ( A, reg, invMap, info, front, B, meta,
             relTol, restart, maxIts, relTolRefine, maxRefineIts,
             progress, time );
}

template<typename F>
Int FGMRESSolveAfter
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        ldl::DistMultiVecNodeMeta& meta,
        Base<F> relTol,
        Int restart,
        Int maxIts,
        Base<F> relTolRefine,
        Int maxRefineIts, 
        bool progress,
        bool time )
{
    DEBUG_CSE

    auto applyA =
      [&]( F alpha, const DistMultiVec<F>& X, F beta, DistMultiVec<F>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( DistMultiVec<F>& W )
      {
        RegularizedSolveAfter
        ( A, reg, d, invMap, info, front, W, meta,
          relTolRefine, maxRefineIts, progress );
      };

    return FGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename F>
Int FGMRESSolveAfter
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        Base<F> relTol,
        Int restart,
        Int maxIts,
        Base<F> relTolRefine,
        Int maxRefineIts, 
        bool progress,
        bool time )
{
    DEBUG_CSE
    ldl::DistMultiVecNodeMeta meta;
    return FGMRESSolveAfter
           ( A, reg, d, invMap, info, front, B, meta,
             relTol, restart, maxIts, relTolRefine, maxRefineIts,
             progress, time );
}

// TODO: Add RGMRES

template<typename F>
Int SolveAfter
( const SparseMatrix<F>& A,
  const Matrix<Base<F>>& reg,
  const vector<Int>& invMap,
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& B,
  const RegSolveCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    switch( ctrl.alg )
    {
    case REG_SOLVE_FGMRES:
        return FGMRESSolveAfter
        ( A, reg, invMap, info, front, B, 
          ctrl.relTol, ctrl.restart, ctrl.maxIts,
          ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress, ctrl.time );
    case REG_SOLVE_LGMRES:
        return LGMRESSolveAfter
        ( A, reg, invMap, info, front, B, 
          ctrl.relTol, ctrl.restart, ctrl.maxIts,
          ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress );
    default:
        LogicError("Invalid refinement algorithm");
        return -1;
    }
}

template<typename F>
Int SolveAfter
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d,
  const vector<Int>& invMap, 
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& B,
  const RegSolveCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    switch( ctrl.alg )
    {
    case REG_SOLVE_FGMRES:
        return FGMRESSolveAfter
        ( A, reg, d, invMap, info, front, B, 
          ctrl.relTol, ctrl.restart, ctrl.maxIts,
          ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress, ctrl.time );
    case REG_SOLVE_LGMRES:
        return LGMRESSolveAfter
        ( A, reg, d, invMap, info, front, B, 
          ctrl.relTol, ctrl.restart, ctrl.maxIts,
          ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress );
    default:
        LogicError("Invalid refinement algorithm");
        return -1;
    }
}

template<typename F>
Int SolveAfter
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        ldl::DistMultiVecNodeMeta& meta,
  const RegSolveCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    switch( ctrl.alg )
    {
    case REG_SOLVE_FGMRES:
        return FGMRESSolveAfter
        ( A, reg, invMap, info, front, B, meta,
          ctrl.relTol, ctrl.restart, ctrl.maxIts,
          ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress, ctrl.time );
    case REG_SOLVE_LGMRES:
        return LGMRESSolveAfter
        ( A, reg, invMap, info, front, B, meta,
          ctrl.relTol, ctrl.restart, ctrl.maxIts,
          ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress );
    default:
        LogicError("Invalid refinement algorithm");
        return -1;
    }
}

template<typename F>
Int SolveAfter
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
  const RegSolveCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    ldl::DistMultiVecNodeMeta meta;
    return SolveAfter( A, reg, invMap, info, front, B, meta, ctrl );
}

template<typename F>
Int SolveAfter
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
        ldl::DistMultiVecNodeMeta& meta,
  const RegSolveCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    switch( ctrl.alg )
    {
    case REG_SOLVE_FGMRES:
        return FGMRESSolveAfter
        ( A, reg, d, invMap, info, front, B, meta,
          ctrl.relTol, ctrl.restart, ctrl.maxIts,
          ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress, ctrl.time );
    case REG_SOLVE_LGMRES:
        return LGMRESSolveAfter
        ( A, reg, d, invMap, info, front, B, meta,
          ctrl.relTol, ctrl.restart, ctrl.maxIts,
          ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress );
    default:
        LogicError("Invalid refinement algorithm");
        return -1;
    }
}

template<typename F>
Int SolveAfter
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& B,
  const RegSolveCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    ldl::DistMultiVecNodeMeta meta;
    return SolveAfter( A, reg, d, invMap, info, front, B, meta, ctrl );
}

#define PROTO(F) \
  template Int RegularizedSolveAfter \
  ( const SparseMatrix<F>& A, \
    const Matrix<Base<F>>& reg, \
    const vector<Int>& invMap, \
    const ldl::NodeInfo& info, \
    const ldl::Front<F>& front, \
          Matrix<F>& B, \
    Base<F> relTol, Int maxRefineIts, bool progress, bool time ); \
  template Int RegularizedSolveAfter \
  ( const SparseMatrix<F>& A, \
    const Matrix<Base<F>>& reg, \
    const Matrix<Base<F>>& d, \
    const vector<Int>& invMap, \
    const ldl::NodeInfo& info, \
    const ldl::Front<F>& front, \
          Matrix<F>& B, \
    Base<F> relTol, Int maxRefineIts, bool progress, bool time ); \
  template Int RegularizedSolveAfter \
  ( const DistSparseMatrix<F>& A, \
    const DistMultiVec<Base<F>>& reg, \
    const DistMap& invMap, \
    const ldl::DistNodeInfo& info, \
    const ldl::DistFront<F>& front, \
          DistMultiVec<F>& B, \
    Base<F> relTol, Int maxRefineIts, bool progress, bool time ); \
  template Int RegularizedSolveAfter \
  ( const DistSparseMatrix<F>& A, \
    const DistMultiVec<Base<F>>& reg, \
    const DistMap& invMap, \
    const ldl::DistNodeInfo& info, \
    const ldl::DistFront<F>& front, \
          DistMultiVec<F>& B, \
          ldl::DistMultiVecNodeMeta& meta, \
    Base<F> relTol, Int maxRefineIts, bool progress, bool time ); \
  template Int RegularizedSolveAfter \
  ( const DistSparseMatrix<F>& A, \
    const DistMultiVec<Base<F>>& reg, \
    const DistMultiVec<Base<F>>& d, \
    const DistMap& invMap, \
    const ldl::DistNodeInfo& info, \
    const ldl::DistFront<F>& front, \
          DistMultiVec<F>& B, \
    Base<F> relTol, Int maxRefineIts, bool progress, bool time ); \
  template Int RegularizedSolveAfter \
  ( const DistSparseMatrix<F>& A, \
    const DistMultiVec<Base<F>>& reg, \
    const DistMultiVec<Base<F>>& d, \
    const DistMap& invMap, \
    const ldl::DistNodeInfo& info, \
    const ldl::DistFront<F>& front, \
          DistMultiVec<F>& B, \
          ldl::DistMultiVecNodeMeta& meta, \
    Base<F> relTol, Int maxRefineIts, bool progress, bool time ); \
  template Int SolveAfter \
  ( const SparseMatrix<F>& A, \
    const Matrix<Base<F>>& reg, \
    const vector<Int>& invMap, \
    const ldl::NodeInfo& info, \
    const ldl::Front<F>& front, \
          Matrix<F>& B, \
    const RegSolveCtrl<Base<F>>& ctrl ); \
  template Int SolveAfter \
  ( const SparseMatrix<F>& A, \
    const Matrix<Base<F>>& reg, \
    const Matrix<Base<F>>& d, \
    const vector<Int>& invMap, \
    const ldl::NodeInfo& info, \
    const ldl::Front<F>& front, \
          Matrix<F>& B, \
    const RegSolveCtrl<Base<F>>& ctrl ); \
  template Int SolveAfter \
  ( const DistSparseMatrix<F>& A, \
    const DistMultiVec<Base<F>>& reg, \
    const DistMap& invMap, \
    const ldl::DistNodeInfo& info, \
    const ldl::DistFront<F>& front, \
          DistMultiVec<F>& B, \
    const RegSolveCtrl<Base<F>>& ctrl ); \
  template Int SolveAfter \
  ( const DistSparseMatrix<F>& A, \
    const DistMultiVec<Base<F>>& reg, \
    const DistMap& invMap, \
    const ldl::DistNodeInfo& info, \
    const ldl::DistFront<F>& front, \
          DistMultiVec<F>& B, \
          ldl::DistMultiVecNodeMeta& meta, \
    const RegSolveCtrl<Base<F>>& ctrl ); \
  template Int SolveAfter \
  ( const DistSparseMatrix<F>& A, \
    const DistMultiVec<Base<F>>& reg, \
    const DistMultiVec<Base<F>>& d, \
    const DistMap& invMap, \
    const ldl::DistNodeInfo& info, \
    const ldl::DistFront<F>& front, \
          DistMultiVec<F>& B, \
    const RegSolveCtrl<Base<F>>& ctrl ); \
  template Int SolveAfter \
  ( const DistSparseMatrix<F>& A, \
    const DistMultiVec<Base<F>>& reg, \
    const DistMultiVec<Base<F>>& d, \
    const DistMap& invMap, \
    const ldl::DistNodeInfo& info, \
    const ldl::DistFront<F>& front, \
          DistMultiVec<F>& B, \
          ldl::DistMultiVecNodeMeta& meta, \
    const RegSolveCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace reg_ldl
} // namespace El
