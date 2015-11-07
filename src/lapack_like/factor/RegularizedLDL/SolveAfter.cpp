/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int IterativeRefinement
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<F>& b,
        Base<F> relTol,
        Int maxRefineIts, 
        bool progress )
{
    DEBUG_ONLY(
      CSE cse("IterativeRefinement");
      if( b.Width() != 1 )    
          LogicError("Expected a single right-hand side");
    )
    auto bOrig = b;
    const Base<F> bNorm = MaxNorm( b );

    // Compute the initial guess
    // =========================
    auto x = b;
    applyAInv( x );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        Matrix<F> dx, xCand, y;
        applyA( x, y );
        b -= y;
        Base<F> errorNorm = MaxNorm( b );
        if( progress )
            Output("original rel error: ",errorNorm/bNorm);
 
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress )
                    Output(errorNorm/bNorm," <= ",relTol);
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            dx = b;
            applyAInv( dx );
            xCand = x;
            xCand += dx;

            // Check the new residual
            // ----------------------
            b = bOrig;
            applyA( xCand, y );
            b -= y;
            auto newErrorNorm = MaxNorm( b );
            if( progress )
                Output("refined rel error: ",newErrorNorm/bNorm);

            if( newErrorNorm < errorNorm )
                x = xCand;
            else
                break;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                break;
        }
    }
    // Store the final result
    // ======================
    b = x;
    return refineIt;
}

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int PromotedIterativeRefinement
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<F>& b,
        Base<F> relTol,
        Int maxRefineIts, 
        bool progress )
{
    DEBUG_ONLY(
      CSE cse("PromotedIterativeRefinement");
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    typedef Base<F> Real;
    typedef Promote<Real> PReal;
    typedef Promote<F> PF;

    Matrix<PF> bProm, bOrigProm;
    Copy( b, bProm );
    Copy( b, bOrigProm );
    const PReal bNorm = MaxNorm( bOrigProm );

    // Compute the initial guess
    // =========================
    applyAInv( b );
    Matrix<PF> xProm;
    Copy( b, xProm );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        Matrix<PF> dxProm, xCandProm, yProm;
        applyA( xProm, yProm );
        bProm -= yProm;
        auto errorNorm = MaxNorm( bProm );
        if( progress )
            Output("original rel error: ",errorNorm/bNorm);
 
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress )
                    Output(errorNorm/bNorm," <= ",relTol);
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            Copy( bProm, b );
            applyAInv( b );
            Copy( b, dxProm );
            xCandProm = xProm;
            xCandProm += dxProm;

            // Check the new residual
            // ----------------------
            applyA( xCandProm, yProm );
            bProm = bOrigProm;
            bProm -= yProm;
            auto newErrorNorm = MaxNorm( bProm );
            if( progress )
                Output("refined rel error: ",newErrorNorm/bNorm);

            if( newErrorNorm < errorNorm )
                xProm = xCandProm;
            else
                break;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                break;
        }
    }
    // Store the final result
    // ======================
    Copy( xProm, b );
    return refineIt;
}

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int IterativeRefinement
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<F>& b,
        Base<F> relTol,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_ONLY(
      CSE cse("IterativeRefinement");
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    mpi::Comm comm = b.Comm();
    const int commRank = mpi::Rank(comm);

    auto bOrig = b;
    const Base<F> bNorm = MaxNorm( b );

    // Compute the initial guess
    // =========================
    auto x = b;
    applyAInv( x );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        DistMultiVec<F> dx(comm), xCand(comm), y(comm);
        applyA( x, y );
        b -= y;
        Base<F> errorNorm = MaxNorm( b );
        if( progress && commRank == 0 )
            Output("original rel error: ",errorNorm/bNorm);

        const Int indent = PushIndent();
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress && commRank == 0 )
                    Output(errorNorm/bNorm," <= ",relTol);
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            dx = b;
            applyAInv( dx );
            xCand = x;
            xCand += dx;

            // Compute the new residual
            // ------------------------
            b = bOrig;
            applyA( xCand, y );
            b -= y;
            Base<F> newErrorNorm = MaxNorm( b );
            if( progress && commRank == 0 )
                Output("refined rel error: ",newErrorNorm/bNorm);

            if( newErrorNorm < errorNorm )
                x = xCand;
            else
                break;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                break;
        }
        SetIndent( indent );
    }
    b = x;
    return refineIt;
}

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int PromotedIterativeRefinement
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<F>& b,
        Base<F> relTol,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_ONLY(
      CSE cse("PromotedIterativeRefinement");
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    typedef Base<F> Real;
    typedef Promote<Real> PReal;
    typedef Promote<F> PF;
    mpi::Comm comm = b.Comm();
    const int commRank = mpi::Rank(comm);

    DistMultiVec<PF> bProm(comm), bOrigProm(comm);
    Copy( b, bProm ); 
    Copy( b, bOrigProm );
    const auto bNorm = Nrm2( bProm );

    // Compute the initial guess
    // =========================
    applyAInv( b );
    DistMultiVec<PF> xProm(comm);
    Copy( b, xProm );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        DistMultiVec<PF> dxProm(comm), xCandProm(comm), yProm(comm);
        applyA( xProm, yProm );
        bProm -= yProm;
        auto errorNorm = Nrm2( bProm );
        if( progress && commRank == 0 )
            Output("original rel error: ",errorNorm/bNorm);

        const Int indent = PushIndent();
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress && commRank == 0 )
                    Output(errorNorm/bNorm," <= ",relTol);
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            Copy( bProm, b );
            applyAInv( b );
            Copy( b, dxProm );
            xCandProm = xProm;
            xCandProm += dxProm;

            // Check the new residual
            // ----------------------
            applyA( xCandProm, yProm );
            bProm = bOrigProm;
            bProm -= yProm;
            auto newErrorNorm = Nrm2( bProm );
            if( progress && commRank == 0 )
                Output("refined rel error: ",newErrorNorm/bNorm);

            if( newErrorNorm < errorNorm )
                xProm = xCandProm;
            else
                break;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                break;
        }
        SetIndent( indent );
    }
    Copy( xProm, b );
    return refineIt;
}

namespace reg_ldl {

// TODO: Switch to returning the relative residual of the refined solution

// TODO: Do not accept iterative refinements which increase the residual norm

// TODO: Implement multi-RHS version of LGMRES

// TODO: Implement multi-RHS version of FGMRES

// TODO: Start with batch version of RegularizedSolveAfter

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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfterNoPromote"))

    // TODO: Use time in these lambdas
    auto applyA =
      [&]( const Matrix<F>& x, Matrix<F>& y )
      {
        y = x;
        DiagonalScale( LEFT, NORMAL, reg, y ); 
        Multiply( NORMAL, F(1), A, x, F(1), y );
      };
    auto applyAInv = 
      [&]( Matrix<F>& y )
      {
        ldl::MatrixNode<F> yNodal( invMap, info, y );
        ldl::SolveAfter( info, front, yNodal );
        yNodal.Push( invMap, info, y );
      };

    // TODO: Perform these in a batch instead
    Int mostRefineIts = 0;
    const Int width = B.Width();
    for( Int j=0; j<width; ++j )
    {
        auto b = B( ALL, IR(j) );
        const Int refineIts = 
          IterativeRefinement
          ( applyA, applyAInv, b, relTol, maxRefineIts, progress );
        mostRefineIts = Max(mostRefineIts,refineIts);
    }
    return mostRefineIts;
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfterNoPromote"))

    // TODO: Use time in these lambdas
    auto applyA =
      [&]( const Matrix<F>& x, Matrix<F>& y )
      {
        y = x;
        DiagonalScale( LEFT, NORMAL, reg, y ); 
        Multiply( NORMAL, F(1), A, x, F(1), y );
      };
    auto applyAInv = 
      [&]( Matrix<F>& y )
      {
        DiagonalSolve( LEFT, NORMAL, d, y );
        ldl::MatrixNode<F> yNodal( invMap, info, y );
        ldl::SolveAfter( info, front, yNodal );
        yNodal.Push( invMap, info, y );
        DiagonalSolve( LEFT, NORMAL, d, y );
      };

    // TODO: Perform these in a batch instead
    Int mostRefineIts = 0;
    const Int width = B.Width();
    for( Int j=0; j<width; ++j )
    {
        auto b = B( ALL, IR(j) );
        const Int refineIts =
          IterativeRefinement
          ( applyA, applyAInv, b, relTol, maxRefineIts, progress );
        mostRefineIts = Max(mostRefineIts,refineIts);
    }
    return mostRefineIts;
}

template<typename F>
inline Int RegularizedSolveAfterPromote
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfterPromote"))
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
      [&]( const Matrix<PF>& xProm, Matrix<PF>& yProm )
      {
        yProm = xProm; 
        DiagonalScale( LEFT, NORMAL, regProm, yProm ); 
        Multiply( NORMAL, PF(1), AProm, xProm, PF(1), yProm );
      }; 
    auto applyAInv =  
      [&]( Matrix<F>& y )
      {
        ldl::MatrixNode<F> yNodal( invMap, info, y );
        ldl::SolveAfter( info, front, yNodal );
        yNodal.Push( invMap, info, y );
      };

    // TODO: Perform these in a batch instead
    Int mostRefineIts = 0;
    const Int width = B.Width();
    for( Int j=0; j<width; ++j )
    {
        auto b = B( ALL, IR(j) );
        const Int refineIts =
          PrommotedIterativeRefinement
          ( applyA, applyAInv, b, relTol, maxRefineIts, progress );
        mostRefineIts = Max(mostRefineIts,refineIts);
    }
    return mostRefineIts;
}

template<typename F>
inline Int RegularizedSolveAfterPromote
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfterPromote"))
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
      [&]( const Matrix<PF>& xProm, Matrix<PF>& yProm )
      {
        yProm = xProm; 
        DiagonalScale( LEFT, NORMAL, regProm, yProm ); 
        Multiply( NORMAL, PF(1), AProm, xProm, PF(1), yProm );
      }; 
    auto applyAInv =  
      [&]( Matrix<F>& y )
      {
        DiagonalSolve( LEFT, NORMAL, d, y );
        ldl::MatrixNode<F> yNodal( invMap, info, y );
        ldl::SolveAfter( info, front, yNodal );
        yNodal.Push( invMap, info, y );
        DiagonalSolve( LEFT, NORMAL, d, y );
      };

    // TODO: Perform these in a batch instead
    Int mostRefineIts = 0;
    const Int width = B.Width();
    for( Int j=0; j<width; ++j )
    {
        auto b = B( ALL, IR(j) );
        const Int refineIts =
          PrommotedIterativeRefinement
          ( applyA, applyAInv, b, relTol, maxRefineIts, progress );
        mostRefineIts = Max(mostRefineIts,refineIts);
    }
    return mostRefineIts;
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfter"))
#ifdef EL_HAVE_QUAD
    return RegularizedSolveAfterPromote
           ( A, reg, invMap, info, front, B, relTol, maxRefineIts, 
             progress, time );
#else
    return RegularizedSolveAfterNoPromote
           ( A, reg, invMap, info, front, B, relTol, maxRefineIts, 
             progress, time );
#endif
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfter"))
#ifdef EL_HAVE_QUAD
    return RegularizedSolveAfterPromote
           ( A, reg, d, invMap, info, front, 
             B, relTol, maxRefineIts, progress, time );
#else
    return RegularizedSolveAfterNoPromote
           ( A, reg, d, invMap, info, front, 
             B, relTol, maxRefineIts, progress, time );
#endif
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfterNoPromote"))

    // TODO: Use time in these lambdas
    auto applyA =
      [&]( const DistMultiVec<F>& x, DistMultiVec<F>& y )
      {
        y = x;
        DiagonalScale( LEFT, NORMAL, reg, y ); 
        Multiply( NORMAL, F(1), A, x, F(1), y );
      };
    auto applyAInv = 
      [&]( DistMultiVec<F>& y )
      {
        ldl::DistMultiVecNode<F> yNodal;
        yNodal.Pull( invMap, info, y, meta );
        ldl::SolveAfter( info, front, yNodal );
        yNodal.Push( invMap, info, y, meta );
      };

    // TODO: Perform these in a batch instead
    const Int height = B.Height();
    const Int width = B.Width();
    Int mostRefineIts = 0;
    DistMultiVec<F> u(B.Comm());
    Zeros( u, height, 1 );
    auto& BLoc = B.Matrix();
    auto& uLoc = u.Matrix();
    for( Int j=0; j<width; ++j )
    {
        auto bLoc = BLoc( ALL, IR(j) );
        Copy( bLoc, uLoc );
        const Int refineIts =
          IterativeRefinement
          ( applyA, applyAInv, u, relTol, maxRefineIts, progress );
        Copy( uLoc, bLoc );
        mostRefineIts = Max(mostRefineIts,refineIts);
    }

    return mostRefineIts;
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfterNoPromote"))
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfterNoPromote"))

    // TODO: Use time in these lambdas
    auto applyA =
      [&]( const DistMultiVec<F>& x, DistMultiVec<F>& y )
      {
        y = x;
        DiagonalScale( LEFT, NORMAL, reg, y ); 
        Multiply( NORMAL, F(1), A, x, F(1), y );
      };
    auto applyAInv = 
      [&]( DistMultiVec<F>& y )
      {
        DiagonalSolve( LEFT, NORMAL, d, y );
        ldl::DistMultiVecNode<F> yNodal;
        yNodal.Pull( invMap, info, y, meta );
        ldl::SolveAfter( info, front, yNodal );
        yNodal.Push( invMap, info, y, meta );
        DiagonalSolve( LEFT, NORMAL, d, y );
      };

    // TODO: Perform these in a batch instead
    const Int height = B.Height();
    const Int width = B.Width();
    Int mostRefineIts = 0;
    DistMultiVec<F> u(B.Comm());
    Zeros( u, height, 1 );
    auto& BLoc = B.Matrix();
    auto& uLoc = u.Matrix();
    for( Int j=0; j<width; ++j )
    {
        auto bLoc = BLoc( ALL, IR(j) );
        Copy( bLoc, uLoc );
        const Int refineIts =
          IterativeRefinement
          ( applyA, applyAInv, u, relTol, maxRefineIts, progress );
        Copy( uLoc, bLoc );
        mostRefineIts = Max(mostRefineIts,refineIts);
    }

    return mostRefineIts;
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfterNoPromote"))
    ldl::DistMultiVecNodeMeta meta;
    return RegularizedSolveAfterNoPromote
           ( A, reg, d, invMap, info, front, B, meta,
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
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfterPromote"))
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
      [&]( const DistMultiVec<PF>& xProm, DistMultiVec<PF>& yProm )
      {
        yProm = xProm;
        DiagonalScale( LEFT, NORMAL, regProm, yProm ); 
        Multiply( NORMAL, PF(1), AProm, xProm, PF(1), yProm );
      };
    auto applyAInv = 
      [&]( DistMultiVec<F>& y )
      {
        ldl::DistMultiVecNode<F> yNodal( invMap, info, y, meta );
        ldl::SolveAfter( info, front, yNodal );
        yNodal.Push( invMap, info, y, meta );
      };

    // TODO: Perform these in a batch instead
    const Int height = B.Height();
    const Int width = B.Width();
    Int mostRefineIts = 0;
    DistMultiVec<F> u(B.Comm());
    Zeros( u, height, 1 );
    auto& BLoc = B.Matrix();
    auto& uLoc = u.Matrix();
    for( Int j=0; j<width; ++j )
    {
        auto bLoc = BLoc( ALL, IR(j) );
        Copy( bLoc, uLoc );
        const Int refineIts =
          PromotedIterativeRefinement
          ( applyA, applyAInv, u, relTol, maxRefineIts, progress );
        Copy( uLoc, bLoc );
        mostRefineIts = Max(mostRefineIts,refineIts);
    }

    return mostRefineIts;
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfterPromote"))
    ldl::DistMultiVecNodeMeta meta;
    return RegularizedSolveAfterPromote
           ( A, reg, invMap, info, front, B, meta,
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
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfterPromote"))
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
      [&]( const DistMultiVec<PF>& xProm, DistMultiVec<PF>& yProm )
      {
        yProm = xProm;
        DiagonalScale( LEFT, NORMAL, regProm, yProm ); 
        Multiply( NORMAL, PF(1), AProm, xProm, PF(1), yProm );
      };
    auto applyAInv = 
      [&]( DistMultiVec<F>& y )
      {
        DiagonalSolve( LEFT, NORMAL, d, y );
        ldl::DistMultiVecNode<F> yNodal( invMap, info, y, meta );
        ldl::SolveAfter( info, front, yNodal );
        yNodal.Push( invMap, info, y, meta );
        DiagonalSolve( LEFT, NORMAL, d, y );
      };

    // TODO: Perform these in a batch instead
    const Int height = B.Height();
    const Int width = B.Width();
    Int mostRefineIts = 0;
    DistMultiVec<F> u(B.Comm());
    Zeros( u, height, 1 );
    auto& BLoc = B.Matrix();
    auto& uLoc = u.Matrix();
    for( Int j=0; j<width; ++j )
    {
        auto bLoc = BLoc( ALL, IR(j) );
        Copy( bLoc, uLoc );
        const Int refineIts =
          PromotedIterativeRefinement
          ( applyA, applyAInv, u, relTol, maxRefineIts, progress );
        Copy( uLoc, bLoc );
        mostRefineIts = Max(mostRefineIts,refineIts);
    }

    return mostRefineIts;
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfterPromote"))
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfter"))
#ifdef EL_HAVE_QUAD
    return RegularizedSolveAfterPromote
    ( A, reg, invMap, info, front, B, meta,
      relTol, maxRefineIts, progress, time );
#else
    return RegularizedSolveAfterNoPromote
    ( A, reg, invMap, info, front, B, meta,
      relTol, maxRefineIts, progress, time );
#endif
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfter"))
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfter"))
#ifdef EL_HAVE_QUAD
    return RegularizedSolveAfterPromote
    ( A, reg, d, invMap, info, front, B, meta,
      relTol, maxRefineIts, progress, time );
#else
    return RegularizedSolveAfterNoPromote
    ( A, reg, d, invMap, info, front, B, meta,
      relTol, maxRefineIts, progress, time );
#endif
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
    DEBUG_ONLY(CSE cse("reg_ldl::RegularizedSolveAfter"))
    ldl::DistMultiVecNodeMeta meta;
    return RegularizedSolveAfter
    ( A, reg, d, invMap, info, front, B, meta,
      relTol, maxRefineIts, progress, time );
}

// LEFT OFF HERE

template<typename F>
inline Int LGMRESSolveAfterSingle
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& reg,
  const vector<Int>& invMap, 
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& b,
  Base<F> relTol,
  Int restart,
  Int maxIts,
  Base<F> relTolRefine,
  Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(
      CSE cse("reg_ldl::LGMRESSolveAfterSingle");
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    typedef Base<F> Real;
    const Int n = A.Height();

    // x := 0
    // ======
    Matrix<F> x;
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    Matrix<F> w; 
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( origResidNorm == Real(0) )
        return 0;

    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    Matrix<F> x0, q, V;
    while( !converged )
    {
        if( progress )
            Output("Starting GMRES iteration ",iter);
        const Int indent = PushIndent();

        Zeros( cs, restart, 1 );
        Zeros( sn, restart, 1 );
        Zeros( H,  restart, restart );
        Zeros( V, n, restart );
        
        // x0 := x
        // =======
        x0 = x;

        // w := inv(M) w
        // =============
        Int refineIts = RegularizedSolveAfter
        ( A, reg, invMap, info, front, w, 
          relTolRefine, maxRefineIts, progress );
        maxLargeRefines = Max( refineIts, maxLargeRefines );

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0 = V( ALL, IR(0) );
        v0 = w;
        v0 *= 1/beta;

        // t := beta e_0
        // =============
        Zeros( t, restart+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(restart)
        // =========================
        for( Int j=0; j<restart; ++j )
        {
            if( progress )
                Output("Starting inner GMRES iteration ",j);
            const Int innerIndent = PushIndent();

            // w := A v_j
            // ----------
            Multiply( NORMAL, F(1), A, V(ALL,IR(j)), F(0), w );

            // w := inv(M) w
            // -------------
            Int refineIts = RegularizedSolveAfter
            ( A, reg, invMap, info, front, w, 
              relTolRefine, maxRefineIts, progress );
            maxLargeRefines = Max( refineIts, maxLargeRefines );

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                auto vi = V( ALL, IR(i) );
                H.Set( i, j, Dot(vi,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), vi, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                restart = j+1;
            if( j+1 != restart )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto vjp1 = V( ALL, IR(j+1) );
                vjp1 = w;
                vjp1 *= 1/delta;
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Vj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                const F eta_i = y.Get(i,0);
                Axpy( eta_i, V( ALL, IR(i) ), x );
            }

            // w := b - A x
            // ------------
            w = b;
            Multiply( NORMAL, F(-1), A, x, F(1), w );

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress )
                    Output("converged with relative tolerance: ",relResidNorm);
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress )
                    Output
                    ("finished iteration ",iter," with relResidNorm=",
                     relResidNorm);
            }
            ++iter;
            if( iter == maxIts )
                RuntimeError("LGMRES did not converge");
            SetIndent( innerIndent );
        }
        SetIndent( indent );
    }
    b = x;
    return maxLargeRefines;
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
    DEBUG_ONLY(CSE cse("reg_ldl::LGMRESSolveAfter"))
    Int mostRefineIts = 0;
    const Int width = B.Width();
    for( Int j=0; j<width; ++j )
    {
        auto b = B( ALL, IR(j) );
        const Int refineIts = 
          LGMRESSolveAfterSingle
          (A,reg,invMap,info,front,b,
           relTol,restart,maxIts,relTolRefine,maxRefineIts,progress);
        mostRefineIts = Max(mostRefineIts,refineIts);
    }
    return mostRefineIts;
}

template<typename F>
inline Int LGMRESSolveAfterSingle
( const SparseMatrix<F>& A,
  const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d,
  const vector<Int>& invMap,
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& b,
  Base<F> relTol,
  Int restart,
  Int maxIts,
  Base<F> relTolRefine,
  Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(
      CSE cse("reg_ldl::LGMRESSolveAfterSingle");
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    typedef Base<F> Real;
    const Int n = A.Height();

    // x := 0
    // ======
    Matrix<F> x;
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    Matrix<F> w; 
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( origResidNorm == Real(0) )
        return 0;

    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    Matrix<F> x0, q, V;
    while( !converged )
    {
        if( progress )
            Output("Starting GMRES iteration ",iter);
        const Int indent = PushIndent();

        Zeros( cs, restart, 1 );
        Zeros( sn, restart, 1 );
        Zeros( H,  restart, restart );
        Zeros( V, n, restart );
        
        // x0 := x
        // =======
        x0 = x;

        // w := inv(M) w
        // =============
        Int refineIts = RegularizedSolveAfter
        ( A, reg, d, invMap, info, front, w, 
          relTolRefine, maxRefineIts, progress );
        maxLargeRefines = Max( refineIts, maxLargeRefines );

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0 = V( ALL, IR(0) );
        v0 = w;
        v0 *= 1/beta;

        // t := beta e_0
        // =============
        Zeros( t, restart+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(restart)
        // =========================
        for( Int j=0; j<restart; ++j )
        {
            if( progress )
                Output("Starting inner GMRES iteration ",j);
            const Int innerIndent = PushIndent();

            // w := A v_j
            // ----------
            Multiply( NORMAL, F(1), A, V(ALL,IR(j)), F(0), w );

            // w := inv(M) w
            // -------------
            Int refineIts = RegularizedSolveAfter
            ( A, reg, d, invMap, info, front, w, 
              relTolRefine, maxRefineIts, progress );
            maxLargeRefines = Max( refineIts, maxLargeRefines );

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                auto vi = V( ALL, IR(i) );
                H.Set( i, j, Dot(vi,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), vi, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                restart = j+1;
            if( j+1 != restart )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto vjp1 = V( ALL, IR(j+1) );
                vjp1 = w;
                vjp1 *= 1/delta;
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Vj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                const F eta_i = y.Get(i,0);
                Axpy( eta_i, V( ALL, IR(i) ), x );
            }

            // w := b - A x
            // ------------
            w = b;
            Multiply( NORMAL, F(-1), A, x, F(1), w );

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress )
                    Output("converged with relative tolerance: ",relResidNorm);
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress )
                    Output
                    ("finished iteration ",iter," with relResidNorm=",
                     relResidNorm);
            }
            ++iter;
            if( iter == maxIts )
                RuntimeError("LGMRES did not converge");
            SetIndent( innerIndent );
        }
        SetIndent( indent );
    }
    b = x;
    return maxLargeRefines;
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
    DEBUG_ONLY(CSE cse("reg_ldl::LGMRESSolveAfter"))
    Int mostRefineIts = 0;
    const Int width = B.Width();
    for( Int j=0; j<width; ++j )
    {
        auto b = B( ALL, IR(j) );
        const Int refineIts = 
          LGMRESSolveAfterSingle
          (A,reg,d,invMap,info,front,b,
           relTol,restart,maxIts,relTolRefine,maxRefineIts,progress);
        mostRefineIts = Max(mostRefineIts,refineIts);
    }
    return mostRefineIts;
}

template<typename F>
inline Int LGMRESSolveAfterSingle
( const DistSparseMatrix<F>& A,
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& b,
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int restart,
  Int maxIts,
  Base<F> relTolRefine,
  Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(
      CSE cse("reg_ldl::LGMRESSolveAfterSingle");
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    typedef Base<F> Real;
    const Int n = A.Height();
    mpi::Comm comm = A.Comm();
    const int commRank = mpi::Rank(comm);

    // x := 0
    // ======
    DistMultiVec<F> x(comm);
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    DistMultiVec<F> w(comm); 
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( origResidNorm == Real(0) )
        return 0;

    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    DistMultiVec<F> x0(comm), q(comm), V(comm);
    while( !converged )
    {
        if( progress && commRank == 0 )
            Output("Starting GMRES iteration ",iter);
        const Int indent = PushIndent();

        Zeros( cs, restart, 1 );
        Zeros( sn, restart, 1 );
        Zeros( H,  restart, restart );
        Zeros( V, n, restart );
        // TODO: Extend DistMultiVec so that it can be directly manipulated
        //       rather than requiring access to the local Matrix and staging
        //       through the temporary vector q
        auto& VLoc = V.Matrix();
        Zeros( q, n, 1 );
        
        // x0 := x
        // =======
        x0 = x;

        // w := inv(M) w
        // =============
        Int refineIts = RegularizedSolveAfter
        ( A, reg, invMap, info, front, w, meta,
          relTolRefine, maxRefineIts, progress );
        maxLargeRefines = Max( refineIts, maxLargeRefines );

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0Loc = VLoc( ALL, IR(0) );
        v0Loc = w.Matrix();
        v0Loc *= 1/beta;

        // t := beta e_0
        // =============
        Zeros( t, restart+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(restart)
        // ===============================
        for( Int j=0; j<restart; ++j )
        {
            if( progress && commRank == 0 )
                Output("Starting inner GMRES iteration ",j);
            const Int innerIndent = PushIndent();

            // w := A v_j
            // ----------
            q.Matrix() = VLoc( ALL, IR(j) );
            Multiply( NORMAL, F(1), A, q, F(0), w );

            // w := inv(M) w
            // -------------
            Int refineIts = RegularizedSolveAfter
            ( A, reg, invMap, info, front, w, meta,
              relTolRefine, maxRefineIts, progress );
            maxLargeRefines = Max( refineIts, maxLargeRefines );

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                q.Matrix() = VLoc( ALL, IR(i) );
                H.Set( i, j, Dot(q,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), q, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                restart = j+1;
            if( j+1 != restart )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto v_jp1Loc = VLoc( ALL, IR(j+1) );
                v_jp1Loc = w.Matrix();
                v_jp1Loc *= 1/delta;
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Vj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                const F eta_i = y.Get(i,0);
                Axpy( eta_i, VLoc( ALL, IR(i) ), x.Matrix() );
            }

            // w := b - A x
            // ------------
            w = b;
            Multiply( NORMAL, F(-1), A, x, F(1), w );

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress && commRank == 0 )
                    Output("converged with relative tolerance: ",relResidNorm);
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress && commRank == 0 )
                    Output
                    ("finished iteration ",iter," with relResidNorm=",
                     relResidNorm);
            }
            ++iter;
            if( iter == maxIts )
                RuntimeError("LGMRES did not converge");
            SetIndent( innerIndent );
        }
        SetIndent( indent );
    }
    b = x;
    return maxLargeRefines;
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
    DEBUG_ONLY(CSE cse("reg_ldl::LGMRESSolveAfter"))
    const Int m = B.Height();
    const Int n = B.Width();
    mpi::Comm comm = A.Comm();

    Int mostRefineIts = 0;
    DistMultiVec<F> u(comm);
    Zeros( u, m, 1 );
    auto& BLoc = B.Matrix();
    auto& uLoc = u.Matrix();
    for( Int j=0; j<n; ++j )
    {
        auto bLoc = BLoc( ALL, IR(j) );
        Copy( bLoc, uLoc );
        const Int refineIts = LGMRESSolveAfterSingle
          ( A, reg, invMap, info, front, u, meta,
            relTol, restart, maxIts, relTolRefine, maxRefineIts, progress );
        Copy( uLoc, bLoc );
        mostRefineIts = Max(mostRefineIts,refineIts);
    }
    return mostRefineIts;
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
    DEBUG_ONLY(CSE cse("reg_ldl::LGMRESSolveAfter"))
    ldl::DistMultiVecNodeMeta meta;
    return LGMRESSolveAfter
           ( A, reg, invMap, info, front, B, meta,
             relTol, restart, maxIts, relTolRefine, maxRefineIts, progress );
}

template<typename F>
inline Int LGMRESSolveAfterSingle
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& b,
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int restart,
  Int maxIts,
  Base<F> relTolRefine,
  Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(
      CSE cse("reg_ldl::LGMRESSolveAfterSingle");
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    typedef Base<F> Real;
    const Int n = A.Height();
    mpi::Comm comm = A.Comm();
    const int commRank = mpi::Rank(comm);

    // x := 0
    // ======
    DistMultiVec<F> x(comm);
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    DistMultiVec<F> w(comm); 
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( origResidNorm == Real(0) )
        return 0;

    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    DistMultiVec<F> x0(comm), q(comm), V(comm);
    while( !converged )
    {
        if( progress && commRank == 0 )
            Output("Starting GMRES iteration ",iter);
        const Int indent = PushIndent();

        Zeros( cs, restart, 1 );
        Zeros( sn, restart, 1 );
        Zeros( H,  restart, restart );
        Zeros( V, n, restart );
        // TODO: Extend DistMultiVec so that it can be directly manipulated
        //       rather than requiring access to the local Matrix and staging
        //       through the temporary vector q
        auto& VLoc = V.Matrix();
        Zeros( q, n, 1 );
        
        // x0 := x
        // =======
        x0 = x;

        // w := inv(M) w
        // =============
        Int refineIts = RegularizedSolveAfter
        ( A, reg, d, invMap, info, front, w, meta,
          relTolRefine, maxRefineIts, progress );
        maxLargeRefines = Max( refineIts, maxLargeRefines );

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0Loc = VLoc( ALL, IR(0) );
        v0Loc = w.Matrix();
        v0Loc *= 1/beta;

        // t := beta e_0
        // =============
        Zeros( t, restart+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(restart)
        // ===============================
        for( Int j=0; j<restart; ++j )
        {
            if( progress && commRank )
                Output("Starting inner GMRES iteration ",j);
            const Int innerIndent = PushIndent();

            // w := A v_j
            // ----------
            q.Matrix() = VLoc( ALL, IR(j) );
            Multiply( NORMAL, F(1), A, q, F(0), w );

            // w := inv(M) w
            // -------------
            Int refineIts = RegularizedSolveAfter
            ( A, reg, d, invMap, info, front, w, meta,
              relTolRefine, maxRefineIts, progress );
            maxLargeRefines = Max( refineIts, maxLargeRefines );

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                q.Matrix() = VLoc( ALL, IR(i) );
                H.Set( i, j, Dot(q,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), q, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                restart = j+1;
            if( j+1 != restart )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto v_jp1Loc = VLoc( ALL, IR(j+1) );
                v_jp1Loc = w.Matrix();
                v_jp1Loc *= 1/delta;
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Vj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                const F eta_i = y.Get(i,0);
                Axpy( eta_i, VLoc( ALL, IR(i) ), x.Matrix() );
            }

            // w := b - A x
            // ------------
            w = b;
            Multiply( NORMAL, F(-1), A, x, F(1), w );

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress && commRank == 0 )
                    Output("converged with relative tolerance: ",relResidNorm);
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress && commRank == 0 )
                    Output
                    ("finished iteration ",iter," with relResidNorm=",
                     relResidNorm);
            }
            ++iter;
            if( iter == maxIts )
                RuntimeError("LGMRES did not converge");
            SetIndent( innerIndent );
        }
        SetIndent( indent );
    }
    b = x;
    return maxLargeRefines;
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
    DEBUG_ONLY(CSE cse("reg_ldl::LGMRESSolveAfter"))
    const Int m = B.Height();
    const Int n = B.Width();
    mpi::Comm comm = A.Comm();

    Int mostRefineIts = 0;
    DistMultiVec<F> u(comm);
    Zeros( u, m, 1 );
    auto& BLoc = B.Matrix();
    auto& uLoc = u.Matrix();
    for( Int j=0; j<n; ++j )
    {
        auto bLoc = BLoc( ALL, IR(j) );
        Copy( bLoc, uLoc );
        const Int refineIts = LGMRESSolveAfterSingle
          ( A, reg, d, invMap, info, front, u, meta,
            relTol, restart, maxIts, relTolRefine, maxRefineIts, progress );
        Copy( uLoc, bLoc );
        mostRefineIts = Max(mostRefineIts,refineIts);
    }
    return mostRefineIts;
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
    DEBUG_ONLY(CSE cse("reg_ldl::LGMRESSolveAfter"))
    ldl::DistMultiVecNodeMeta meta;
    return LGMRESSolveAfter
           ( A, reg, d, invMap, info, front, B, meta,
             relTol, restart, maxIts, relTolRefine, maxRefineIts, progress );
}

// The pseudocode for Flexible GMRES can be found in "Algorithm 2.2" in
//   Youcef Saad
//   "A flexible inner-outer preconditioned GMRES algorithm"
//   SIAM J. Sci. Comput., Vol. 14, No. 2, pp. 461--469, 1993.

template<typename F>
inline Int FGMRESSolveAfterSingle
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& reg,
  const vector<Int>& invMap, 
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& b,
  Base<F> relTol,
  Int restart,
  Int maxIts,
  Base<F> relTolRefine,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    DEBUG_ONLY(
      CSE cse("reg_ldl::FGMRESSolveAfterSingle");
      if( b.Width() != 1 )    
          LogicError("Expected a single right-hand side");
    )

    // Avoid half of the matrix-vector products by keeping the results of 
    // A z_j and A x_0
    const bool saveProducts = true;

    typedef Base<F> Real;
    const Int n = A.Height();
    Timer iterTimer, timer;

    // x := 0
    // ======
    Matrix<F> x;
    Zeros( x, n, 1 );

    Matrix<F> Ax0;
    if( saveProducts )
    {
        // A x_0 := 0
        // ==========
        Zeros( Ax0, n, 1 );
    }

    // w := b (= b - A x_0)
    // ====================
    auto w = b;
    const Real origResidNorm = Nrm2( w );
    if( progress )
        Output("origResidNorm: ",origResidNorm);
    if( origResidNorm == Real(0) )
        return 0;

    // TODO: Constrain the maximum number of iterations

    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    Matrix<F> x0, V, Z, AZ, q;
    while( !converged )
    {
        if( progress )
            Output("Starting FGMRES iteration ",iter);
        const Int indent = PushIndent();

        Zeros( cs, restart, 1 );
        Zeros( sn, restart, 1 );
        Zeros( H,  restart, restart );
        Zeros( V, n, restart );
        Zeros( Z, n, restart );
        if( saveProducts )
            Zeros( AZ, n, restart );
        
        // x0 := x
        // =======
        x0 = x;
        if( saveProducts && iter != 0 )
            Ax0 = q;

        // NOTE: w = b - A x already

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0 = V( ALL, IR(0) );
        v0 = w;
        v0 *= 1/beta;

        // t := beta e_0
        // =============
        Zeros( t, restart+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(restart)
        // ===============================
        for( Int j=0; j<restart; ++j )
        {
            if( progress )
                Output("Starting inner FGMRES iteration ",j);
            if( time )
                iterTimer.Start();
            const Int innerIndent = PushIndent();

            // z_j := inv(M) v_j
            // =================
            if( time )
                timer.Start();
            auto vj = V( ALL, IR(j) );
            auto zj = Z( ALL, IR(j) );
            zj = vj;
            Int refineIts = RegularizedSolveAfter
            ( A, reg, invMap, info, front, zj, 
              relTolRefine, maxRefineIts, progress, time );
            maxLargeRefines = Max( refineIts, maxLargeRefines );
            if( time )
                Output("solve took ",timer.Stop()," secs");

            // w := A z_j
            // ----------
            if( time )
                timer.Start();
            Multiply( NORMAL, F(1), A, zj, F(0), w );
            if( time )
                Output("mat-vec took ",timer.Stop()," secs");
            if( saveProducts )
            {
                auto Azj = AZ( ALL, IR(j) );
                Azj = w;
            }

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                auto vi = V( ALL, IR(i) );
                H.Set( i, j, Dot(vi,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), vi, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                restart = j+1;
            if( j+1 != restart )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto vjp1 = V( ALL, IR(j+1) );
                vjp1 = w;
                vjp1 *= 1/delta;
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Zj y
            // ^^^^^^^^^^^^^^
            x = x0;
            auto Zj = Z( ALL, IR(0,j+1) );
            auto yj = y( IR(0,j+1), ALL );
            Gemv( NORMAL, F(1), Zj, yj, F(1), x );

            // w := b - A x
            // ------------
            if( time )
                timer.Start();
            w = b;
            if( saveProducts )
            {
                // q := Ax = Ax0 + A Z_j y_j
                // ^^^^^^^^^^^^^^^^^^^^^^^^^
                q = Ax0;
                auto AZj = AZ( ALL, IR(0,j+1) );
                Gemv( NORMAL, F(1), AZj, yj, F(1), q );

                // w := b - A x
                // ^^^^^^^^^^^^
                w -= q;
            }
            else
            {
                Multiply( NORMAL, F(-1), A, x, F(1), w );
            }
            if( time )
                Output("mat-vec took ",timer.Stop()," secs");

            if( time )
                Output("iter took ",iterTimer.Stop()," secs");

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress )
                    Output("converged with relative tolerance: ",relResidNorm);
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress )
                    Output
                    ("finished iteration ",iter," with relResidNorm=",
                     relResidNorm);
            }
            ++iter;
            if( iter == maxIts )
                RuntimeError("FGMRES did not converge");
            SetIndent( innerIndent );
        }
        SetIndent( indent );
    }
    b = x;
    return maxLargeRefines;
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
    DEBUG_ONLY(CSE cse("reg_ldl::FGMRESSolveAfter"))
    Int mostRefineIts = 0;
    const Int width = B.Width();
    for( Int j=0; j<width; ++j )
    {
        auto b = B( ALL, IR(j) );
        const Int refineIts = 
          FGMRESSolveAfterSingle
          (A,reg,invMap,info,front,b,
           relTol,restart,maxIts,relTolRefine,maxRefineIts,progress,time);
        mostRefineIts = Max(mostRefineIts,refineIts);
    }
    return mostRefineIts;
}

template<typename F>
inline Int FGMRESSolveAfterSingle
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d,
  const vector<Int>& invMap, 
  const ldl::NodeInfo& info,
  const ldl::Front<F>& front, 
        Matrix<F>& b,
  Base<F> relTol,
  Int restart,
  Int maxIts,
  Base<F> relTolRefine,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    DEBUG_ONLY(
      CSE cse("reg_ldl::FGMRESSolveAfterSingle");
      if( b.Width() != 1 )    
          LogicError("Expected a single right-hand side");
    )

    // Avoid half of the matrix-vector products by keeping the results of
    // A z_j and A x_0
    const bool saveProducts = true;

    typedef Base<F> Real;
    const Int n = A.Height();

    // x := 0
    // ======
    Matrix<F> x;
    Zeros( x, n, 1 );

    Matrix<F> Ax0;
    if( saveProducts )
    {
        // A x_0 := 0
        // ==========
        Zeros( Ax0, n, 1 );
    }

    // w := b (= b - A x_0)
    // ====================
    auto w = b;
    const Real origResidNorm = Nrm2( w );
    if( progress )
        Output("origResidNorm: ",origResidNorm);
    if( origResidNorm == Real(0) )
        return 0;

    // TODO: Constrain the maximum number of iterations
    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    Matrix<F> x0, V, Z, AZ, q;
    while( !converged )
    {
        if( progress )
            Output("Starting FGMRES iteration ",iter);
        const Int indent = PushIndent();

        Zeros( cs, restart, 1 );
        Zeros( sn, restart, 1 );
        Zeros( H,  restart, restart );
        Zeros( V, n, restart );
        Zeros( Z, n, restart );
        if( saveProducts )
            Zeros( AZ, n, restart );
        
        // x0 := x
        // =======
        x0 = x;
        if( saveProducts && iter != 0 )
            Ax0 = q;

        // NOTE: w = b - A x already

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0 = V( ALL, IR(0) );
        v0 = w;
        v0 *= 1/beta;

        // t := beta e_0
        // =============
        Zeros( t, restart+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(restart)
        // ===============================
        for( Int j=0; j<restart; ++j )
        {
            if( progress )
                Output("Starting inner FGMRES iteration ",j);
            const Int innerIndent = PushIndent();

            // z_j := inv(M) v_j
            // =================
            auto vj = V( ALL, IR(j) );
            auto zj = Z( ALL, IR(j) );
            zj = vj;
            Int refineIts = RegularizedSolveAfter
            ( A, reg, d, invMap, info, front, zj, 
              relTolRefine, maxRefineIts, progress, time );
            maxLargeRefines = Max( refineIts, maxLargeRefines );

            // w := A z_j
            // ----------
            Multiply( NORMAL, F(1), A, zj, F(0), w );
            if( saveProducts )
            {
                auto Azj = AZ( ALL, IR(j) );
                Azj = w;
            }

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                auto vi = V( ALL, IR(i) );
                H.Set( i, j, Dot(vi,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), vi, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                restart = j+1; 
            if( j+1 != restart )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto vjp1 = V( ALL, IR(j+1) );
                vjp1 = w;
                vjp1 *= 1/delta;
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Zj y
            // ^^^^^^^^^^^^^^
            x = x0;
            auto Zj = Z( ALL, IR(0,j+1) );
            auto yj = y( IR(0,j+1), ALL );
            Gemv( NORMAL, F(1), Zj, yj, F(1), x );

            // w := b - A x
            // ------------
            w = b;
            if( saveProducts )
            {
                // q := A x = Ax0 + A Z_j y_j
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                q = Ax0;
                auto AZj = AZ( ALL, IR(0,j+1) );
                Gemv( NORMAL, F(1), AZj, yj, F(1), q );

                // w := b - A x
                // ^^^^^^^^^^^^
                w -= q;
            }
            else
            {
                Multiply( NORMAL, F(-1), A, x, F(1), w );
            }

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress )
                    Output("converged with relative tolerance: ",relResidNorm);
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress )
                    Output
                    ("finished iteration ",iter," with relResidNorm=",
                     relResidNorm);
            }
            ++iter;
            if( iter == maxIts )
                RuntimeError("FGMRES did not converge");
            SetIndent( innerIndent );
        }
        SetIndent( indent );
    }
    b = x;
    return maxLargeRefines;
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
    DEBUG_ONLY(CSE cse("reg_ldl::FGMRESSolveAfter"))
    Int mostRefineIts = 0;
    const Int width = B.Width();
    for( Int j=0; j<width; ++j )
    {
        auto b = B( ALL, IR(j) );
        const Int refineIts = 
          FGMRESSolveAfterSingle
          (A,reg,d,invMap,info,front,b,
           relTol,restart,maxIts,relTolRefine,maxRefineIts,progress,time);
        mostRefineIts = Max(mostRefineIts,refineIts);
    }
    return mostRefineIts;
}

template<typename F>
inline Int FGMRESSolveAfterSingle
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& b,
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int restart,
  Int maxIts,
  Base<F> relTolRefine,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    DEBUG_ONLY(
      CSE cse("reg_ldl::FGMRESSolveAfterSingle");
      if( b.Width() != 1 )    
          LogicError("Expected a single right-hand side");
    )

    // Avoid half of the matrix-vector products by keeping the results of
    // A z_j and A x_0
    const bool saveProducts = true;

    typedef Base<F> Real;
    const Int n = A.Height();
    mpi::Comm comm = A.Comm();
    const int commRank = mpi::Rank(comm);
    Timer iterTimer, timer;

    // x := 0
    // ======
    DistMultiVec<F> x(comm);
    Zeros( x, n, 1 );

    DistMultiVec<F> Ax0(comm);
    if( saveProducts )
    {
        // A x_0 := 0
        // ==========
        Zeros( Ax0, n, 1 );
    }

    // w := b (= b - A x_0)
    // ====================
    DistMultiVec<F> w(comm); 
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( progress && commRank == 0 )
        Output("origResidNorm: ",origResidNorm);
    if( origResidNorm == Real(0) )
        return 0;

    // TODO: Constrain the maximum number of iterations
    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    DistMultiVec<F> x0(comm), q(comm), V(comm), Z(comm), AZ(comm);
    while( !converged )
    {
        if( progress && commRank == 0 )
            Output("Starting FGMRES iteration ",iter);
        const Int indent = PushIndent();

        // x0 := x
        // =======
        x0 = x;
        if( saveProducts && iter != 0 )
            Ax0 = q;

        Zeros( cs, restart, 1 );
        Zeros( sn, restart, 1 );
        Zeros( H,  restart, restart );
        Zeros( V, n, restart );
        Zeros( Z, n, restart );
        if( saveProducts )
            Zeros( AZ, n, restart );         

        // TODO: Extend DistMultiVec so that it can be directly manipulated
        //       rather than requiring access to the local Matrix and staging
        //       through the temporary vector q
        auto& VLoc = V.Matrix();
        auto& ZLoc = Z.Matrix();
        Zeros( q, n, 1 );
        
        // NOTE: w = b - A x already

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0Loc = VLoc( ALL, IR(0) );
        v0Loc = w.Matrix();
        v0Loc *= 1/beta;

        // t := beta e_0
        // =============
        Zeros( t, restart+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(restart)
        // ===============================
        for( Int j=0; j<restart; ++j )
        {
            if( progress && commRank == 0 )
                Output("Starting inner FGMRES iteration ",j);
            if( time && commRank == 0 )
                iterTimer.Start();
            const Int innerIndent = PushIndent();

            // z_j := inv(M) v_j
            // =================
            if( time && commRank == 0 )
                timer.Start();
            auto vjLoc = VLoc( ALL, IR(j) );
            auto zjLoc = ZLoc( ALL, IR(j) );
            q.Matrix() = vjLoc;
            Int refineIts = RegularizedSolveAfter
            ( A, reg, invMap, info, front, q, meta,
              relTolRefine, maxRefineIts, progress, time );
            maxLargeRefines = Max( refineIts, maxLargeRefines );
            zjLoc = q.Matrix();
            if( time && commRank == 0 )
                Output("solve took ",timer.Stop()," secs");

            // w := A z_j
            // ----------
            if( time && commRank == 0 )
                timer.Start();
            // NOTE: q currently contains z_j
            Multiply( NORMAL, F(1), A, q, F(0), w );
            if( time && commRank == 0 )
                Output("mat-vec took ",timer.Stop()," secs");
            if( saveProducts )
            {
                auto& AZLoc = AZ.Matrix();
                auto AzjLoc = AZLoc( ALL, IR(j) );
                AzjLoc = w.LockedMatrix();
            }

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                q.Matrix() = VLoc( ALL, IR(i) );
                H.Set( i, j, Dot(q,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), q, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                restart = j+1;
            if( j+1 != restart )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto v_jp1Loc = VLoc( ALL, IR(j+1) );
                v_jp1Loc = w.Matrix();
                v_jp1Loc *= 1/delta;
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Zj y
            // ^^^^^^^^^^^^^^
            x = x0;
            auto ZjLoc = ZLoc( ALL, IR(0,j+1) );
            auto yj = y( IR(0,j+1), ALL );
            Gemv( NORMAL, F(1), ZjLoc, yj, F(1), x.Matrix() );

            // w := b - A x
            // ------------
            if( time && commRank == 0 )
                timer.Start();
            w = b;
            if( saveProducts )
            {
                // q := Ax = Ax0 + A Z_j y_j
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^^
                q = Ax0;
                const auto& AZLoc = AZ.LockedMatrix();
                auto AZjLoc = AZLoc( ALL, IR(0,j+1) );
                Gemv( NORMAL, F(1), AZjLoc, yj, F(1), q.Matrix() );

                // w := b - A x
                // ^^^^^^^^^^^^
                w -= q;
            }
            else
            {
                Multiply( NORMAL, F(-1), A, x, F(1), w );
            }
            if( time && commRank == 0 )
                Output("mat-vec took ",timer.Stop()," secs");

            if( time && commRank == 0 )
                Output("iter took ",iterTimer.Stop()," secs");

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress && commRank == 0 )
                    Output("converged with relative tolerance: ",relResidNorm);
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress && commRank == 0 )
                    Output
                    ("finished iteration ",iter," with relResidNorm=",
                     relResidNorm);
            }
            ++iter;
            if( iter == maxIts )
                RuntimeError("FGMRES did not converge");
            SetIndent( innerIndent );
        }
        SetIndent( indent );
    }
    b = x;
    return maxLargeRefines;
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
    DEBUG_ONLY(CSE cse("reg_ldl::FGMRESSolveAfter"))
    const Int m = B.Height();
    const Int n = B.Width();
    mpi::Comm comm = A.Comm();

    Int mostRefineIts = 0;
    DistMultiVec<F> u(comm);
    Zeros( u, m, 1 );
    auto& BLoc = B.Matrix();
    auto& uLoc = u.Matrix();
    for( Int j=0; j<n; ++j )
    {
        auto bLoc = BLoc( ALL, IR(j) );
        Copy( bLoc, uLoc );
        const Int refineIts = FGMRESSolveAfterSingle
          ( A, reg, invMap, info, front, u, meta,
            relTol, restart, maxIts, relTolRefine, maxRefineIts,
            progress, time );
        Copy( uLoc, bLoc );
        mostRefineIts = Max(mostRefineIts,refineIts);
    }
    return mostRefineIts;
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
    DEBUG_ONLY(CSE cse("reg_ldl::FGMRESSolveAfter"))
    ldl::DistMultiVecNodeMeta meta;
    return FGMRESSolveAfter
           ( A, reg, invMap, info, front, B, meta,
             relTol, restart, maxIts, relTolRefine, maxRefineIts,
             progress, time );
}

template<typename F>
inline Int FGMRESSolveAfterSingle
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap, 
  const ldl::DistNodeInfo& info,
  const ldl::DistFront<F>& front, 
        DistMultiVec<F>& b,
        ldl::DistMultiVecNodeMeta& meta,
  Base<F> relTol,
  Int restart,
  Int maxIts,
  Base<F> relTolRefine,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    DEBUG_ONLY(
      CSE cse("reg_ldl::FGMRESSolveAfterSingle");
      if( b.Width() != 1 )    
          LogicError("Expected a single right-hand side");
    )

    // Avoid half of the matrix-vector products by keeping the results of
    // A z_j and A x_0
    const bool saveProducts = true;

    typedef Base<F> Real;
    const Int n = A.Height();
    mpi::Comm comm = A.Comm();
    const int commRank = mpi::Rank(comm);
    Timer iterTimer, timer;

    // x := 0
    // ======
    DistMultiVec<F> x(comm);
    Zeros( x, n, 1 );

    DistMultiVec<F> Ax0(comm);
    if( saveProducts )
    {
        // A x_0 := 0
        // ==========
        Zeros( Ax0, n, 1 );
    }

    // w := b (= b - A x_0)
    // ====================
    DistMultiVec<F> w(comm); 
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( progress && commRank == 0 )
        Output("origResidNorm: ",origResidNorm);
    if( origResidNorm == Real(0) )
        return 0;

    // TODO: Constrain the maximum number of iterations
    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    DistMultiVec<F> x0(comm), q(comm), V(comm), Z(comm), AZ(comm);
    while( !converged )
    {
        if( progress && commRank == 0 )
            Output("Starting FGMRES iteration ",iter);
        const Int indent = PushIndent();

        // x0 := x
        // =======
        x0 = x;
        if( saveProducts && iter != 0 )
            Ax0 = q;

        Zeros( cs, restart, 1 );
        Zeros( sn, restart, 1 );
        Zeros( H,  restart, restart );
        Zeros( V, n, restart );
        Zeros( Z, n, restart );
        if( saveProducts )
            Zeros( AZ, n, restart );

        // TODO: Extend DistMultiVec so that it can be directly manipulated
        //       rather than requiring access to the local Matrix and staging
        //       through the temporary vector q
        auto& VLoc = V.Matrix();
        auto& ZLoc = Z.Matrix();
        Zeros( q, n, 1 );
        
        // NOTE: w = b - A x already

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0Loc = VLoc( ALL, IR(0) );
        v0Loc = w.Matrix();
        v0Loc *= 1/beta;

        // t := beta e_0
        // =============
        Zeros( t, restart+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(restart)
        // ===============================
        for( Int j=0; j<restart; ++j )
        {
            if( progress && commRank == 0 )
                Output("Starting inner FGMRES iteration ",j);
            if( time && commRank == 0 )
                iterTimer.Start();
            const Int innerIndent = PushIndent();

            // z_j := inv(M) v_j
            // =================
            if( time && commRank == 0 )
                timer.Start();
            auto vjLoc = VLoc( ALL, IR(j) );
            auto zjLoc = ZLoc( ALL, IR(j) );
            q.Matrix() = vjLoc;
            Int refineIts = RegularizedSolveAfter
            ( A, reg, d, invMap, info, front, q, meta,
              relTolRefine, maxRefineIts, progress, time );
            maxLargeRefines = Max( refineIts, maxLargeRefines );
            zjLoc = q.Matrix();
            if( time && commRank == 0 )
                Output("solve took ",timer.Stop()," secs");

            // w := A z_j
            // ----------
            if( time && commRank == 0 )
                timer.Start();
            // NOTE: q currently contains z_j
            Multiply( NORMAL, F(1), A, q, F(0), w );
            if( time && commRank == 0 )
                Output("mat-vec took ",timer.Stop()," secs");
            if( saveProducts )
            {
                auto& AZLoc = AZ.Matrix();
                auto AzjLoc = AZLoc( ALL, IR(j) );
                AzjLoc = w.LockedMatrix();
            }

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                q.Matrix() = VLoc( ALL, IR(i) );
                H.Set( i, j, Dot(q,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), q, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                restart = j+1;
            if( j+1 != restart )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto v_jp1Loc = VLoc( ALL, IR(j+1) );
                v_jp1Loc = w.Matrix();
                v_jp1Loc *= 1/delta;
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Zj y
            // ^^^^^^^^^^^^^^
            x = x0;
            auto ZjLoc = ZLoc( ALL, IR(0,j+1) );
            auto yj = y( IR(0,j+1), ALL );
            Gemv( NORMAL, F(1), ZjLoc, yj, F(1), x.Matrix() );

            // w := b - A x
            // ------------
            if( time && commRank == 0 )
                timer.Start();
            w = b;
            if( saveProducts )
            {
                // q := Ax = Ax0 + A Z_j y_j
                // ^^^^^^^^^^^^^^^^^^^^^^^^^
                q = Ax0;
                const auto& AZLoc = AZ.LockedMatrix();
                auto AZjLoc = AZLoc( ALL, IR(0,j+1) );
                Gemv( NORMAL, F(1), AZjLoc, yj, F(1), q.Matrix() );

                // w := b - A x
                // ^^^^^^^^^^^^
                w -= q;
            }
            else
            {
                Multiply( NORMAL, F(-1), A, x, F(1), w );
            }
            if( time && commRank == 0 )
                Output("mat-vec took ",timer.Stop()," secs");

            if( time && commRank == 0 )
                Output("iter took ",iterTimer.Stop()," secs");

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress && commRank == 0 )
                    Output("converged with relative tolerance: ",relResidNorm);
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress && commRank == 0 )
                    Output
                    ("finished iteration ",iter," with relResidNorm=",
                     relResidNorm);
            }
            ++iter;
            if( iter == maxIts )
                RuntimeError("FGMRES did not converge");
            SetIndent( innerIndent );
        }
        SetIndent( indent );
    }
    b = x;
    return maxLargeRefines;
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
    DEBUG_ONLY(CSE cse("reg_ldl::FGMRESSolveAfter"))
    const Int m = B.Height();
    const Int n = B.Width();
    mpi::Comm comm = A.Comm();

    Int mostRefineIts = 0;
    DistMultiVec<F> u(comm);
    Zeros( u, m, 1 );
    // TODO: Batch solve
    auto& BLoc = B.Matrix();
    auto& uLoc = u.Matrix();
    for( Int j=0; j<n; ++j )
    {
        auto bLoc = BLoc( ALL, IR(j) );
        Copy( bLoc, uLoc );
        const Int refineIts = FGMRESSolveAfterSingle
          ( A, reg, d, invMap, info, front, u, meta,
            relTol, restart, maxIts, relTolRefine, maxRefineIts,
            progress, time );
        Copy( uLoc, bLoc );
        mostRefineIts = Max(mostRefineIts,refineIts);
    }
    return mostRefineIts;
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
    DEBUG_ONLY(CSE cse("reg_ldl::FGMRESSolveAfter"))
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
    DEBUG_ONLY(CSE cse("reg_ldl::SolveAfter"))
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
    DEBUG_ONLY(CSE cse("reg_ldl::SolveAfter"))
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
    DEBUG_ONLY(CSE cse("reg_ldl::SolveAfter"))
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
    DEBUG_ONLY(CSE cse("reg_ldl::SolveAfter"))
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
    DEBUG_ONLY(CSE cse("reg_ldl::SolveAfter"))
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
    DEBUG_ONLY(CSE cse("reg_ldl::SolveAfter"))
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
#include "El/macros/Instantiate.h"

} // namespace reg_ldl
} // namespace El
