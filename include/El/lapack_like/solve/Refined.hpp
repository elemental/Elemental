/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SOLVE_REFINED_HPP
#define EL_SOLVE_REFINED_HPP

namespace El {

namespace refined_solve {

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int Single
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<F>& b,
        Base<F> relTol,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_ONLY(
      CSE cse("refined_solve::Single");
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    if( maxRefineIts <= 0 )
    {
        applyAInv( b );
        return 0;
    }

    auto bOrig = b;
    const Base<F> bNorm = MaxNorm( b );

    // Compute the initial guess
    // =========================
    auto x = b;
    applyAInv( x );

    Int refineIt = 0;
    Matrix<F> dx, xCand, y;
    Zeros( y, x.Height(), 1 );
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
        applyA( xCand, y );
        b = bOrig;
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
    b = x;
    return refineIt;
}

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int Batch
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<F>& B,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_ONLY(CSE cse("refined_solve::Batch"))
    if( maxRefineIts <= 0 )
    {
        applyAInv( B );    
        return 0;
    }

    // TODO: Allow for early exits

    // Compute the initial guesses
    // ===========================
    auto BOrig = B;
    auto X = B;
    applyAInv( X );

    Int refineIt = 0;
    Matrix<F> dX, Y;
    Zeros( Y, X.Height(), X.Width() );
    applyA( X, Y );
    B -= Y;

    while( true )
    {
        // Compute the updates to the solutions
        // ------------------------------------
        dX = B;
        applyAInv( dX );
        X += dX;

        ++refineIt;
        if( refineIt < maxRefineIts )
        {
            // Compute the new residual
            // ------------------------
            applyA( X, Y );
            B = BOrig;
            B -= Y;
        }
        else
            break;
    }
    B = X;
    return refineIt;
}

} // namespace refined_solve

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int RefinedSolve
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<F>& B,
        Base<F> relTol,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_ONLY(CSE cse("RefinedSolve"))
    if( B.Width() == 1 )
        return refined_solve::Single
               ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
    else
        return refined_solve::Batch
               ( applyA, applyAInv, B, maxRefineIts, progress );
}

namespace refined_solve {

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int PromotedSingle
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<F>& b,
        Base<F> relTol,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_ONLY(
      CSE cse("refined_solve::PromotedSingle");
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    if( maxRefineIts <= 0 )
    {
        applyAInv( b );
        return 0;
    }
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
    Matrix<PF> dxProm, xCandProm, yProm;
    Zeros( yProm, xProm.Height(), 1 );
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
    Copy( xProm, b );
    return refineIt;
}

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int PromotedBatch
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<F>& B,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_ONLY(CSE cse("refined_solve::PromotedBatch"))
    if( maxRefineIts <= 0 )
    {
        applyAInv( B );
        return 0;
    }
    typedef Promote<F> PF;

    Matrix<PF> BProm, BOrigProm;
    Copy( B, BProm );
    Copy( B, BOrigProm );

    // Compute the initial guesses
    // ===========================
    applyAInv( B );
    Matrix<PF> XProm;
    Copy( B, XProm );

    Int refineIt = 0;
    Matrix<PF> dXProm, YProm;
    Zeros( YProm, XProm.Height(), XProm.Width() );
    applyA( XProm, YProm );
    BProm -= YProm;

    while( true )
    {
        // Update the solutions
        // --------------------
        Copy( BProm, B );
        applyAInv( B );
        Copy( B, dXProm );
        XProm += dXProm;

        ++refineIt;
        if( refineIt < maxRefineIts )
        {
            // Form the new residuals
            // ----------------------
            applyA( XProm, YProm );
            BProm = BOrigProm;
            BProm -= YProm;
        }
        else
            break;
    }
    Copy( XProm, B );
    return refineIt;
}

} // namespace refined_solve

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int PromotedRefinedSolve
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<F>& B,
        Base<F> relTol,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_ONLY(CSE cse("PromotedRefinedSolve"))
    if( B.Width() == 1 )
        return refined_solve::PromotedSingle
               ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
    else
        return refined_solve::PromotedBatch
               ( applyA, applyAInv, B, maxRefineIts, progress );
}

namespace refined_solve {

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int Single
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<F>& b,
        Base<F> relTol,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_ONLY(
      CSE cse("refined_solve::Single");
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    if( maxRefineIts <= 0 )
    {
        applyAInv( b );
        return 0;
    }
    mpi::Comm comm = b.Comm();
    const int commRank = mpi::Rank(comm);

    auto bOrig = b;
    const Base<F> bNorm = MaxNorm( b );

    // Compute the initial guess
    // =========================
    auto x = b;
    applyAInv( x );

    Int refineIt = 0;
    DistMultiVec<F> dx(comm), xCand(comm), y(comm);
    Zeros( y, x.Height(), 1 );
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
        applyA( xCand, y );
        b = bOrig;
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
    b = x;
    return refineIt;
}

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int Batch
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<F>& B,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_ONLY(CSE cse("refined_solve::Batch"))
    if( maxRefineIts <= 0 )
    {
        applyAInv( B );
        return 0;
    }
    mpi::Comm comm = B.Comm();

    auto BOrig = B;

    // Compute the initial guess
    // =========================
    auto X = B;
    applyAInv( X );

    Int refineIt = 0;
    DistMultiVec<F> dX(comm), Y(comm);
    Zeros( Y, X.Height(), X.Width() );
    applyA( X, Y );
    B -= Y;

    const Int indent = PushIndent();
    while( true )
    {
        // Compute the proposed updates to the solutions
        // ---------------------------------------------
        dX = B;
        applyAInv( dX );
        X += dX;

        ++refineIt;
        if( refineIt < maxRefineIts )
        {
            // Compute the new residual
            // ------------------------
            applyA( X, Y );
            B = BOrig;
            B -= Y;
        }
        else
            break;
    }
    SetIndent( indent );
    B = X;
    return refineIt;
}

} // namespace refined_solve

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int RefinedSolve
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<F>& B,
        Base<F> relTol,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_ONLY(CSE cse("RefinedSolve"))
    if( B.Width() == 1 )
        return refined_solve::Single
               ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
    else
        return refined_solve::Batch
               ( applyA, applyAInv, B, maxRefineIts, progress );
}

namespace refined_solve {

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int PromotedSingle
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<F>& b,
        Base<F> relTol,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_ONLY(
      CSE cse("refined_solve::PromotedSingle");
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    if( maxRefineIts <= 0 )
    {
        applyAInv( b );
        return 0;
    }
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
    DistMultiVec<PF> dxProm(comm), xCandProm(comm), yProm(comm);
    Zeros( yProm, xProm.Height(), 1 );
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
    Copy( xProm, b );
    return refineIt;
}

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int PromotedBatch
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<F>& B,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_ONLY(CSE cse("refined_solve::PromotedBatch"))
    if( maxRefineIts <= 0 )
    {
        applyAInv( B );
        return 0;
    }
    typedef Promote<F> PF;
    mpi::Comm comm = B.Comm();
    const int commRank = mpi::Rank(comm);

    DistMultiVec<PF> BProm(comm), BOrigProm(comm);
    Copy( B, BProm );
    Copy( B, BOrigProm );

    // Compute the initial guess
    // =========================
    applyAInv( B );
    DistMultiVec<PF> XProm(comm);
    Copy( B, XProm );

    Int refineIt = 0;
    DistMultiVec<PF> dXProm(comm), YProm(comm);
    Zeros( YProm, XProm.Height(), XProm.Width() );
    applyA( XProm, YProm );
    BProm -= YProm;

    const Int indent = PushIndent();
    while( true )
    {
        // Compute the proposed updates to the solutions
        // ---------------------------------------------
        Copy( BProm, B );
        applyAInv( B );
        Copy( B, dXProm );
        XProm += dXProm;

        ++refineIt;
        if( refineIt < maxRefineIts )
        {
            // Form the new residuals
            // ----------------------
            applyA( XProm, YProm );
            BProm = BOrigProm;
            BProm -= YProm;
        }
        else
            break;
    }
    SetIndent( indent );
    Copy( XProm, B );
    return refineIt;
}

} // namespace refined_solve

template<typename F,class ApplyAType,class ApplyAInvType>
inline Int PromotedRefinedSolve
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<F>& B,
        Base<F> relTol,
        Int maxRefineIts,
        bool progress )
{
    DEBUG_ONLY(CSE cse("PromotedRefinedSolve"))
    if( B.Width() == 1 )
        return refined_solve::PromotedSingle
               ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
    else
        return refined_solve::PromotedBatch
               ( applyA, applyAInv, B, maxRefineIts, progress );
}

} // namespace El

#endif // ifndef EL_SOLVE_REFINED_HPP
