/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SOLVE_REFINED_HPP
#define EL_SOLVE_REFINED_HPP

// The following routines each have three underlying implementations:
//
//  1) a single right-hand side,
//  2) a pair of right-hand sides, and
//  3) an arbitrary number of right-hand sides.
//
// In the first and second cases, proper iterative refinement is performed.
// In the third case, "batch" iterative refinement is performed for a
// specified number of iterations. While batch refinement will occasionally
// result in higher residual norms than necessary (as monotonicity must be
// specifically enforced), and more iterations than necessary may be performed,
// it is a compromise between performance and accuracy.

// TODO(poulson): DistMatrix implementations
// TODO(poulson): Simplify once DistMultiVec is eliminated
// TODO(poulson): Allow for a choice between max and two norms?

namespace El {

template<typename Real>
struct RefinedSolveInfo
{
    Int numIts=0;
    Real relTol=1;
    bool metRequestedTol=false;
};

namespace refined_solve {

// In what follows, 'applyA' should be a function of the form
//
//   void applyA( const Matrix<Field>& x, Matrix<Field>& y )
//
// and overwrite y with A x. However, 'applyAInv' should have the form
//
//   void applyAInv( Matrix<Field>& b )
//
// and overwrite b with inv(A) b.
//

template<typename Field,class ApplyAType,class ApplyAInvType>
RefinedSolveInfo<Base<Field>>
Single
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<Field>& b,
        Base<Field> relTol,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    typedef Base<Field> Real;
    RefinedSolveInfo<Real> info;

    auto bOrig = b;
    const Real bNorm = MaxNorm( b );

    // Compute the initial guess
    // =========================
    auto x = b;
    applyAInv( x );

    Matrix<Field> dx, xCand, y;
    Zeros( y, x.Height(), 1 );

    applyA( x, y );
    b -= y;
    Real errorNorm = MaxNorm( b );
    info.relTol = errorNorm / bNorm;
    if( progress )
        Output("original relative accuracy: ",info.relTol);

    while( true )
    {
        if( info.relTol <= relTol )
        {
            if( progress )
                Output(info.relTol," <= ",relTol);
            info.metRequestedTol = true;
            break;
        }
        if( info.numIts >= maxRefineIts )
        {
            if( progress )
                Output("WARNING: Iterative refinement did not converge");
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
            Output("refined relative accuracy: ",newErrorNorm/bNorm);

        ++info.numIts;
        if( newErrorNorm < errorNorm )
        {
            x = xCand;
            errorNorm = newErrorNorm;
            info.relTol = newErrorNorm/bNorm;
        }
        else
        {
            if( progress )
                Output("WARNING: Not accepting lower quality iterate");
            break;
        }
    }
    b = x;
    return info;
}

template<typename Field,class ApplyAType,class ApplyAInvType>
RefinedSolveInfo<Base<Field>>
Pair
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<Field>& B,
        Base<Field> relTol,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( B.Width() != 2 )
          LogicError("Expected a pair of right-hand sides");
    )
    typedef Base<Field> Real;
    RefinedSolveInfo<Real> info;

    auto bL = B( ALL, IR(0) );
    auto bR = B( ALL, IR(1) );
    const Real bLNorm = MaxNorm( bL );
    const Real bRNorm = MaxNorm( bR );

    auto BOrig = B;
    auto bOrigL = BOrig( ALL, IR(0) );
    auto bOrigR = BOrig( ALL, IR(1) );

    // Compute the initial guess
    // =========================
    auto X = B;
    auto xL = X( ALL, IR(0) );
    auto xR = X( ALL, IR(1) );
    applyAInv( X );

    Matrix<Field> dX, XCand, Y;
    Zeros( dX, X.Height(), X.Width() );
    Zeros( XCand, X.Height(), X.Width() );
    Zeros( Y, X.Height(), X.Width() );

    auto dxL = dX( ALL, IR(0) );
    auto dxR = dX( ALL, IR(1) );
    auto xLCand = XCand( ALL, IR(0) );
    auto xRCand = XCand( ALL, IR(1) );
    auto yL = Y( ALL, IR(0) );
    auto yR = Y( ALL, IR(1) );

    applyA( X, Y );
    B -= Y;
    Real errorNormL = MaxNorm( bL );
    Real errorNormR = MaxNorm( bR );
    if( progress )
        Output
        ("original rel errors: ",errorNormL/bLNorm," and ",errorNormR/bRNorm);

    bool leftConv=false, rightConv=false;
    bool leftMetTol=false, rightMetTol=false;
    while( true )
    {
        const Real relErrorL = errorNormL/bLNorm;
        const Real relErrorR = errorNormR/bRNorm;
        info.relTol = Max( relErrorL, relErrorR );
        if( !leftConv && relErrorL <= relTol )
        {
            if( progress )
                Output("Left converged with ",relErrorL," <= ",relTol);
            leftConv = true;
            leftMetTol = true;
        }
        if( !rightConv && relErrorR <= relTol )
        {
            if( progress )
                Output("Right converged with ",relErrorR," <= ",relTol);
            rightConv = true;
            rightMetTol = true;
        }
        if( leftConv && rightConv )
        {
            info.metRequestedTol = leftMetTol && rightMetTol;
            break;
        }
        if( info.numIts >= maxRefineIts )
        {
            if( progress )
                Output("WARNING: Iterative refinement did not converge");
            break;
        }

        if( leftConv )
        {
            // Compute the proposed update to the solution
            // -------------------------------------------
            dxR = bR;
            applyAInv( dxR );
            xRCand = xR;
            xRCand += dxR;

            // Check the new residual
            // ----------------------
            applyA( xRCand, yR );
            bR = bOrigR;
            bR -= yR;
            auto newErrorNormR = MaxNorm( bR );
            if( progress )
                Output("Right refined rel error: ",newErrorNormR/bRNorm);

            ++info.numIts;
            if( newErrorNormR < errorNormR )
            {
                xR = xRCand;
                errorNormR = newErrorNormR;
            }
            else
            {
                if( progress )
                    Output
                    ("WARNING: Not accepting lower quality right iterate");
                rightConv = true;
            }
        }
        else if( rightConv )
        {
            // Compute the proposed update to the solution
            // -------------------------------------------
            dxL = bL;
            applyAInv( dxL );
            xLCand = xL;
            xLCand += dxL;

            // Check the new residual
            // ----------------------
            applyA( xLCand, yL );
            bL = bOrigL;
            bL -= yL;
            auto newErrorNormL = MaxNorm( bL );
            if( progress )
                Output("Left refined rel error: ",newErrorNormL/bLNorm);

            ++info.numIts;
            if( newErrorNormL < errorNormL )
            {
                xL = xLCand;
                errorNormL = newErrorNormL;
            }
            else
            {
                if( progress )
                    Output("WARNING: Not accepting lower quality left iterate");
                leftConv = true;
            }
        }
        else
        {
            // Compute the proposed update to the solution
            // -------------------------------------------
            dX = B;
            applyAInv( dX );
            XCand = X;
            XCand += dX;

            // Check the new residual
            // ----------------------
            applyA( XCand, Y );
            B = BOrig;
            B -= Y;
            auto newErrorNormL = MaxNorm( bL );
            auto newErrorNormR = MaxNorm( bR );
            if( progress )
            {
                Output("Left refined rel error: ",newErrorNormL/bLNorm);
                Output("Right refined rel error: ",newErrorNormR/bRNorm);
            }

            ++info.numIts;
            if( newErrorNormL < errorNormL )
            {
                xL = xLCand;
                errorNormL = newErrorNormL;
            }
            else
            {
                if( progress )
                    Output("WARNING: Not accepting lower quality left iterate");
                leftConv = true;
            }

            if( newErrorNormR < errorNormR )
            {
                xR = xRCand;
                errorNormR = newErrorNormR;
            }
            else
            {
                if( progress )
                    Output
                    ("WARNING: Not accepting lower quality right iterate");
                rightConv = true;
            }
        }
    }
    B = X;
    return info;
}

template<typename Field,class ApplyAType,class ApplyAInvType>
RefinedSolveInfo<Base<Field>>
Batch
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<Field>& B,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    RefinedSolveInfo<Real> info;

    // TODO(poulson): Allow for early exits

    // Compute the initial guesses
    // ===========================
    auto BOrig = B;
    auto X = B;
    applyAInv( X );

    Matrix<Field> dX, Y;
    Zeros( Y, X.Height(), X.Width() );
    applyA( X, Y );
    B -= Y;

    // TODO(poulson): Compute original column norms of B.

    while( true )
    {
        // Compute the updates to the solutions
        // ------------------------------------
        dX = B;
        applyAInv( dX );
        X += dX;

        ++info.numIts;
        if( info.numIts < maxRefineIts )
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

    // TODO(poulson): Properly check convergence tolerances.
    info.relTol = 0;
    info.metRequestedTol = true;

    return info;
}

} // namespace refined_solve

template<typename Field,class ApplyAType,class ApplyAInvType>
RefinedSolveInfo<Base<Field>>
RefinedSolve
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<Field>& B,
        Base<Field> relTol,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    if( B.Width() == 1 )
        return refined_solve::Single
               ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
    else if( B.Width() == 2 )
        return refined_solve::Pair
               ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
    else
        return refined_solve::Batch
               ( applyA, applyAInv, B, maxRefineIts, progress );
}

namespace refined_solve {

template<typename Field,class ApplyAType,class ApplyAInvType>
RefinedSolveInfo<Base<Field>>
PromotedSingle
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<Field>& b,
        Base<Field> relTol,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    typedef Base<Field> Real;
    typedef Promote<Real> PReal;
    typedef Promote<Field> PField;
    RefinedSolveInfo<Real> info;

    Matrix<PField> bProm, bOrigProm;
    Copy( b, bProm );
    Copy( b, bOrigProm );
    const PReal bNorm = MaxNorm( bOrigProm );

    // Compute the initial guess
    // =========================
    applyAInv( b );
    Matrix<PField> xProm;
    Copy( b, xProm );

    Matrix<PField> dxProm, xCandProm, yProm;
    Zeros( yProm, xProm.Height(), 1 );

    applyA( xProm, yProm );
    bProm -= yProm;
    auto errorNorm = MaxNorm( bProm );
    info.relTol = Real(errorNorm/bNorm);
    if( progress )
        Output("original rel error: ",info.relTol);

    while( true )
    {
        if( info.relTol <= relTol )
        {
            if( progress )
                Output(info.relTol," <= ",relTol);
            info.metRequestedTol = true;
            break;
        }
        if( info.numIts >= maxRefineIts )
        {
            if( progress )
                Output("WARNING: Iterative refinement did not converge");
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

        ++info.numIts;
        if( newErrorNorm < errorNorm )
        {
            xProm = xCandProm;
            errorNorm = newErrorNorm;
            info.relTol = Real(newErrorNorm/bNorm);
        }
        else
        {
            if( progress )
                Output("WARNING: Not accepting lower quality iterate");
            break;
        }
    }
    Copy( xProm, b );
    return info;
}

template<typename Field,class ApplyAType,class ApplyAInvType>
RefinedSolveInfo<Base<Field>>
PromotedPair
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<Field>& B,
        Base<Field> relTol,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( B.Width() != 2 )
          LogicError("Expected a pair of right-hand sides");
    )
    typedef Base<Field> Real;
    typedef Promote<Real> PReal;
    typedef Promote<Field> PField;
    RefinedSolveInfo<Real> info;

    Matrix<PField> BProm, BOrigProm;
    Copy( B, BProm );
    Copy( B, BOrigProm );

    auto bL = B( ALL, IR(0) );
    auto bR = B( ALL, IR(1) );
    auto bLProm = BProm( ALL, IR(0) );
    auto bRProm = BProm( ALL, IR(1) );
    auto bLOrigProm = BOrigProm( ALL, IR(0) );
    auto bROrigProm = BOrigProm( ALL, IR(1) );

    const PReal bLNorm = MaxNorm( bLOrigProm );
    const PReal bRNorm = MaxNorm( bROrigProm );

    // Compute the initial guess
    // =========================
    applyAInv( B );
    Matrix<PField> XProm;
    Copy( B, XProm );

    auto xLProm = XProm( ALL, IR(0) );
    auto xRProm = XProm( ALL, IR(1) );

    Matrix<PField> dXProm, XCandProm, YProm;
    Zeros( dXProm, XProm.Height(), XProm.Width() );
    Zeros( XCandProm, XProm.Height(), XProm.Width() );
    Zeros( YProm, XProm.Height(), XProm.Width() );

    auto dxLProm = dXProm( ALL, IR(0) );
    auto dxRProm = dXProm( ALL, IR(1) );
    auto xLCandProm = XCandProm( ALL, IR(0) );
    auto xRCandProm = XCandProm( ALL, IR(1) );
    auto yLProm = YProm( ALL, IR(0) );
    auto yRProm = YProm( ALL, IR(1) );

    applyA( XProm, YProm );
    BProm -= YProm;
    auto errorNormL = MaxNorm( bLProm );
    auto errorNormR = MaxNorm( bRProm );
    if( progress )
        Output
        ("original rel errors: ",errorNormL/bLNorm," and ",errorNormR/bRNorm);

    bool leftConv=false, rightConv=false;
    bool leftMetTol=false, rightMetTol=false;
    while( true )
    {
        const PReal relErrorL = errorNormL/bLNorm;
        const PReal relErrorR = errorNormR/bRNorm;
        info.relTol = Real( Max( relErrorL, relErrorR ) );
        if( !leftConv && relErrorL <= relTol )
        {
            if( progress )
                Output("Left converged with ",relErrorL," <= ",relTol);
            leftConv = true;
            leftMetTol = true;
        }
        if( !rightConv && relErrorR <= relTol )
        {
            if( progress )
                Output("Right converged with ",relErrorR," <= ",relTol);
            rightConv = true;
            rightMetTol = true;
        }
        if( leftConv && rightConv )
        {
            info.metRequestedTol = leftMetTol && rightMetTol;
            break;
        }
        if( info.numIts >= maxRefineIts )
        {
            if( progress )
                Output("WARNING: Iterative refinement did not converge");
            break;
        }

        if( leftConv )
        {
            // Compute the proposed update to the solution
            // -------------------------------------------
            Copy( bRProm, bR );
            applyAInv( bR );
            Copy( bR, dxRProm );
            xRCandProm = xRProm;
            xRCandProm += dxRProm;

            // Check the new residual
            // ----------------------
            applyA( xRCandProm, yRProm );
            bRProm = bROrigProm;
            bRProm -= yRProm;
            auto newErrorNormR = MaxNorm( bRProm );
            if( progress )
                Output("Right refined rel error: ",newErrorNormR/bRNorm);

            ++info.numIts;
            if( newErrorNormR < errorNormR )
            {
                xRProm = xRCandProm;
                errorNormR = newErrorNormR;
            }
            else
            {
                if( progress )
                    Output
                    ("WARNING: Not accepting lower quality right iterate");
                rightConv = true;
            }
        }
        else if( rightConv )
        {
            // Compute the proposed update to the solution
            // -------------------------------------------
            Copy( bLProm, bL );
            applyAInv( bL );
            Copy( bL, dxLProm );
            xLCandProm = xLProm;
            xLCandProm += dxLProm;

            // Check the new residual
            // ----------------------
            applyA( xLCandProm, yLProm );
            bLProm = bLOrigProm;
            bLProm -= yLProm;
            auto newErrorNormL = MaxNorm( bLProm );
            if( progress )
                Output("Left refined rel error: ",newErrorNormL/bLNorm);

            ++info.numIts;
            if( newErrorNormL < errorNormL )
            {
                xLProm = xLCandProm;
                errorNormL = newErrorNormL;
            }
            else
            {
                if( progress )
                    Output("WARNING: Not accepting lower quality left iterate");
                leftConv = true;
            }
        }
        else
        {
            // Compute the proposed update to the solution
            // -------------------------------------------
            Copy( BProm, B );
            applyAInv( B );
            Copy( B, dXProm );
            XCandProm = XProm;
            XCandProm += dXProm;

            // Check the new residual
            // ----------------------
            applyA( XCandProm, YProm );
            BProm = BOrigProm;
            BProm -= YProm;
            auto newErrorNormL = MaxNorm( bLProm );
            auto newErrorNormR = MaxNorm( bRProm );
            if( progress )
            {
                Output("Left refined rel error: ",newErrorNormL/bLNorm);
                Output("Right refined rel error: ",newErrorNormR/bRNorm);
            }

            ++info.numIts;
            if( newErrorNormL < errorNormL )
            {
                xLProm = xLCandProm;
                errorNormL = newErrorNormL;
            }
            else
            {
                if( progress )
                    Output("WARNING: Not accepting lower quality left iterate");
                leftConv = true;
            }

            if( newErrorNormR < errorNormR )
            {
                xRProm = xRCandProm;
                errorNormR = newErrorNormR;
            }
            else
            {
                if( progress )
                    Output
                    ("WARNING: Not accepting lower quality right iterate");
                rightConv = true;
            }
        }
    }
    Copy( XProm, B );
    return info;
}

template<typename Field,class ApplyAType,class ApplyAInvType>
RefinedSolveInfo<Base<Field>>
PromotedBatch
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<Field>& B,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    typedef Promote<Field> PField;
    typedef Base<Field> Real;
    RefinedSolveInfo<Real> info;

    Matrix<PField> BProm, BOrigProm;
    Copy( B, BProm );
    Copy( B, BOrigProm );

    // Compute the initial guesses
    // ===========================
    applyAInv( B );
    Matrix<PField> XProm;
    Copy( B, XProm );

    Matrix<PField> dXProm, YProm;
    Zeros( YProm, XProm.Height(), XProm.Width() );
    applyA( XProm, YProm );
    BProm -= YProm;

    // TODO(poulson): Compute original error norms.

    while( true )
    {
        // Update the solutions
        // --------------------
        Copy( BProm, B );
        applyAInv( B );
        Copy( B, dXProm );
        XProm += dXProm;

        ++info.numIts;
        if( info.numIts < maxRefineIts )
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

    // TODO(poulson): Properly compute error norms.
    info.relTol = 0;
    info.metRequestedTol = true;

    return info;
}

} // namespace refined_solve

template<typename Field,class ApplyAType,class ApplyAInvType>
DisableIf<IsSame<Field,Promote<Field>>,RefinedSolveInfo<Base<Field>>>
PromotedRefinedSolve
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<Field>& B,
        Base<Field> relTol,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    if( B.Width() == 1 )
        return refined_solve::PromotedSingle
               ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
    else if( B.Width() == 2 )
        return refined_solve::PromotedPair
               ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
    else
        return refined_solve::PromotedBatch
               ( applyA, applyAInv, B, maxRefineIts, progress );
}

template<typename Field,class ApplyAType,class ApplyAInvType>
EnableIf<IsSame<Field,Promote<Field>>,RefinedSolveInfo<Base<Field>>>
PromotedRefinedSolve
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        Matrix<Field>& B,
        Base<Field> relTol,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    return RefinedSolve( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

namespace refined_solve {

template<typename Field,class ApplyAType,class ApplyAInvType>
RefinedSolveInfo<Base<Field>>
Single
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<Field>& b,
        Base<Field> relTol,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    typedef Base<Field> Real;
    const Grid& grid = b.Grid();
    const int commRank = grid.Rank();
    RefinedSolveInfo<Real> info;

    auto bOrig = b;
    const Real bNorm = MaxNorm( b );

    // Compute the initial guess
    // =========================
    auto x = b;
    applyAInv( x );

    DistMultiVec<Field> dx(grid), xCand(grid), y(grid);
    Zeros( y, x.Height(), 1 );
    applyA( x, y );
    b -= y;
    Real errorNorm = MaxNorm( b );
    info.relTol = errorNorm/bNorm;
    if( progress && commRank == 0 )
        Output("original rel error: ",info.relTol);

    const Int indent = PushIndent();
    while( true )
    {
        if( info.relTol <= relTol )
        {
            if( progress && commRank == 0 )
                Output(info.relTol," <= ",relTol);
            info.metRequestedTol = true;
            break;
        }
        if( info.numIts >= maxRefineIts )
        {
            if( progress && commRank == 0 )
                Output("WARNING: Iterative refinement did not converge");
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
        Real newErrorNorm = MaxNorm( b );
        if( progress && commRank == 0 )
            Output("refined rel error: ",newErrorNorm/bNorm);

        ++info.numIts;
        if( newErrorNorm < errorNorm )
        {
            x = xCand;
            errorNorm = newErrorNorm;
            info.relTol = errorNorm/bNorm;
        }
        else
        {
            if( progress && commRank == 0 )
                Output("WARNING: Not accepting lower quality iterate");
            break;
        }
    }
    SetIndent( indent );
    b = x;
    return info;
}

template<typename Field,class ApplyAType,class ApplyAInvType>
RefinedSolveInfo<Base<Field>>
Batch
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<Field>& B,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Grid& grid = B.Grid();
    RefinedSolveInfo<Real> info;

    auto BOrig = B;

    // Compute the initial guess
    // =========================
    auto X = B;
    applyAInv( X );

    DistMultiVec<Field> dX(grid), Y(grid);
    Zeros( Y, X.Height(), X.Width() );
    applyA( X, Y );
    B -= Y;

    // TODO(poulson): Compute original error norms.

    const Int indent = PushIndent();
    while( true )
    {
        // Compute the proposed updates to the solutions
        // ---------------------------------------------
        dX = B;
        applyAInv( dX );
        X += dX;

        ++info.numIts;
        if( info.numIts < maxRefineIts )
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

    // TODO(poulson): Properly compute errors.
    info.relTol = 0;
    info.metRequestedTol = true;

    return info;
}

} // namespace refined_solve

template<typename Field,class ApplyAType,class ApplyAInvType>
RefinedSolveInfo<Base<Field>>
RefinedSolve
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<Field>& B,
        Base<Field> relTol,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    if( B.Width() == 1 )
        return refined_solve::Single
               ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
    else
        return refined_solve::Batch
               ( applyA, applyAInv, B, maxRefineIts, progress );
}

namespace refined_solve {

template<typename Field,class ApplyAType,class ApplyAInvType>
RefinedSolveInfo<Base<Field>>
PromotedSingle
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<Field>& b,
        Base<Field> relTol,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    typedef Base<Field> Real;
    typedef Promote<Field> PField;
    const Grid& grid = b.Grid();
    const int commRank = grid.Rank();
    RefinedSolveInfo<Real> info;

    DistMultiVec<PField> bProm(grid), bOrigProm(grid);
    Copy( b, bProm );
    Copy( b, bOrigProm );
    const auto bNorm = MaxNorm( bProm );

    // Compute the initial guess
    // =========================
    applyAInv( b );
    DistMultiVec<PField> xProm(grid);
    Copy( b, xProm );

    DistMultiVec<PField> dxProm(grid), xCandProm(grid), yProm(grid);
    Zeros( yProm, xProm.Height(), 1 );
    applyA( xProm, yProm );
    bProm -= yProm;
    auto errorNorm = MaxNorm( bProm );
    info.relTol = Real(errorNorm/bNorm);
    if( progress && commRank == 0 )
        Output("original rel error: ",info.relTol);

    const Int indent = PushIndent();
    while( true )
    {
        if( info.relTol <= relTol )
        {
            if( progress && commRank == 0 )
                Output(info.relTol," <= ",relTol);
            info.metRequestedTol = true;
            break;
        }
        if( info.numIts >= maxRefineIts )
        {
            if( progress && commRank == 0 )
                Output("WARNING: Iterative refinement did not converge");
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
        if( progress && commRank == 0 )
            Output("refined rel error: ",newErrorNorm/bNorm);

        ++info.numIts;
        if( newErrorNorm < errorNorm )
        {
            xProm = xCandProm;
            errorNorm = newErrorNorm;
            info.relTol = Real(errorNorm/bNorm);
        }
        else
        {
            if( progress && commRank == 0 )
                Output("WARNING: Not accepting lower quality iterate");
            break;
        }
    }
    SetIndent( indent );
    Copy( xProm, b );
    return info;
}

template<typename Field,class ApplyAType,class ApplyAInvType>
RefinedSolveInfo<Base<Field>>
PromotedBatch
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<Field>& B,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Promote<Field> PField;
    const Grid& grid = B.Grid();
    RefinedSolveInfo<Real> info;

    DistMultiVec<PField> BProm(grid), BOrigProm(grid);
    Copy( B, BProm );
    Copy( B, BOrigProm );

    // Compute the initial guess
    // =========================
    applyAInv( B );
    DistMultiVec<PField> XProm(grid);
    Copy( B, XProm );

    DistMultiVec<PField> dXProm(grid), YProm(grid);
    Zeros( YProm, XProm.Height(), XProm.Width() );
    applyA( XProm, YProm );
    BProm -= YProm;

    // TODO(poulson): Compute original error norms.

    const Int indent = PushIndent();
    while( true )
    {
        // Compute the proposed updates to the solutions
        // ---------------------------------------------
        Copy( BProm, B );
        applyAInv( B );
        Copy( B, dXProm );
        XProm += dXProm;

        ++info.numIts;
        if( info.numIts < maxRefineIts )
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

    // TODO(poulson): Properly compute error norms.
    info.relTol = 0;
    info.metRequestedTol = true;

    return info;
}

} // namespace refined_solve

template<typename Field,class ApplyAType,class ApplyAInvType>
DisableIf<IsSame<Field,Promote<Field>>,RefinedSolveInfo<Base<Field>>>
PromotedRefinedSolve
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<Field>& B,
        Base<Field> relTol,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    if( B.Width() == 1 )
        return refined_solve::PromotedSingle
               ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
    else
        return refined_solve::PromotedBatch
               ( applyA, applyAInv, B, maxRefineIts, progress );
}

template<typename Field,class ApplyAType,class ApplyAInvType>
EnableIf<IsSame<Field,Promote<Field>>,RefinedSolveInfo<Base<Field>>>
PromotedRefinedSolve
( const ApplyAType& applyA,
  const ApplyAInvType& applyAInv,
        DistMultiVec<Field>& B,
        Base<Field> relTol,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    return RefinedSolve( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

} // namespace El

#endif // ifndef EL_SOLVE_REFINED_HPP
