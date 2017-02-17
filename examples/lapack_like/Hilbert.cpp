/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

template<typename Real>
void SolveHilbert( El::Int n, El::Int maxRefineIts, bool refineProgress )
{
    El::Output("Attempting to solve Hilbert system with ",El::TypeName<Real>());

    // Form an n x n Hilbert matrix
    El::Matrix<Real> A;
    El::Hilbert( A, n );

    // Form a uniform random vector
    El::Matrix<Real> x;
    El::Uniform( x, n, 1 );
    const Real xFrob = El::FrobeniusNorm( x );

    // Form b := A x
    El::Matrix<Real> b;
    El::Gemv( El::NORMAL, Real(1), A, x, b );
    const Real bFrob = El::FrobeniusNorm( b );

    // Form xComp := inv(A) b. Rather than calling El::HPDSolve, we will instead
    // call El::hpd_solve::Overwrite so that we can reuse the Cholesky factor.
    El::Matrix<Real> xComp(b);
    El::Matrix<Real> L(A);
    El::hpd_solve::Overwrite( El::LOWER, El::NORMAL, L, xComp );

    // Form r := b - A x
    El::Matrix<Real> r(b);
    El::Gemv( El::NORMAL, Real(-1), A, xComp, Real(1), r );
    const Real rFrob = El::FrobeniusNorm( r );
    El::Output("|| r ||_2 / || b ||_2 = ",rFrob/bFrob);

    // Form e := x - xComp
    El::Matrix<Real> e(x);
    e -= xComp;
    const Real eFrob = El::FrobeniusNorm( e );
    El::Output("|| e ||_2 / || x ||_2 = ",eFrob/xFrob);

    // Run iterative refinement to try to increase the solution accuracy
    // (unfortunately, for Hilbert matrices the residual is already small)
    auto applyA = [&]( const El::Matrix<Real>& y, El::Matrix<Real>& Ay ) {
        El::Gemv( El::NORMAL, Real(1), A, y, Ay );
    };
    auto applyAInv = [&]( El::Matrix<Real>& y ) {
        El::cholesky::SolveAfter( El::LOWER, El::NORMAL, L, y );
    };
    const Real relTol = El::Pow( El::limits::Epsilon<Real>(), Real(0.9) );
    xComp = b;
    El::RefinedSolve
    ( applyA, applyAInv, xComp, relTol, maxRefineIts, refineProgress );

    // Form r := b - A x
    r = b;
    El::Gemv( El::NORMAL, Real(-1), A, xComp, Real(1), r );
    const Real rRefinedFrob = El::FrobeniusNorm( r );
    El::Output("|| rRefined ||_2 / || b ||_2 = ",rRefinedFrob/bFrob);

    // Form e := x - xComp
    e = x;
    e -= xComp;
    const Real eRefinedFrob = El::FrobeniusNorm( e );
    El::Output("|| eRefined ||_2 / || x ||_2 = ",eRefinedFrob/xFrob);

    // Run iterative refinement in increased precision to try to increase the
    // solution accuracy.
    typedef El::Promote<Real> RealProm;
    El::Matrix<RealProm> AProm;
    El::Copy( A, AProm );
    auto applyAProm = [&]( const El::Matrix<RealProm>& y,
                                 El::Matrix<RealProm>& Ay ) {
        El::Gemv( El::NORMAL, RealProm(1), AProm, y, Ay );
    };
    const Real relTolProm =
      Real(El::Pow( El::limits::Epsilon<RealProm>(), RealProm(0.9) ));
    El::Output
    ("Promoting from ",El::TypeName<Real>()," to ",El::TypeName<RealProm>(),
     " and setting relative tolerance to ",relTolProm);
    xComp = b;
    El::PromotedRefinedSolve
    ( applyAProm, applyAInv, xComp, relTolProm, maxRefineIts, refineProgress );

    // Form r := b - A x
    r = b;
    El::Gemv( El::NORMAL, Real(-1), A, xComp, Real(1), r );
    const Real rPromRefinedFrob = El::FrobeniusNorm( r );
    El::Output("|| rPromRefined ||_2 / || b ||_2 = ",rPromRefinedFrob/bFrob);

    // Form e := x - xComp
    e = x;
    e -= xComp;
    const Real ePromRefinedFrob = El::FrobeniusNorm( e );
    El::Output("|| ePromRefined ||_2 / || x ||_2 = ",ePromRefinedFrob/xFrob);

    El::Output("");
}

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );
    try
    {
        const El::Int n = El::Input("--n","matrix width",20);
        const El::Int maxRefineIts =
          El::Input("--maxRefineIts","max IR steps",3);
        const bool refineProgress =
          El::Input("--refineProgress","show IR progress?",false);
#ifdef EL_HAVE_MPC
        mpfr_prec_t prec = El::Input("--prec","MPFR precision",256);
#endif
        El::ProcessInput();

#ifdef EL_HAVE_MPC
        El::mpfr::SetPrecision( prec );
#endif

        try { SolveHilbert<float>( n, maxRefineIts, refineProgress ); }
        catch( std::exception& e ) { El::ReportException(e); }
        try { SolveHilbert<double>( n, maxRefineIts, refineProgress ); }
        catch( std::exception& e ) { El::ReportException(e); }
#ifdef EL_HAVE_QD
        try
        { SolveHilbert<El::DoubleDouble>( n, maxRefineIts, refineProgress ); }
        catch( std::exception& e ) { El::ReportException(e); }
        try { SolveHilbert<El::QuadDouble>( n, maxRefineIts, refineProgress ); }
        catch( std::exception& e ) { El::ReportException(e); }
#endif
#ifdef EL_HAVE_QUAD
        try { SolveHilbert<El::Quad>( n, maxRefineIts, refineProgress ); }
        catch( std::exception& e ) { El::ReportException(e); }
#endif
#ifdef EL_HAVE_MPC
        try { SolveHilbert<El::BigFloat>( n, maxRefineIts, refineProgress ); }
        catch( std::exception& e ) { El::ReportException(e); }
#endif
    } catch( std::exception& e ) { El::ReportException(e); }
    return 0;
}
