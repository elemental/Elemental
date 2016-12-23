/*
   Copyright (c) 2009-2016, Jack Poulson and Tim Moon
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename F>
void TestSafeMultiShiftTrsm
( bool print,
  LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation,
  Int m, Int n,
  F alpha,
  const Grid& g )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();

    typedef Base<F> Real;
    DistMatrix<F> A(g), B(g), X(g);
    DistMatrix<F,VR,STAR> shifts(g), scales(g);

    // Generate test matrices
    if( side == LEFT )
    {
        HermitianUniformSpectrum( A, m, -5, 5 );
        Uniform( shifts, n, 1, F(0), Real(10) );
    }
    else
    {
        HermitianUniformSpectrum( A, n, -5, 5 );
        Uniform( shifts, m, 1, F(0), Real(1) );
    }
    MakeTrapezoidal( uplo, A );
    Uniform( B, m, n );
    if( print )
    {
        Print( A, "A" );
        Print( shifts, "shifts" );
        OutputFromRoot(g.Comm(),"alpha=",alpha);
        Print( B, "B" );
    }

    auto modShifts( shifts );
    if( orientation == ADJOINT )
        Conjugate( modShifts );

    // Perform SafeMultiShiftTrsm
    X = B;
    OutputFromRoot(g.Comm(),"Starting SafeMultiShiftTrsm");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    SafeMultiShiftTrsm
    ( side, uplo, orientation,
      alpha, A, shifts, X, scales );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    OutputFromRoot(g.Comm(),"Finished after ",runTime," seconds");
    if( print )
    {
        Print( X, "X" );
        Print( scales, "scales" );
    }

    // Compute residual
    DistMatrix<F> R(B);
    if( side == LEFT )
    {
        DiagonalScale( RIGHT, NORMAL, scales, R );
        Gemm( orientation, NORMAL, F(1), A, X, -alpha, R );
        for( Int j=0; j<n; ++j )
        {
            auto xj = X( ALL, IR(j) );
            auto rj = R( ALL, IR(j) );
            Axpy( -modShifts.Get(j,0), xj, rj );
        }
    }
    else
    {
        DiagonalScale( LEFT, NORMAL, scales, R );
        Gemm( NORMAL, orientation, F(1), X, A, -alpha, R );
        for( Int i=0; i<m; ++i )
        {
            auto xi = X( IR(i), ALL );
            auto ri = R( IR(i), ALL );
            Axpy( -modShifts.Get(i,0), xi, ri );
        }
    }

    // Compute maximum relative error
    const DistMatrix<F> diag = GetDiagonal( A );
    DistMatrix<F> shiftedDiag( diag );
    FillDiagonal( A, F(0) );
    const Real AOffFrob = FrobeniusNorm( A );
    SetDiagonal( A, diag );
    Real maxErr = 0;
    Real maxRelErr = 0;
    Real minScales = 1;
    Real maxNorm, minNorm;
    if( side == LEFT )
    {
        for( Int j=0; j<n; ++j )
        {
            const Real RjNrm2 = Nrm2( R(ALL,IR(j)) );
            const Real XjNrm2 = Nrm2( X(ALL,IR(j)) );
            shiftedDiag = diag;
            for( Int i=0; i<m; ++i )
            {
              shiftedDiag.Update( i, 0, shifts.Get(j,0) );
            }
            const Real AjDiagFrob = Nrm2( shiftedDiag );
            const Real AjFrob =
              Sqrt( AOffFrob*AOffFrob + AjDiagFrob*AjDiagFrob );
            const Real currErr = RjNrm2/(AjFrob*XjNrm2);
            maxErr = Max( maxErr, RjNrm2 );
            maxRelErr = Max( maxRelErr, currErr );
            minScales = Min( minScales, Abs(scales.Get(j,0)) );
            minNorm = j>0 ? Min( minNorm, XjNrm2 ) : XjNrm2;
            maxNorm = j>0 ? Max( maxNorm, XjNrm2 ) : XjNrm2;
        }
    }
    else
        LogicError("side=RIGHT is not yet supported");

    const Real AFrob = FrobeniusNorm( A );
    OutputFromRoot
    (g.Comm(),
     "|| A ||_F = ",AFrob,"\n",Indent(),
     "max( || Aj xj - bj sj ||_2 ) = ",maxErr,"\n",Indent(),
     "max( || Aj xj - bj sj ||_2 / || Aj ||_F || xj ||_2 ) = ",maxRelErr,
     "\n",Indent(),
     "max( || xj ||_2 ) = ",maxNorm,"\n",Indent(),
     "min( || xj ||_2 ) = ",minNorm,"\n",Indent(),
     "min( sj ) = ",minScales);

    PopIndent();
}

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        int gridHeight = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char sideChar = Input("--side","side to solve from: L/R",'L');
        const char uploChar = Input
            ("--uplo","lower or upper quasi-triangular: L/U",'U');
        const char transChar = Input
            ("--trans","orientation of quasi-triangular matrix: N/T/C",'N');
        const Int m = Input("--m","height of result",100);
        const Int n = Input("--n","width of result",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( gridHeight == 0 )
            gridHeight = Grid::DefaultHeight( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );
        const LeftOrRight side = CharToLeftOrRight( sideChar );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        OutputFromRoot
        (comm,"Will test SafeMultiShiftTrsm ",sideChar,uploChar,transChar);

        TestSafeMultiShiftTrsm<float>
        ( print, side, uplo, orientation, m, n, float(3), g );
        TestSafeMultiShiftTrsm<Complex<float>>
        ( print, side, uplo, orientation, m, n, Complex<float>(3), g );

        TestSafeMultiShiftTrsm<double>
        ( print, side, uplo, orientation, m, n, double(3), g );
        TestSafeMultiShiftTrsm<Complex<double>>
        ( print, side, uplo, orientation, m, n, Complex<double>(3), g );

#ifdef EL_HAVE_QUAD
        TestSafeMultiShiftTrsm<Quad>
        ( print, side, uplo, orientation, m, n, Quad(3), g );
        TestSafeMultiShiftTrsm<Complex<Quad>>
        ( print, side, uplo, orientation, m, n, Complex<Quad>(3), g );
#endif

#ifdef EL_HAVE_QD
        TestSafeMultiShiftTrsm<DoubleDouble>
        ( print, side, uplo, orientation, m, n, DoubleDouble(3), g );
        TestSafeMultiShiftTrsm<QuadDouble>
        ( print, side, uplo, orientation, m, n, QuadDouble(3), g );

        TestSafeMultiShiftTrsm<Complex<DoubleDouble>>
        ( print, side, uplo, orientation, m, n, Complex<DoubleDouble>(3), g );
        TestSafeMultiShiftTrsm<Complex<QuadDouble>>
        ( print, side, uplo, orientation, m, n, Complex<QuadDouble>(3), g );
#endif

#ifdef EL_HAVE_MPC
        TestSafeMultiShiftTrsm<BigFloat>
        ( print, side, uplo, orientation, m, n, BigFloat(3), g );
        TestSafeMultiShiftTrsm<Complex<BigFloat>>
        ( print, side, uplo, orientation, m, n, Complex<BigFloat>(3), g );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
