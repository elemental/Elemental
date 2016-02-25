/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template<typename F>
void MakeQuasiTriangular( UpperOrLower uplo, DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeQuasiTriangular"))
    const Int n = A.Height();
    if( n < 3 )
        return;
    if( uplo == LOWER )
    {
        MakeTrapezoidal( LOWER, A, 1 );
        auto dSup = GetDiagonal(A,+1);
        DistMatrix<F,STAR,STAR> dSup_STAR_STAR( dSup );
        for( Int j=0; j<n-2; ++j )
        {
            const F thisSup = dSup_STAR_STAR.Get(j,  0);
            const F nextSup = dSup_STAR_STAR.Get(j+1,0);
            if( thisSup != F(0) && nextSup != F(0) )
            {
                A.Set(j+1,j+2,0);
                dSup_STAR_STAR.Set(j+1,0,0);
            }
        }
    }
    else
    {
        MakeTrapezoidal( UPPER, A, -1 );
        auto dSub = GetDiagonal(A,-1);
        DistMatrix<F,STAR,STAR> dSub_STAR_STAR( dSub );
        for( Int j=0; j<n-2; ++j )
        {
            const F thisSub = dSub_STAR_STAR.Get(j,  0);
            const F nextSub = dSub_STAR_STAR.Get(j+1,0);
            if( thisSub != F(0) && nextSub != F(0) )
            {
                A.Set(j+2,j+1,0);
                dSub_STAR_STAR.Set(j+1,0,0);
            }
        }
    }
}

template<typename F> 
void TestMultiShiftQuasiTrsm
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation, 
  Int m,
  Int n,
  F alpha,
  const Grid& g,
  bool print )
{
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());
    typedef Base<F> Real;
    DistMatrix<F> H(g), X(g);
    DistMatrix<F,VR,STAR> shifts(g);

    if( side == LEFT )
    {
        HermitianUniformSpectrum( H, m, 1, 10 );
        Uniform( shifts, n, 1, F(0), Real(0.5) );
    }
    else
    {
        HermitianUniformSpectrum( H, n, 1, 10 );
        Uniform( shifts, m, 1, F(0), Real(0.5) );
    }
    MakeQuasiTriangular( uplo, H );

    auto modShifts(shifts);
    if( orientation == ADJOINT )
        Conjugate( modShifts );

    Uniform( X, m, n );
    DistMatrix<F> Y(g);
    Zeros( Y, m, n );
    if( side == LEFT )
    {
        Gemm( orientation, NORMAL, F(1)/alpha, H, X, F(1), Y );
        for( Int j=0; j<n; ++j )
        {
            auto x = X( ALL, IR(j) );
            auto y = Y( ALL, IR(j) );
            Axpy( -modShifts.Get(j,0)/alpha, x, y );
        }
    }
    else
    {
        Gemm( NORMAL, orientation, F(1)/alpha, X, H, F(1), Y );
        for( Int i=0; i<m; ++i )
        {
            auto x = X( IR(i), ALL );
            auto y = Y( IR(i), ALL );
            Axpy( -modShifts.Get(i,0)/alpha, x, y );
        }
    }

    if( print )
    {
        Print( H, "H" );
        Print( shifts, "shifts" );
        Print( X, "X" );
        Print( Y, "Y" );
    }
    if( g.Rank() == 0 )
        Output("  Starting MultiShiftQuasiTrsm");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    MultiShiftQuasiTrsm( side, uplo, orientation, alpha, H, shifts, Y );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 
        ( side==LEFT ? double(m)*double(m)*double(n)
                     : double(m)*double(n)*double(n) ) /(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
        Output("  Finished after ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
        Print( Y, "Y after solve" );
    Y -= X;
    const auto HFrob = FrobeniusNorm( H );
    const auto XFrob = FrobeniusNorm( X );
    const auto EFrob = FrobeniusNorm( Y );
    if( g.Rank() == 0 )
    {
        Output("  || H ||_F = ",HFrob);
        Output("  || X ||_F = ",XFrob);
        Output("  || E ||_F = ",EFrob);
    }
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--r","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char sideChar = Input("--side","side to solve from: L/R",'L');
        const char uploChar = Input
            ("--uplo","lower or upper quasi-triangular: L/U",'L');
        const char transChar = Input
            ("--trans","orientation of quasi-triangular matrix: N/T/C",'N');
        const Int m = Input("--m","height of result",100);
        const Int n = Input("--n","width of result",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const LeftOrRight side = CharToLeftOrRight( sideChar );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        if( commRank == 0 )
            Output("Testing MultiShiftQuasiTrsm ",sideChar,uploChar,transChar);

        TestMultiShiftQuasiTrsm<float>
        ( side, uplo, orientation,
          m, n,
          float(3),
          g, print );
        TestMultiShiftQuasiTrsm<Complex<float>>
        ( side, uplo, orientation,
          m, n,
          Complex<float>(3),
          g, print );

        TestMultiShiftQuasiTrsm<double>
        ( side, uplo, orientation,
          m, n,
          double(3),
          g, print );
        TestMultiShiftQuasiTrsm<Complex<double>>
        ( side, uplo, orientation,
          m, n,
          Complex<double>(3),
          g, print );

#ifdef EL_HAVE_QD
        TestMultiShiftQuasiTrsm<DoubleDouble>
        ( side, uplo, orientation,
          m, n,
          DoubleDouble(3),
          g, print );
        TestMultiShiftQuasiTrsm<QuadDouble>
        ( side, uplo, orientation,
          m, n,
          QuadDouble(3),
          g, print );
#endif

#ifdef EL_HAVE_QUAD
        TestMultiShiftQuasiTrsm<Quad>
        ( side, uplo, orientation,
          m, n,
          Quad(3),
          g, print );
        TestMultiShiftQuasiTrsm<Complex<Quad>>
        ( side, uplo, orientation,
          m, n,
          Complex<Quad>(3),
          g, print );
#endif

#ifdef EL_HAVE_MPC
        TestMultiShiftQuasiTrsm<BigFloat>
        ( side, uplo, orientation,
          m, n,
          BigFloat(3),
          g, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
