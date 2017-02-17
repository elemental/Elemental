/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename F>
void MakeQuasiTriangular( UpperOrLower uplo, DistMatrix<F>& A )
{
    EL_DEBUG_ONLY(CallStackEntry cse("MakeQuasiTriangular"))
    const Int n = A.Height();
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
void TestQuasiTrsm
( bool print,
  LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation,
  Int m,
  Int n,
  F alpha,
  const Grid& g )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();

    DistMatrix<F> H(g), X(g);

    if( side == LEFT )
        HermitianUniformSpectrum( H, m, 1, 10 );
    else
        HermitianUniformSpectrum( H, n, 1, 10 );
    MakeQuasiTriangular( uplo, H );

    Uniform( X, m, n );
    DistMatrix<F> Y(g);
    if( side == LEFT )
        Gemm( orientation, NORMAL, F(1)/alpha, H, X, Y );
    else
        Gemm( NORMAL, orientation, F(1)/alpha, X, H, Y );

    if( print )
    {
        Print( H, "H" );
        Print( X, "X" );
        Print( Y, "Y" );
    }
    OutputFromRoot(g.Comm(),"Starting QuasiTrsm");

    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    QuasiTrsm( side, uplo, orientation, alpha, H, Y );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    const double realGFlops =
        ( side==LEFT ? double(m)*double(m)*double(n)
                     : double(m)*double(n)*double(n) ) /(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    OutputFromRoot
    (g.Comm(),"Finished in ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
        Print( Y, "Y after solve" );
    Y -= X;
    const auto HFrob = FrobeniusNorm( H );
    const auto XFrob = FrobeniusNorm( X );
    const auto EFrob = FrobeniusNorm( Y );
    OutputFromRoot
    (g.Comm(),
     "|| H ||_F = ",HFrob,"\n",Indent(),
     "|| X ||_F = ",XFrob,"\n",Indent(),
     "|| E ||_F = ",EFrob);

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
            ("--uplo","lower or upper quasi-triangular: L/U",'L');
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
        OutputFromRoot(comm,"Will test QuasiTrsm ",sideChar,uploChar,transChar);

        TestQuasiTrsm<float>
        ( print, side, uplo, orientation, m, n, float(3), g );
        TestQuasiTrsm<Complex<float>>
        ( print, side, uplo, orientation, m, n, Complex<float>(3), g );

        TestQuasiTrsm<double>
        ( print, side, uplo, orientation, m, n, double(3), g );
        TestQuasiTrsm<Complex<double>>
        ( print, side, uplo, orientation, m, n, Complex<double>(3), g );

#ifdef EL_HAVE_QUAD
        TestQuasiTrsm<Quad>
        ( print, side, uplo, orientation, m, n, Quad(3), g );
        TestQuasiTrsm<Complex<Quad>>
        ( print, side, uplo, orientation, m, n, Complex<Quad>(3), g );
#endif

#ifdef EL_HAVE_QD
        TestQuasiTrsm<DoubleDouble>
        ( print, side, uplo, orientation, m, n, DoubleDouble(3), g );
        TestQuasiTrsm<QuadDouble>
        ( print, side, uplo, orientation, m, n, QuadDouble(3), g );

        TestQuasiTrsm<Complex<DoubleDouble>>
        ( print, side, uplo, orientation, m, n, Complex<DoubleDouble>(3), g );
        TestQuasiTrsm<Complex<QuadDouble>>
        ( print, side, uplo, orientation, m, n, Complex<QuadDouble>(3), g );
#endif

#ifdef EL_HAVE_MPC
        TestQuasiTrsm<BigFloat>
        ( print, side, uplo, orientation, m, n, BigFloat(3), g );
        TestQuasiTrsm<Complex<BigFloat>>
        ( print, side, uplo, orientation, m, n, Complex<BigFloat>(3), g );

#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
