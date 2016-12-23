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
void TestTrsm
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation,
  UnitOrNonUnit diag,
  Int m,
  Int n,
  F alpha,
  const Grid& g,
  bool print )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();

    DistMatrix<F> A(g), X(g);

    if( side == LEFT )
        HermitianUniformSpectrum( A, m, 1, 10 );
    else
        HermitianUniformSpectrum( A, n, 1, 10 );
    auto S( A );
    MakeTrapezoidal( uplo, S );
    if( diag == UNIT )
        FillDiagonal( S, F(1) );

    Uniform( X, m, n );
    DistMatrix<F> Y(g);
    if( side == LEFT )
        Gemm( orientation, NORMAL, F(1)/alpha, S, X, Y );
    else
        Gemm( NORMAL, orientation, F(1)/alpha, X, S, Y );

    if( print )
    {
        Print( A, "A" );
        Print( S, "S" );
        Print( X, "X" );
        Print( Y, "Y" );
    }
    OutputFromRoot(g.Comm(),"Starting Trsm");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    Trsm( side, uplo, orientation, diag, alpha, A, Y );
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
    const auto SFrob = FrobeniusNorm( S );
    const auto XFrob = FrobeniusNorm( X );
    const auto EFrob = FrobeniusNorm( Y );
    OutputFromRoot
    (g.Comm(),
     "|| S ||_F = ",SFrob,"\n",Indent(),
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
            ("--uplo","lower or upper triangular: L/U",'L');
        const char transChar = Input
            ("--trans","orientation of triangular matrix: N/T/C",'N');
        const char diagChar = Input("--diag","(non-)unit diagonal: N/U",'N');
        const Int m = Input("--m","height of result",100);
        const Int n = Input("--n","width of result",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( gridHeight == 0 )
            gridHeight = Grid::DefaultHeight( mpi::Size(comm) );
        const GridOrder order = colMajor ? COLUMN_MAJOR : ROW_MAJOR;
        const Grid g( comm, gridHeight, order );
        const LeftOrRight side = CharToLeftOrRight( sideChar );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        OutputFromRoot
        (comm,"Will test Trsm ",sideChar,uploChar,transChar,diagChar);

        TestTrsm<float>
        ( side, uplo, orientation, diag,
          m, n,
          float(3),
          g, print );
        TestTrsm<Complex<float>>
        ( side, uplo, orientation, diag,
          m, n,
          Complex<float>(3),
          g, print );

        TestTrsm<double>
        ( side, uplo, orientation, diag,
          m, n,
          double(3),
          g, print );
        TestTrsm<Complex<double>>
        ( side, uplo, orientation, diag,
          m, n,
          Complex<double>(3),
          g, print );

#ifdef EL_HAVE_QD
        TestTrsm<DoubleDouble>
        ( side, uplo, orientation, diag,
          m, n,
          DoubleDouble(3),
          g, print );
        TestTrsm<QuadDouble>
        ( side, uplo, orientation, diag,
          m, n,
          QuadDouble(3),
          g, print );

        TestTrsm<Complex<DoubleDouble>>
        ( side, uplo, orientation, diag,
          m, n,
          Complex<DoubleDouble>(3),
          g, print );
        TestTrsm<Complex<QuadDouble>>
        ( side, uplo, orientation, diag,
          m, n,
          Complex<QuadDouble>(3),
          g, print );
#endif

#ifdef EL_HAVE_QUAD
        TestTrsm<Quad>
        ( side, uplo, orientation, diag,
          m, n,
          Quad(3),
          g, print );
        TestTrsm<Complex<Quad>>
        ( side, uplo, orientation, diag,
          m, n,
          Complex<Quad>(3),
          g, print );
#endif

#ifdef EL_HAVE_MPC
        TestTrsm<BigFloat>
        ( side, uplo, orientation, diag,
          m, n,
          BigFloat(3),
          g, print );
        TestTrsm<Complex<BigFloat>>
        ( side, uplo, orientation, diag,
          m, n,
          Complex<BigFloat>(3),
          g, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
