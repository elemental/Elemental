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
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());
    DistMatrix<F> A(g), X(g);

    if( side == LEFT )
        HermitianUniformSpectrum( A, m, 1, 10 );
    else
        HermitianUniformSpectrum( A, n, 1, 10 );
    auto S( A );
    MakeTrapezoidal( uplo, S );

    Uniform( X, m, n );
    DistMatrix<F> Y(g);
    Gemm( NORMAL, NORMAL, F(1)/alpha, S, X, Y );

    if( print )
    {
        Print( A, "A" );
        Print( S, "S" );
        Print( X, "X" );
        Print( Y, "Y" );
    }
    if( g.Rank() == 0 )
        Output("  Starting Trsm");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Trsm( side, uplo, orientation, diag, alpha, A, Y );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 
      ( side==LEFT ? double(m)*double(m)*double(n)
                   : double(m)*double(n)*double(n) ) /(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
        Output("  Finished in ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
        Print( Y, "Y after solve" );
    Y -= X;
    const auto SFrob = FrobeniusNorm( S );
    const auto XFrob = FrobeniusNorm( X );
    const auto EFrob = FrobeniusNorm( Y );
    if( g.Rank() == 0 )
    {
        Output("  || S ||_F = ",SFrob);
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

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const LeftOrRight side = CharToLeftOrRight( sideChar );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        if( commRank == 0 )
            Output("Will test Trsm ",sideChar,uploChar,transChar,diagChar);

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
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
