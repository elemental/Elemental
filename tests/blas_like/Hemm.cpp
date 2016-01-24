/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template<typename T> 
void TestHemm
( LeftOrRight side,
  UpperOrLower uplo,
  Int m,
  Int n,
  T alpha,
  T beta,
  const Grid& g,
  bool print )
{
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<T>());
    DistMatrix<T> A(g), B(g), C(g);

    if( side == LEFT )
        HermitianUniformSpectrum( A, m, -10, 10 );
    else
        HermitianUniformSpectrum( A, n, -10, 10 );
    Uniform( B, m, n );
    Uniform( C, m, n );
    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
        Print( C, "C" );
    }

    if( g.Rank() == 0 )
        Output("  Starting Hemm...");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Hemm( side, uplo, alpha, A, B, beta, C );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double mD = double(m);
    const double nD = double(n);
    const double realGFlops = 
      ( side==LEFT ? 2.*mD*mD*nD : 2.*mD*nD*nD ) / (1.e9*runTime);
    const double gFlops = ( IsComplex<T>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
      Output("  Finished after ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        if( side == LEFT )
            Print( C, BuildString("C := ",alpha," Herm(A) B + ",beta," C") );
        else
            Print( C, BuildString("C := ",alpha," B Herm(A) + ",beta," C") );
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
        const char sideChar = Input("--side","side to apply from: L/R",'L');
        const char uploChar = Input("--uplo","lower/upper storage: L/U",'L');
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
        SetBlocksize( nb );

        ComplainIfDebug();
        if( commRank == 0 )
            Output("Will test Hemm ",sideChar,uploChar);

        TestHemm<float>
        ( side, uplo,
          m, n,
          float(3), float(4),
          g, print );
        TestHemm<Complex<float>>
        ( side, uplo,
          m, n,
          Complex<float>(3), Complex<float>(4),
          g, print );

        TestHemm<double>
        ( side, uplo,
          m, n,
          double(3), double(4),
          g, print );
        TestHemm<Complex<double>>
        ( side, uplo,
          m, n,
          Complex<double>(3), Complex<double>(4),
          g, print );

#ifdef EL_HAVE_QD
        TestHemm<DoubleDouble>
        ( side, uplo,
          m, n,
          DoubleDouble(3), DoubleDouble(4),
          g, print );
        TestHemm<QuadDouble>
        ( side, uplo,
          m, n,
          QuadDouble(3), QuadDouble(4),
          g, print );
#endif

#ifdef EL_HAVE_QUAD
        TestHemm<Quad>
        ( side, uplo,
          m, n,
          Quad(3), Quad(4),
          g, print );
        TestHemm<Complex<Quad>>
        ( side, uplo,
          m, n,
          Complex<Quad>(3), Complex<Quad>(4),
          g, print );
#endif

#ifdef EL_HAVE_MPC
        TestHemm<BigFloat>
        ( side, uplo,
          m, n,
          BigFloat(3), BigFloat(4),
          g, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
