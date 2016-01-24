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
void TestSymm
( LeftOrRight side,
  UpperOrLower uplo,
  Int m,
  Int n,
  T alpha,
  T beta,
  bool print,
  const Grid& g )
{
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<T>());
    DistMatrix<T> A(g), B(g), C(g);

    if( side == LEFT )
        Uniform( A, m, m );
    else
        Uniform( A, n, n );
    Uniform( B, m, n );
    Uniform( C, m, n );
    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
        Print( C, "C" );
    }

    // Test Symm
    if( g.Rank() == 0 )
        Output("  Starting Symm");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Symm( side, uplo, alpha, A, B, beta, C );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double mD = double(m);
    const double nD = double(n);
    const double realGFlops = 
      ( side==LEFT ? 2.*mD*mD*nD : 2.*mD*nD*nD ) / (1.e9*runTime);
    const double gFlops = ( IsComplex<T>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
        Output("  Finished in ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        if( side == LEFT )
            Print( C, BuildString("C := ",alpha," Symm(A) B + ",beta," C") );
        else
            Print( C, BuildString("C := ",alpha," B Symm(A) + ",beta," C") );
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
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
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
            Output("Will test Symm ",sideChar,uploChar);

        TestSymm<float>
        ( side, uplo, m, n,
          float(3), float(4),
          print, g ); 
        TestSymm<Complex<float>>
        ( side, uplo, m, n,
          Complex<float>(3), Complex<float>(4),
          print, g ); 

        TestSymm<double>
        ( side, uplo, m, n,
          double(3), double(4),
          print, g ); 
        TestSymm<Complex<double>>
        ( side, uplo, m, n,
          Complex<double>(3), Complex<double>(4),
          print, g ); 

#ifdef EL_HAVE_QD
        TestSymm<DoubleDouble>
        ( side, uplo, m, n,
          DoubleDouble(3), DoubleDouble(4),
          print, g ); 
        TestSymm<QuadDouble>
        ( side, uplo, m, n,
          QuadDouble(3), QuadDouble(4),
          print, g ); 
#endif

#ifdef EL_HAVE_QUAD
        TestSymm<Quad>
        ( side, uplo, m, n,
          Quad(3), Quad(4),
          print, g ); 
        TestSymm<Complex<Quad>>
        ( side, uplo, m, n,
          Complex<Quad>(3), Complex<Quad>(4),
          print, g ); 
#endif

#ifdef EL_HAVE_MPC
        TestSymm<BigFloat>
        ( side, uplo, m, n,
          BigFloat(3), BigFloat(4),
          print, g ); 
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
