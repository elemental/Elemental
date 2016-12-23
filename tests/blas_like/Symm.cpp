/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename T>
void TestSymm
( bool conjugate,
  LeftOrRight side,
  UpperOrLower uplo,
  Int m,
  Int n,
  T alpha,
  T beta,
  const Grid& g,
  bool print )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<T>());
    PushIndent();

    DistMatrix<T> A(g), B(g), C(g);

    if( side == LEFT )
        Uniform( A, m, m );
    else
        Uniform( A, n, n );
    MakeSymmetric( uplo, A, conjugate );
    Uniform( B, m, n );
    Uniform( C, m, n );
    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
        Print( C, "C" );
    }

    const Int numRHS = 100;
    DistMatrix<T> X(g), Y(g);
    Uniform( X, m, numRHS );
    if( print )
        Print( X, "X" );
    if( side == LEFT )
    {
        // Y := alpha Symm(A) (B X) + beta C X
        DistMatrix<T> Z(g);
        Gemm( NORMAL, NORMAL, T(1), B, X, Z );
        Gemm( NORMAL, NORMAL, alpha, A, Z, Y );
        Gemm( NORMAL, NORMAL, beta, C, X, T(1), Y );
        if( print )
            Print( Y, "Y := alpha Symm(A) (B X) + beta C X" );
    }
    else
    {
        // Y := alpha B (Symm(A) X) + beta C X
        DistMatrix<T> Z(g);
        Gemm( NORMAL, NORMAL, T(1), A, X, Z );
        Gemm( NORMAL, NORMAL, alpha, B, Z, Y );
        Gemm( NORMAL, NORMAL, beta, C, X, T(1), Y );
        if( print )
            Print( Y, "Y := alpha B (Symm(A) X) + beta C X" );
    }
    const Base<T> YFrobNorm = FrobeniusNorm( Y );

    Timer timer;

    OutputFromRoot(g.Comm(),"Starting Symm...");
    mpi::Barrier( g.Comm() );
    timer.Start();
    Symm( side, uplo, alpha, A, B, beta, C, conjugate );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    const double mD = double(m);
    const double nD = double(n);
    const double realGFlops =
      ( side==LEFT ? 2.*mD*mD*nD : 2.*mD*nD*nD ) / (1.e9*runTime);
    const double gFlops = ( IsComplex<T>::value ? 4*realGFlops : realGFlops );
    OutputFromRoot
    (g.Comm(),"Finished after ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        if( side == LEFT )
            Print( C, BuildString("C := ",alpha," Symm(A) B + ",beta," C") );
        else
            Print( C, BuildString("C := ",alpha," B Symm(A) + ",beta," C") );
    }

    Gemm( NORMAL, NORMAL, T(-1), C, X, T(1), Y );
    const Base<T> EFrobNorm = FrobeniusNorm( Y );
    if( print )
        Print( Y, "E" );
    OutputFromRoot
    (g.Comm(),"|| E ||_F / || Y ||_F = ",
     EFrobNorm,"/",YFrobNorm,"=",EFrobNorm/YFrobNorm);

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
        const bool conjugate = Input("--conjugate","conjugate Symm?",false);
        const char sideChar = Input("--side","side to apply from: L/R",'L');
        const char uploChar = Input("--uplo","lower/upper storage: L/U",'L');
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
        SetBlocksize( nb );

        ComplainIfDebug();
        OutputFromRoot(comm,"Will test Symm ",sideChar,uploChar);

        TestSymm<float>
        ( conjugate, side, uplo,
          m, n,
          float(3), float(4),
          g, print );
        TestSymm<Complex<float>>
        ( conjugate, side, uplo,
          m, n,
          Complex<float>(3), Complex<float>(4),
          g, print );

        TestSymm<double>
        ( conjugate, side, uplo,
          m, n,
          double(3), double(4),
          g, print );
        TestSymm<Complex<double>>
        ( conjugate, side, uplo,
          m, n,
          Complex<double>(3), Complex<double>(4),
          g, print );

#ifdef EL_HAVE_QD
        TestSymm<DoubleDouble>
        ( conjugate, side, uplo,
          m, n,
          DoubleDouble(3), DoubleDouble(4),
          g, print );
        TestSymm<QuadDouble>
        ( conjugate, side, uplo,
          m, n,
          QuadDouble(3), QuadDouble(4),
          g, print );

        TestSymm<Complex<DoubleDouble>>
        ( conjugate, side, uplo,
          m, n,
          Complex<DoubleDouble>(3), Complex<DoubleDouble>(4),
          g, print );
        TestSymm<Complex<QuadDouble>>
        ( conjugate, side, uplo,
          m, n,
          Complex<QuadDouble>(3), Complex<QuadDouble>(4),
          g, print );
#endif

#ifdef EL_HAVE_QUAD
        TestSymm<Quad>
        ( conjugate, side, uplo,
          m, n,
          Quad(3), Quad(4),
          g, print );
        TestSymm<Complex<Quad>>
        ( conjugate, side, uplo,
          m, n,
          Complex<Quad>(3), Complex<Quad>(4),
          g, print );
#endif

#ifdef EL_HAVE_MPC
        TestSymm<BigFloat>
        ( conjugate, side, uplo,
          m, n,
          BigFloat(3), BigFloat(4),
          g, print );
        TestSymm<Complex<BigFloat>>
        ( conjugate, side, uplo,
          m, n,
          Complex<BigFloat>(3), Complex<BigFloat>(4),
          g, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
