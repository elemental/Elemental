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
void TestCorrectness
( bool pivot,
  UpperOrLower uplo,
  const Matrix<F>& A,
  const Permutation& p,
  const Matrix<F>& AOrig,
        Int numRHS=100 )
{
    typedef Base<F> Real;
    const Int n = AOrig.Height();
    const Real eps = limits::Epsilon<Real>();

    // Test correctness by multiplying a random set of vectors by A, then
    // using the Cholesky factorization to solve.
    Matrix<F> X, Y;
    Uniform( X, n, numRHS );
    Zeros( Y, n, numRHS );
    Hemm( LEFT, uplo, F(1), AOrig, X, F(0), Y );
    const Real oneNormY = OneNorm( Y );

    if( pivot )
        cholesky::SolveAfter( uplo, NORMAL, A, p, Y );
    else
        cholesky::SolveAfter( uplo, NORMAL, A, Y );
    X -= Y;
    const Real infNormE = InfinityNorm( X );
    const Real relErr = infNormE / (eps*n*oneNormY);

    Output("||X - A \\ Y ||_oo / (eps n || Y ||_1) = ",relErr);
    // TODO(poulson): Use more refined failure criteria
    if( relErr > Real(100) )
        LogicError("Relative error was unacceptably large");
}

template<typename F>
void TestCorrectness
( bool pivot,
  UpperOrLower uplo,
  const DistMatrix<F>& A,
  const DistPermutation& p,
  const DistMatrix<F>& AOrig,
        Int numRHS=100 )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = AOrig.Height();
    const Real eps = limits::Epsilon<Real>();

    // Test correctness by multiplying a random set of vectors by A, then
    // using the Cholesky factorization to solve.
    DistMatrix<F> X(g), Y(g);
    Uniform( X, n, numRHS );
    Zeros( Y, n, numRHS );
    Hemm( LEFT, uplo, F(1), AOrig, X, F(0), Y );
    const Real oneNormY = OneNorm( Y );

    if( pivot )
        cholesky::SolveAfter( uplo, NORMAL, A, p, Y );
    else
        cholesky::SolveAfter( uplo, NORMAL, A, Y );
    X -= Y;
    const Real infNormE = InfinityNorm( X );
    const Real relErr = infNormE / (eps*n*oneNormY);

    OutputFromRoot
    (g.Comm(), "||X - A \\ Y ||_oo / (eps n || Y ||_1) = ",relErr);
    // TODO(poulson): Use more refined failure criteria
    if( relErr > Real(100) )
        LogicError("Relative error was unacceptably large");
}

template<typename F>
void TestSequentialCholesky
( UpperOrLower uplo,
  bool pivot,
  Int m,
  bool print,
  bool printDiag,
  bool correctness )
{
    Output("Testing sequential Cholesky with ",TypeName<F>());
    PushIndent();
    Matrix<F> A, AOrig;
    Permutation p;

    HermitianUniformSpectrum( A, m, 1e-9, 10 );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    Output("Cholesky...");
    Timer timer;
    timer.Start();
    if( pivot )
        Cholesky( uplo, A, p );
    else
        Cholesky( uplo, A );
    const double runTime = timer.Stop();
    const double realGFlops = (1./3.)*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    Output(runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        if( pivot )
        {
            Matrix<Int> P;
            p.ExplicitMatrix( P );
            Print( P, "P" );
        }
    }
    if( printDiag )
        Print( GetRealPartOfDiagonal(A), "diag(A)" );
    if( correctness )
        TestCorrectness( pivot, uplo, A, p, AOrig );
    PopIndent();
}

template<typename F>
void TestCholesky
( const Grid& g,
  UpperOrLower uplo,
  bool pivot,
  Int m,
  Int nbLocal,
  bool print,
  bool printDiag,
  bool correctness,
  bool scalapack )
{
    OutputFromRoot(g.Comm(),"Testing distributed Cholesky with ",TypeName<F>());
    PushIndent();
    DistMatrix<F> A(g), AOrig(g);
    DistPermutation p(g);

    SetLocalTrrkBlocksize<F>( nbLocal );

    HermitianUniformSpectrum( A, m, 1e-9, 10 );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    if( scalapack && !pivot )
        OutputFromRoot
        (g.Comm(),"ScaLAPACK Cholesky (including round-trip conversion)...");
    else
        OutputFromRoot(g.Comm(),"Elemental Cholesky...");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    if( pivot )
        Cholesky( uplo, A, p );
    else
        Cholesky( uplo, A, scalapack );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    const double realGFlops = 1./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    OutputFromRoot(g.Comm(),runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        if( pivot )
        {
            DistMatrix<Int,VC,STAR> P(g);
            p.ExplicitMatrix( P );
            Print( P, "P" );
        }
    }
    if( printDiag )
        Print( GetRealPartOfDiagonal(A), "diag(A)" );
    if( correctness )
        TestCorrectness( pivot, uplo, A, p, AOrig );
    PopIndent();
}

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        Int gridHeight = Input("--gridHeight","process grid height",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int m = Input("--m","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool pivot = Input("--pivot","use pivoting?",false);
        const bool correctness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool printDiag = Input("--printDiag","print diag of fact?",false);
        const bool sequential = Input("--sequential","test sequential?",true);
#ifdef EL_HAVE_SCALAPACK
        const bool scalapack = Input("--scalapack","test ScaLAPACK?",false);
#else
        const bool scalapack = false;
#endif
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        ProcessInput();
        PrintInputReport();

#ifdef EL_HAVE_MPC
        mpfr::SetPrecision( prec );
#endif

        if( gridHeight == 0 )
            gridHeight = Grid::DefaultHeight( mpi::Size(comm) );
        const GridOrder order = colMajor ? COLUMN_MAJOR : ROW_MAJOR;
        const Grid g( comm, gridHeight, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );

        ComplainIfDebug();

        if( sequential && mpi::Rank(comm) == 0 )
        {
            TestSequentialCholesky<float>
            ( uplo, pivot, m, print, printDiag, correctness );
            TestSequentialCholesky<Complex<float>>
            ( uplo, pivot, m, print, printDiag, correctness );
            TestSequentialCholesky<double>
            ( uplo, pivot, m, print, printDiag, correctness );
            TestSequentialCholesky<Complex<double>>
            ( uplo, pivot, m, print, printDiag, correctness );

#ifdef EL_HAVE_QD
            TestSequentialCholesky<DoubleDouble>
            ( uplo, pivot, m, print, printDiag, correctness );
            TestSequentialCholesky<QuadDouble>
            ( uplo, pivot, m, print, printDiag, correctness );

            TestSequentialCholesky<Complex<DoubleDouble>>
            ( uplo, pivot, m, print, printDiag, correctness );
            TestSequentialCholesky<Complex<QuadDouble>>
            ( uplo, pivot, m, print, printDiag, correctness );
#endif

#ifdef EL_HAVE_QUAD
            TestSequentialCholesky<Quad>
            ( uplo, pivot, m, print, printDiag, correctness );
            TestSequentialCholesky<Complex<Quad>>
            ( uplo, pivot, m, print, printDiag, correctness );
#endif

#ifdef EL_HAVE_MPC
            TestSequentialCholesky<BigFloat>
            ( uplo, pivot, m, print, printDiag, correctness );
            TestSequentialCholesky<Complex<BigFloat>>
            ( uplo, pivot, m, print, printDiag, correctness );
#endif
        }

        TestCholesky<float>
        ( g, uplo, pivot, m, nbLocal,
          print, printDiag, correctness, scalapack );
        TestCholesky<Complex<float>>
        ( g, uplo, pivot, m, nbLocal,
          print, printDiag, correctness, scalapack );
        TestCholesky<double>
        ( g, uplo, pivot, m, nbLocal,
          print, printDiag, correctness, scalapack );
        TestCholesky<Complex<double>>
        ( g, uplo, pivot, m, nbLocal,
          print, printDiag, correctness, scalapack );

#ifdef EL_HAVE_QD
        TestCholesky<DoubleDouble>
        ( g, uplo, pivot, m, nbLocal,
          print, printDiag, correctness, scalapack );
        TestCholesky<QuadDouble>
        ( g, uplo, pivot, m, nbLocal,
          print, printDiag, correctness, scalapack );

        TestCholesky<Complex<DoubleDouble>>
        ( g, uplo, pivot, m, nbLocal,
          print, printDiag, correctness, scalapack );
        TestCholesky<Complex<QuadDouble>>
        ( g, uplo, pivot, m, nbLocal,
          print, printDiag, correctness, scalapack );
#endif

#ifdef EL_HAVE_QUAD
        TestCholesky<Quad>
        ( g, uplo, pivot, m, nbLocal,
          print, printDiag, correctness, scalapack );
        TestCholesky<Complex<Quad>>
        ( g, uplo, pivot, m, nbLocal,
          print, printDiag, correctness, scalapack );
#endif

#ifdef EL_HAVE_MPC
        TestCholesky<BigFloat>
        ( g, uplo, pivot, m, nbLocal,
          print, printDiag, correctness, scalapack );
        TestCholesky<Complex<BigFloat>>
        ( g, uplo, pivot, m, nbLocal,
          print, printDiag, correctness, scalapack );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
