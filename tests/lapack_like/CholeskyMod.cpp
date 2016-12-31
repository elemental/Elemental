/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename Field>
void TestCorrectness
(       UpperOrLower uplo,
  const Matrix<Field>& T,
        Base<Field> alpha,
  const Matrix<Field>& V,
  const Matrix<Field>& A,
        Int numRHS=100 )
{
    typedef Base<Field> Real;
    const Int n = V.Height();
    const Real eps = limits::Epsilon<Real>();

    Matrix<Field> B( A );
    Herk( uplo, NORMAL, alpha, V, Real(1), B );

    // Test correctness by multiplying a random set of vectors by
    // A + alpha V V^H, then using the Cholesky factorization to solve.
    Matrix<Field> X, Y;
    Uniform( X, n, numRHS );
    Zeros( Y, n, numRHS );
    Hemm( LEFT, uplo, Field(1), B, X, Field(0), Y );
    const Real oneNormY = OneNorm( Y );

    cholesky::SolveAfter( uplo, NORMAL, T, Y );
    X -= Y;
    const Real infNormE = InfinityNorm( X );
    const Real relError = infNormE / (eps*n*oneNormY);

    Output("|| X - B \\ Y ||_oo / (n eps || Y ||_1) = ",relError);

    // TODO(poulson): Use a more refined failure condition
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large");
}

template<typename Field>
void TestCorrectness
(       UpperOrLower uplo,
  const DistMatrix<Field>& T,
        Base<Field> alpha,
  const DistMatrix<Field>& V,
  const DistMatrix<Field>& A,
        Int numRHS=100 )
{
    typedef Base<Field> Real;
    const Int n = V.Height();
    const Real eps = limits::Epsilon<Real>();
    const Grid& grid = T.Grid();

    DistMatrix<Field> B( A );
    Herk( uplo, NORMAL, alpha, V, Real(1), B );

    // Test correctness by multiplying a random set of vectors by
    // A + alpha V V^H, then using the Cholesky factorization to solve.
    DistMatrix<Field> X(grid), Y(grid);
    Uniform( X, n, numRHS );
    Zeros( Y, n, numRHS );
    Hemm( LEFT, uplo, Field(1), B, X, Field(0), Y );
    const Real oneNormY = OneNorm( Y );

    cholesky::SolveAfter( uplo, NORMAL, T, Y );
    X -= Y;
    const Real infNormE = InfinityNorm( X );
    const Real relError = infNormE / (eps*n*oneNormY);

    OutputFromRoot
    (grid.Comm(),"|| X - B \\ Y ||_oo / (n eps || Y ||_1) = ",relError);

    // TODO(poulson): Use a more refined failure condition
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large");
}

template<typename Field>
void TestCholeskyMod
( UpperOrLower uplo,
  Int m,
  Int n,
  Base<Field> alpha,
  bool correctness,
  bool print )
{
    Output("Testing with ",TypeName<Field>());
    PushIndent();

    Matrix<Field> T, A;
    HermitianUniformSpectrum( T, m, 1e-9, 10 );
    if( correctness )
        A = T;
    if( print )
        Print( T, "A" );

    Output("Starting Cholesky...");
    Timer timer;
    timer.Start();
    Cholesky( uplo, T );
    double runTime = timer.Stop();
    double realGFlops = 1./3.*Pow(double(m),3.)/(1.e9*runTime);
    double gFlops = ( IsComplex<Field>::value ? 4*realGFlops : realGFlops );
    Output(runTime," seconds (",gFlops," GFlop/s)");
    MakeTrapezoidal( uplo, T );
    if( print )
        Print( T, "Cholesky factor" );

    Matrix<Field> V, VMod;
    Uniform( V, m, n );
    V *= Field(1)/Sqrt(Field(m)*Field(n));
    VMod = V;
    if( print )
        Print( V, "V" );

    Output("Starting Cholesky mod...");
    timer.Start();
    CholeskyMod( uplo, T, alpha, VMod );
    runTime = timer.Stop();
    Output(runTime," seconds");
    if( print )
        Print( T, "Modified Cholesky factor" );

    if( correctness )
        TestCorrectness( uplo, T, alpha, V, A );
    PopIndent();
}

template<typename Field>
void TestCholeskyMod
( const Grid& grid,
  UpperOrLower uplo,
  Int m,
  Int n,
  Base<Field> alpha,
  bool correctness,
  bool print )
{
    OutputFromRoot(grid.Comm(),"Testing with ",TypeName<Field>());
    PushIndent();

    DistMatrix<Field> T(grid), A(grid);
    HermitianUniformSpectrum( T, m, 1e-9, 10 );
    if( correctness )
        A = T;
    if( print )
        Print( T, "A" );

    OutputFromRoot(grid.Comm(),"Starting Cholesky...");
    Timer timer;
    timer.Start();
    Cholesky( uplo, T );
    double runTime = timer.Stop();
    double realGFlops = 1./3.*Pow(double(m),3.)/(1.e9*runTime);
    double gFlops = ( IsComplex<Field>::value ? 4*realGFlops : realGFlops );
    OutputFromRoot(grid.Comm(),runTime," seconds (",gFlops," GFlop/s)");
    MakeTrapezoidal( uplo, T );
    if( print )
        Print( T, "Cholesky factor" );

    DistMatrix<Field> V(grid), VMod(grid);
    Uniform( V, m, n );
    V *= Field(1)/Sqrt(Field(m)*Field(n));
    VMod = V;
    if( print )
        Print( V, "V" );

    OutputFromRoot(grid.Comm(),"Starting Cholesky mod...");
    timer.Start();
    CholeskyMod( uplo, T, alpha, VMod );
    runTime = timer.Stop();
    OutputFromRoot(grid.Comm(),runTime," seconds");
    if( print )
        Print( T, "Modified Cholesky factor" );

    if( correctness )
        TestCorrectness( uplo, T, alpha, V, A );
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
        const Int n = Input("--n","rank of update",5);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const double alpha = Input("--alpha","update scaling",3.);
        const bool sequential =
          Input("--sequential","test sequential?",true);
        const bool correctness =
          Input("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
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
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        ComplainIfDebug();

        if( sequential && mpi::Rank() == 0 )
        {
            TestCholeskyMod<float>
            ( uplo, m, n, alpha, correctness, print );
            TestCholeskyMod<Complex<float>>
            ( uplo, m, n, alpha, correctness, print );

            TestCholeskyMod<double>
            ( uplo, m, n, alpha, correctness, print );
            TestCholeskyMod<Complex<double>>
            ( uplo, m, n, alpha, correctness, print );

#ifdef EL_HAVE_QD
            TestCholeskyMod<DoubleDouble>
            ( uplo, m, n, alpha, correctness, print );
            TestCholeskyMod<QuadDouble>
            ( uplo, m, n, alpha, correctness, print );

            TestCholeskyMod<Complex<DoubleDouble>>
            ( uplo, m, n, alpha, correctness, print );
            TestCholeskyMod<Complex<QuadDouble>>
            ( uplo, m, n, alpha, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
            TestCholeskyMod<Quad>
            ( uplo, m, n, alpha, correctness, print );
            TestCholeskyMod<Complex<Quad>>
            ( uplo, m, n, alpha, correctness, print );
#endif

#ifdef EL_HAVE_MPC
            TestCholeskyMod<BigFloat>
            ( uplo, m, n, alpha, correctness, print );
            TestCholeskyMod<Complex<BigFloat>>
            ( uplo, m, n, alpha, correctness, print );
#endif
        }

        TestCholeskyMod<float>
        ( g, uplo, m, n, alpha, correctness, print );
        TestCholeskyMod<Complex<float>>
        ( g, uplo, m, n, alpha, correctness, print );

        TestCholeskyMod<double>
        ( g, uplo, m, n, alpha, correctness, print );
        TestCholeskyMod<Complex<double>>
        ( g, uplo, m, n, alpha, correctness, print );

#ifdef EL_HAVE_QD
        TestCholeskyMod<DoubleDouble>
        ( g, uplo, m, n, alpha, correctness, print );
        TestCholeskyMod<QuadDouble>
        ( g, uplo, m, n, alpha, correctness, print );

        TestCholeskyMod<Complex<DoubleDouble>>
        ( g, uplo, m, n, alpha, correctness, print );
        TestCholeskyMod<Complex<QuadDouble>>
        ( g, uplo, m, n, alpha, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
        TestCholeskyMod<Quad>
        ( g, uplo, m, n, alpha, correctness, print );
        TestCholeskyMod<Complex<Quad>>
        ( g, uplo, m, n, alpha, correctness, print );
#endif

#ifdef EL_HAVE_MPC
        TestCholeskyMod<BigFloat>
        ( g, uplo, m, n, alpha, correctness, print );
        TestCholeskyMod<Complex<BigFloat>>
        ( g, uplo, m, n, alpha, correctness, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
