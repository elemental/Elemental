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
( bool conjugated,
  bool print,
  const Matrix<Field>& A,
  const Matrix<Field>& dSub,
  const Permutation& p,
  const Matrix<Field>& AOrig,
        Int numRHS=100 )
{
    typedef Base<Field> Real;
    const Int m = AOrig.Height();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = HermitianOneNorm( LOWER, AOrig );

    Matrix<Field> X, Y;
    Uniform( X, m, numRHS );
    Y = X;
    const Real oneNormX = OneNorm( X );

    // Test correctness by comparing the application of AOrig against a
    // random set of 100 vectors to the application of tril(A) tril(A)^H
    if( print )
        Print( X, "X" );
    ldl::MultiplyAfter( A, dSub, p, Y, conjugated );
    if( print )
        Print( Y, "P' L B L' P X" );
    Symm( LEFT, LOWER, Field(-1), AOrig, X, Field(1), Y, conjugated );
    if( print )
        Print( Y, "P' L B L' P X - A X" );
    const Real infError = InfinityNorm( Y );
    const Real relError = infError / (m*eps*Max(oneNormA,oneNormX));

    Output
    ("||A X - L D L^[T/H] X||_oo / (eps m Max(||A||_1,||X||_1)) = ",relError);

    // TODO: A more refined failure condition
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably high");
}

template<typename Field>
void TestCorrectness
( bool conjugated,
  bool print,
  const DistMatrix<Field>& A,
  const DistMatrix<Field,MD,STAR>& dSub,
  const DistPermutation& p,
  const DistMatrix<Field>& AOrig,
        Int numRHS=100 )
{
    typedef Base<Field> Real;
    const Grid& grid = A.Grid();
    const Int m = AOrig.Height();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = HermitianOneNorm( LOWER, AOrig );

    DistMatrix<Field> X(grid), Y(grid);
    Uniform( X, m, numRHS );
    Y = X;
    const Real oneNormX = OneNorm( X );

    // Test correctness by comparing the application of AOrig against a
    // random set of 100 vectors to the application of tril(A) tril(A)^H
    if( print )
        Print( X, "X" );
    ldl::MultiplyAfter( A, dSub, p, Y, conjugated );
    if( print )
        Print( Y, "P' L B L' P X" );
    Symm( LEFT, LOWER, Field(-1), AOrig, X, Field(1), Y, conjugated );
    if( print )
        Print( Y, "P' L B L' P X - A X" );
    const Real infError = InfinityNorm( Y );
    const Real relError = infError / (m*eps*Max(oneNormA,oneNormX));

    OutputFromRoot
    (grid.Comm(),
     "||A X - L D L^[T/H] X||_oo / (eps m Max(||A||_1,||X||_1)) = ",relError);

    // TODO(poulson): A more refined failure condition
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably high");
}

template<typename Field>
void TestLDL
( Int m,
  bool conjugated,
  Int nbLocal,
  bool correctness,
  bool print )
{
    Matrix<Field> A, AOrig;
    Output("Testing with ",TypeName<Field>());
    PushIndent();

    SetLocalTrrkBlocksize<Field>( nbLocal );

    if( conjugated )
        HermitianUniformSpectrum( A, m, -100, 100 );
    else
        Uniform( A, m, m );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    Output("Starting LDL^[T/H] factorization...");
    Timer timer;
    timer.Start();
    Matrix<Field> dSub;
    Permutation p;
    LDL( A, dSub, p, conjugated );
    const double runTime = timer.Stop();
    const double realGFlops = 1./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = IsComplex<Field>::value ? 4*realGFlops : realGFlops;
    Output(runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        Matrix<Int> P;
        p.ExplicitMatrix( P );
        Print( P, "P" );
    }
    if( correctness )
        TestCorrectness( conjugated, print, A, dSub, p, AOrig );
    PopIndent();
}

template<typename Field>
void TestLDL
( const Grid& grid,
  Int m,
  bool conjugated,
  Int nbLocal,
  bool correctness,
  bool print )
{
    DistMatrix<Field> A(grid), AOrig(grid);
    OutputFromRoot(grid.Comm(),"Testing with ",TypeName<Field>());
    PushIndent();

    SetLocalTrrkBlocksize<Field>( nbLocal );

    if( conjugated )
        HermitianUniformSpectrum( A, m, -100, 100 );
    else
        Uniform( A, m, m );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    OutputFromRoot(grid.Comm(),"Starting LDL^[T/H] factorization...");
    mpi::Barrier( grid.Comm() );
    Timer timer;
    timer.Start();
    DistMatrix<Field,MD,STAR> dSub(grid);
    DistPermutation p(grid);
    LDL( A, dSub, p, conjugated );
    mpi::Barrier( grid.Comm() );
    const double runTime = timer.Stop();
    const double realGFlops = 1./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = IsComplex<Field>::value ? 4*realGFlops : realGFlops;
    OutputFromRoot(grid.Comm(),runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        DistMatrix<Int> P(grid);
        p.ExplicitMatrix( P );
        Print( P, "P" );
    }
    if( correctness )
        TestCorrectness( conjugated, print, A, dSub, p, AOrig );
    PopIndent();
}

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        int gridHeight = Input("--gridHeight","process grid height",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int m = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool conjugated = Input("--conjugate","conjugate LDL?",false);
        const bool sequential = Input("--sequential","test sequential?",true);
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
        const GridOrder order = colMajor ? COLUMN_MAJOR : ROW_MAJOR;
        const Grid grid( comm, gridHeight, order );
        SetBlocksize( nb );
        ComplainIfDebug();

        if( sequential && mpi::Rank() == 0 )
        {
            TestLDL<float>
            ( m, conjugated, nbLocal, correctness, print );
            TestLDL<Complex<float>>
            ( m, conjugated, nbLocal, correctness, print );

            TestLDL<double>
            ( m, conjugated, nbLocal, correctness, print );
            TestLDL<Complex<double>>
            ( m, conjugated, nbLocal, correctness, print );

#ifdef EL_HAVE_QD
            TestLDL<DoubleDouble>
            ( m, conjugated, nbLocal, correctness, print );
            TestLDL<QuadDouble>
            ( m, conjugated, nbLocal, correctness, print );

            TestLDL<Complex<DoubleDouble>>
            ( m, conjugated, nbLocal, correctness, print );
            TestLDL<Complex<QuadDouble>>
            ( m, conjugated, nbLocal, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
            TestLDL<Quad>
            ( m, conjugated, nbLocal, correctness, print );
            TestLDL<Complex<Quad>>
            ( m, conjugated, nbLocal, correctness, print );
#endif

#ifdef EL_HAVE_MPC
            TestLDL<BigFloat>
            ( m, conjugated, nbLocal, correctness, print );
            TestLDL<Complex<BigFloat>>
            ( m, conjugated, nbLocal, correctness, print );
#endif
        }

        TestLDL<float>
        ( grid, m, conjugated, nbLocal, correctness, print );
        TestLDL<Complex<float>>
        ( grid, m, conjugated, nbLocal, correctness, print );

        TestLDL<double>
        ( grid, m, conjugated, nbLocal, correctness, print );
        TestLDL<Complex<double>>
        ( grid, m, conjugated, nbLocal, correctness, print );

#ifdef EL_HAVE_QD
        TestLDL<DoubleDouble>
        ( grid, m, conjugated, nbLocal, correctness, print );
        TestLDL<QuadDouble>
        ( grid, m, conjugated, nbLocal, correctness, print );

        TestLDL<Complex<DoubleDouble>>
        ( grid, m, conjugated, nbLocal, correctness, print );
        TestLDL<Complex<QuadDouble>>
        ( grid, m, conjugated, nbLocal, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
        TestLDL<Quad>
        ( grid, m, conjugated, nbLocal, correctness, print );
        TestLDL<Complex<Quad>>
        ( grid, m, conjugated, nbLocal, correctness, print );
#endif

#ifdef EL_HAVE_MPC
        TestLDL<BigFloat>
        ( grid, m, conjugated, nbLocal, correctness, print );
        TestLDL<Complex<BigFloat>>
        ( grid, m, conjugated, nbLocal, correctness, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
