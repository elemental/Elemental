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
( bool print,
  const Matrix<Field>& A,
  const Permutation& P,
  const Matrix<Field>& AOrig,
        Int numRHS=100 )
{
    typedef Base<Field> Real;
    const Int n = AOrig.Width();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( AOrig );

    Output("Testing error...");

    // Generate random right-hand sides
    Matrix<Field> X;
    Uniform( X, n, numRHS );
    auto Y( X );
    const Real oneNormY = OneNorm( Y );
    P.PermuteRows( Y );
    lu::SolveAfter( NORMAL, A, Y );

    // Now investigate the residual, ||AOrig Y - X||_oo
    Gemm( NORMAL, NORMAL, Field(-1), AOrig, Y, Field(1), X );
    const Real infError = InfinityNorm( X );
    const Real relError = infError / (eps*n*Max(oneNormA,oneNormY));

    // TODO(poulson): Use a rigorous failure condition
    Output("||A X - Y||_oo / (eps n Max(||A||_1,||Y||_1)) = ",relError);
    if( relError > Real(1000) )
        LogicError("Unacceptably large relative error");
}

template<typename Field>
void TestCorrectness
( bool print,
  const DistMatrix<Field>& A,
  const DistPermutation& P,
  const DistMatrix<Field>& AOrig,
        Int numRHS=100 )
{
    typedef Base<Field> Real;
    const Grid& grid = A.Grid();
    const Int n = AOrig.Width();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( AOrig );

    OutputFromRoot(grid.Comm(),"Testing error...");

    // Generate random right-hand sides
    DistMatrix<Field> X(grid);
    Uniform( X, n, numRHS );
    auto Y( X );
    const Real oneNormY = OneNorm( Y );
    P.PermuteRows( Y );
    lu::SolveAfter( NORMAL, A, Y );

    // Now investigate the residual, ||AOrig Y - X||_oo
    Gemm( NORMAL, NORMAL, Field(-1), AOrig, Y, Field(1), X );
    const Real infError = InfinityNorm( X );
    const Real relError = infError / (eps*n*Max(oneNormA,oneNormY));

    // TODO(poulson): Use a rigorous failure condition
    OutputFromRoot
    (grid.Comm(),"||A X - Y||_oo / (eps n Max(||A||_1,||Y||_1)) = ",relError);
    if( relError > Real(1000) )
        LogicError("Unacceptably large relative error");
}

template<typename Field>
void TestLUMod
( Int m,
  Int updateRank,
  bool conjugate,
  Base<Field> tau,
  bool correctness,
  bool print )
{
    Output("Testing with ",TypeName<Field>());
    PushIndent();
    Matrix<Field> A, AOrig;
    Permutation P;

    Uniform( A, m, m );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    {
        Output("Starting LU factorization...");
        Timer timer;
        timer.Start();
        //P.ReserveSwaps( m+2*m-1 );
        LU( A, P );
        const double runTime = timer.Stop();
        const double realGFlops = 2./3.*Pow(double(m),3.)/(1.e9*runTime);
        const double gFlops =
          ( IsComplex<Field>::value ? 4*realGFlops : realGFlops );
        Output(runTime," seconds (",gFlops," GFlop/s)");
    }

    // TODO(poulson): Print permutation

    // Generate random vectors u and v
    Matrix<Field> U, V;
    Uniform( U, m, updateRank );
    Uniform( V, m, updateRank );
    if( correctness )
    {
        if( conjugate )
            Gemm( NORMAL, ADJOINT, Field(1), U, V, Field(1), AOrig );
        else
            Gemm( NORMAL, TRANSPOSE, Field(1), U, V, Field(1), AOrig );
    }

    {
        Output("Starting low-rank LU modification...");
        Timer timer;
        timer.Start();
        LUMod( A, P, U, V, conjugate, tau );
        const double runTime = timer.Stop();
        Output(runTime," seconds");
    }

    // TODO(poulson): Print permutation

    if( correctness )
        TestCorrectness( print, A, P, AOrig );

    PopIndent();
}

template<typename Field>
void TestLUMod
( const Grid& grid,
  Int m,
  Int updateRank,
  bool conjugate,
  Base<Field> tau,
  bool correctness,
  bool print )
{
    OutputFromRoot(grid.Comm(),"Testing with ",TypeName<Field>());
    PushIndent();
    DistMatrix<Field> A(grid), AOrig(grid);
    DistPermutation P(grid);

    Uniform( A, m, m );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    {
        OutputFromRoot(grid.Comm(),"Starting LU factorization...");
        mpi::Barrier( grid.Comm() );
        Timer timer;
        timer.Start();
        P.ReserveSwaps( m+2*m-1 );
        LU( A, P );
        mpi::Barrier( grid.Comm() );
        const double runTime = timer.Stop();
        const double realGFlops = 2./3.*Pow(double(m),3.)/(1.e9*runTime);
        const double gFlops =
          ( IsComplex<Field>::value ? 4*realGFlops : realGFlops );
        OutputFromRoot(grid.Comm(),runTime," seconds (",gFlops," GFlop/s)");
    }

    // TODO(poulson): Print permutation

    // Generate random vectors u and v
    DistMatrix<Field> U(grid), V(grid);
    Uniform( U, m, updateRank );
    Uniform( V, m, updateRank );
    if( correctness )
    {
        if( conjugate )
            Gemm( NORMAL, ADJOINT, Field(1), U, V, Field(1), AOrig );
        else
            Gemm( NORMAL, TRANSPOSE, Field(1), U, V, Field(1), AOrig );
    }

    {
        OutputFromRoot(grid.Comm(),"Starting low-rank LU modification...");
        mpi::Barrier( grid.Comm() );
        Timer timer;
        timer.Start();
        LUMod( A, P, U, V, conjugate, tau );
        mpi::Barrier( grid.Comm() );
        const double runTime = timer.Stop();
        OutputFromRoot(grid.Comm(),runTime," seconds");
    }

    // TODO(poulson): Print permutation

    if( correctness )
        TestCorrectness( print, A, P, AOrig );
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
        const Int m = Input("--height","height of matrix",100);
        const Int updateRank = Input("--updateRank","rank of LU update",10);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const double tau = Input("--tau","pivot threshold",0.1);
        const bool conjugate = Input("--conjugate","conjugate v?",true);
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
            TestLUMod<float>
            ( m, updateRank, conjugate, tau, correctness, print );
            TestLUMod<Complex<float>>
            ( m, updateRank, conjugate, tau, correctness, print );

            TestLUMod<double>
            ( m, updateRank, conjugate, tau, correctness, print );
            TestLUMod<Complex<double>>
            ( m, updateRank, conjugate, tau, correctness, print );

#ifdef EL_HAVE_QD
            TestLUMod<DoubleDouble>
            ( m, updateRank, conjugate, tau, correctness, print );
            TestLUMod<QuadDouble>
            ( m, updateRank, conjugate, tau, correctness, print );

            TestLUMod<Complex<DoubleDouble>>
            ( m, updateRank, conjugate, tau, correctness, print );
            TestLUMod<Complex<QuadDouble>>
            ( m, updateRank, conjugate, tau, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
            TestLUMod<Quad>
            ( m, updateRank, conjugate, tau, correctness, print );
            TestLUMod<Complex<Quad>>
            ( m, updateRank, conjugate, tau, correctness, print );
#endif

#ifdef EL_HAVE_MPC
            TestLUMod<BigFloat>
            ( m, updateRank, conjugate, tau, correctness, print );
            TestLUMod<Complex<BigFloat>>
            ( m, updateRank, conjugate, tau, correctness, print );
#endif
        }

        TestLUMod<float>
        ( grid, m, updateRank, conjugate, tau, correctness, print );
        TestLUMod<Complex<float>>
        ( grid, m, updateRank, conjugate, tau, correctness, print );

        TestLUMod<double>
        ( grid, m, updateRank, conjugate, tau, correctness, print );
        TestLUMod<Complex<double>>
        ( grid, m, updateRank, conjugate, tau, correctness, print );

#ifdef EL_HAVE_QD
        TestLUMod<DoubleDouble>
        ( grid, m, updateRank, conjugate, tau, correctness, print );
        TestLUMod<QuadDouble>
        ( grid, m, updateRank, conjugate, tau, correctness, print );

        TestLUMod<Complex<DoubleDouble>>
        ( grid, m, updateRank, conjugate, tau, correctness, print );
        TestLUMod<Complex<QuadDouble>>
        ( grid, m, updateRank, conjugate, tau, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
        TestLUMod<Quad>
        ( grid, m, updateRank, conjugate, tau, correctness, print );
        TestLUMod<Complex<Quad>>
        ( grid, m, updateRank, conjugate, tau, correctness, print );
#endif

#ifdef EL_HAVE_MPC
        TestLUMod<BigFloat>
        ( grid, m, updateRank, conjugate, tau, correctness, print );
        TestLUMod<Complex<BigFloat>>
        ( grid, m, updateRank, conjugate, tau, correctness, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
