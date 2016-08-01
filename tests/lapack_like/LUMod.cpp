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
( bool print, 
  const Matrix<F>& A,
  const Permutation& P,
  const Matrix<F>& AOrig,
        Int numRHS=100 )
{
    typedef Base<F> Real;
    const Int n = AOrig.Width();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( AOrig );

    Output("Testing error...");

    // Generate random right-hand sides
    Matrix<F> X;
    Uniform( X, n, numRHS );
    auto Y( X );
    const Real oneNormY = OneNorm( Y );
    P.PermuteRows( Y );
    lu::SolveAfter( NORMAL, A, Y );

    // Now investigate the residual, ||AOrig Y - X||_oo
    Gemm( NORMAL, NORMAL, F(-1), AOrig, Y, F(1), X );
    const Real infError = InfinityNorm( X );
    const Real relError = infError / (eps*n*Max(oneNormA,oneNormY));

    // TODO: Use a rigorous failure condition
    Output("||A X - Y||_oo / (eps n Max(||A||_1,||Y||_1)) = ",relError);
    if( relError > Real(1000) )
        LogicError("Unacceptably large relative error");
}

template<typename F> 
void TestCorrectness
( bool print, 
  const DistMatrix<F>& A,
  const DistPermutation& P,
  const DistMatrix<F>& AOrig,
        Int numRHS=100 )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = AOrig.Width();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( AOrig );

    OutputFromRoot(g.Comm(),"Testing error...");

    // Generate random right-hand sides
    DistMatrix<F> X(g);
    Uniform( X, n, numRHS );
    auto Y( X );
    const Real oneNormY = OneNorm( Y );
    P.PermuteRows( Y );
    lu::SolveAfter( NORMAL, A, Y );

    // Now investigate the residual, ||AOrig Y - X||_oo
    Gemm( NORMAL, NORMAL, F(-1), AOrig, Y, F(1), X );
    const Real infError = InfinityNorm( X );
    const Real relError = infError / (eps*n*Max(oneNormA,oneNormY));

    // TODO: Use a rigorous failure condition
    OutputFromRoot
    (g.Comm(),"||A X - Y||_oo / (eps n Max(||A||_1,||Y||_1)) = ",relError);
    if( relError > Real(1000) )
        LogicError("Unacceptably large relative error");
}

template<typename F> 
void TestLUMod
( Int m,
  bool conjugate,
  Base<F> tau,
  bool correctness,
  bool print )
{
    Output("Testing with ",TypeName<F>());
    PushIndent();
    Matrix<F> A, AOrig;
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
        P.ReserveSwaps( m+2*m-1 );
        LU( A, P );
        const double runTime = timer.Stop();
        const double realGFlops = 2./3.*Pow(double(m),3.)/(1.e9*runTime);
        const double gFlops =
          ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
        Output(runTime," seconds (",gFlops," GFlop/s)");
    }

    // TODO: Print permutation

    // Generate random vectors u and v
    Matrix<F> u, v;
    Uniform( u, m, 1 );
    Uniform( v, m, 1 );
    if( correctness )
    {
        if( conjugate )
            Ger( F(1), u, v, AOrig );
        else
            Geru( F(1), u, v, AOrig );
    }

    { 
        Output("Starting rank-one LU modification...");
        Timer timer;
        timer.Start();
        LUMod( A, P, u, v, conjugate, tau );
        const double runTime = timer.Stop();
        Output(runTime," seconds");
    }

    // TODO: Print permutation

    if( correctness )
        TestCorrectness( print, A, P, AOrig );
    PopIndent();
}

template<typename F> 
void TestLUMod
( const Grid& g,
  Int m,
  bool conjugate,
  Base<F> tau,
  bool correctness,
  bool print )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();
    DistMatrix<F> A(g), AOrig(g);
    DistPermutation P(g);

    Uniform( A, m, m );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    {
        OutputFromRoot(g.Comm(),"Starting LU factorization...");
        mpi::Barrier( g.Comm() );
        Timer timer;
        timer.Start();
        P.ReserveSwaps( m+2*m-1 );
        LU( A, P );
        mpi::Barrier( g.Comm() );
        const double runTime = timer.Stop();
        const double realGFlops = 2./3.*Pow(double(m),3.)/(1.e9*runTime);
        const double gFlops =
          ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
        OutputFromRoot(g.Comm(),runTime," seconds (",gFlops," GFlop/s)");
    }

    // TODO: Print permutation

    // Generate random vectors u and v
    DistMatrix<F> u(g), v(g);
    Uniform( u, m, 1 );
    Uniform( v, m, 1 );
    if( correctness )
    {
        if( conjugate )
            Ger( F(1), u, v, AOrig );
        else
            Geru( F(1), u, v, AOrig );
    }

    { 
        OutputFromRoot(g.Comm(),"Starting rank-one LU modification...");
        mpi::Barrier( g.Comm() );
        Timer timer;
        timer.Start();
        LUMod( A, P, u, v, conjugate, tau );
        mpi::Barrier( g.Comm() );
        const double runTime = timer.Stop();
        OutputFromRoot(g.Comm(),runTime," seconds");
    }

    // TODO: Print permutation

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
            gridHeight = Grid::FindFactor( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );
        SetBlocksize( nb );
        ComplainIfDebug();

        if( sequential && mpi::Rank() == 0 )
        {
            TestLUMod<float>
            ( m, conjugate, tau, correctness, print );
            TestLUMod<Complex<float>>
            ( m, conjugate, tau, correctness, print );

            TestLUMod<double>
            ( m, conjugate, tau, correctness, print );
            TestLUMod<Complex<double>>
            ( m, conjugate, tau, correctness, print );

#ifdef EL_HAVE_QD
            TestLUMod<DoubleDouble>
            ( m, conjugate, tau, correctness, print );
            TestLUMod<QuadDouble>
            ( m, conjugate, tau, correctness, print );

            TestLUMod<Complex<DoubleDouble>>
            ( m, conjugate, tau, correctness, print );
            TestLUMod<Complex<QuadDouble>>
            ( m, conjugate, tau, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
            TestLUMod<Quad>
            ( m, conjugate, tau, correctness, print );
            TestLUMod<Complex<Quad>>
            ( m, conjugate, tau, correctness, print );
#endif

#ifdef EL_HAVE_MPC
            TestLUMod<BigFloat>
            ( m, conjugate, tau, correctness, print );
            TestLUMod<Complex<BigFloat>>
            ( m, conjugate, tau, correctness, print );
#endif
        }

        TestLUMod<float>
        ( g, m, conjugate, tau, correctness, print );
        TestLUMod<Complex<float>>
        ( g, m, conjugate, tau, correctness, print );

        TestLUMod<double>
        ( g, m, conjugate, tau, correctness, print );
        TestLUMod<Complex<double>>
        ( g, m, conjugate, tau, correctness, print );

#ifdef EL_HAVE_QD
        TestLUMod<DoubleDouble>
        ( g, m, conjugate, tau, correctness, print );
        TestLUMod<QuadDouble>
        ( g, m, conjugate, tau, correctness, print );

        TestLUMod<Complex<DoubleDouble>>
        ( g, m, conjugate, tau, correctness, print );
        TestLUMod<Complex<QuadDouble>>
        ( g, m, conjugate, tau, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
        TestLUMod<Quad>
        ( g, m, conjugate, tau, correctness, print );
        TestLUMod<Complex<Quad>>
        ( g, m, conjugate, tau, correctness, print );
#endif

#ifdef EL_HAVE_MPC
        TestLUMod<BigFloat>
        ( g, m, conjugate, tau, correctness, print );
        TestLUMod<Complex<BigFloat>>
        ( g, m, conjugate, tau, correctness, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
