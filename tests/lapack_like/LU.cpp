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
( const Matrix<F>& AOrig,
  const Matrix<F>& A,
  const Permutation& P,
  const Permutation& Q,
  Int pivoting,
  bool print,
  Int numRHS=100 )
{
    typedef Base<F> Real;
    const Int m = AOrig.Height();
    const Int n = AOrig.Width();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormAOrig = OneNorm( AOrig );

    Output("Testing error...");
    PushIndent();

    // Generate random right-hand sides
    Matrix<F> X;
    Uniform( X, m, numRHS );
    auto Y( X );
    const Real oneNormY = OneNorm( Y );
    if( pivoting == 0 )
        lu::SolveAfter( NORMAL, A, Y );
    else if( pivoting == 1 )
        lu::SolveAfter( NORMAL, A, P, Y );
    else
        lu::SolveAfter( NORMAL, A, P, Q, Y );

    // Now investigate the residual, ||AOrig Y - X||_oo
    Gemm( NORMAL, NORMAL, F(-1), AOrig, Y, F(1), X );
    const Real infError = InfinityNorm( X );
    const Real relError = infError / (eps*Max(m,n)*Max(oneNormAOrig,oneNormY));

    Output
    ("|| Y - A X ||_oo / (eps Max(m,n) Max(||A||_1,||Y||_1)) = ",relError);

    // TODO: Use a more refined failure condition
    if( pivoting > 0 && relError > Real(100) )
        LogicError("Relative error was unacceptably large");

    PopIndent();
}

template<typename F> 
void TestCorrectness
( const DistMatrix<F>& AOrig,
  const DistMatrix<F>& A,
  const DistPermutation& P,
  const DistPermutation& Q,
  Int pivoting,
  bool print,
  Int numRHS=100 )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = AOrig.Height();
    const Int n = AOrig.Width();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormAOrig = OneNorm( AOrig );

    OutputFromRoot(g.Comm(),"Testing error...");
    PushIndent();

    // Generate random right-hand sides
    DistMatrix<F> X(g);
    Uniform( X, m, numRHS );
    auto Y( X );
    const Real oneNormY = OneNorm( Y );
    if( pivoting == 0 )
        lu::SolveAfter( NORMAL, A, Y );
    else if( pivoting == 1 )
        lu::SolveAfter( NORMAL, A, P, Y );
    else
        lu::SolveAfter( NORMAL, A, P, Q, Y );

    // Now investigate the residual, ||AOrig Y - X||_oo
    Gemm( NORMAL, NORMAL, F(-1), AOrig, Y, F(1), X );
    const Real infError = InfinityNorm( X );
    const Real relError = infError / (eps*Max(m,n)*Max(oneNormAOrig,oneNormY));

    OutputFromRoot
    (g.Comm(),
     "|| Y - A X ||_oo / (eps Max(m,n) Max(||A||_1,||Y||_1)) = ",relError);

    // TODO: Use a more refined failure condition
    if( pivoting > 0 && relError > Real(100) )
        LogicError("Relative error was unacceptably large");

    PopIndent();
}

template<typename F> 
void TestLU
( Int m,
  Int pivoting, 
  bool correctness,
  bool forceGrowth,
  bool print )
{
    Output("Testing with ",TypeName<F>());
    PushIndent();
    Matrix<F> A, AOrig;
    Permutation P, Q;

    if( forceGrowth )
        GEPPGrowth( A, m );
    else
        Uniform( A, m, m );

    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    Output("Starting LU factorization...");
    Timer timer;
    timer.Start();
    if( pivoting == 0 )
        LU( A );
    else if( pivoting == 1 )
        LU( A, P );
    else if( pivoting == 2 )
        LU( A, P, Q );
    const double runTime = timer.Stop();
    const double realGFlops = 2./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    Output(runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        if( pivoting >= 1 )
        {
            // This needs to be rewritten in light of Permutation
            /*
            Matrix<Int> P;
            Matrix<Int> p;
            Print( P, "rowPiv after factorization" );
            PivotsToPermutation( rowPiv, p );
            Print( p, "p after factorization");
            ExplicitPermutation( p, P );
            Print( P, "P" );
            */
        }
        if( pivoting == 2 )
        {
            // This needs to be rewritten in light of Permutation
            /*
            Matrix<Int> Q;
            Matrix<Int> q;
            Print( colPiv, "colPiv after factorization" );
            PivotsToPermutation( colPiv, q );
            Print( q, "q after factorization");
            ExplicitPermutation( q, Q );
            Print( Q, "Q" );
            */
        }
    }
    if( correctness )
        TestCorrectness( AOrig, A, P, Q, pivoting, print );
    PopIndent();
}

template<typename F> 
void TestLU
( const Grid& g,
  Int m,
  Int pivoting, 
  bool correctness,
  bool forceGrowth,
  bool print )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();
    DistMatrix<F> A(g), AOrig(g);
    DistPermutation P(g), Q(g);

    if( forceGrowth )
        GEPPGrowth( A, m );
    else
        Uniform( A, m, m );

    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    OutputFromRoot(g.Comm(),"Starting LU factorization...");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    if( pivoting == 0 )
        LU( A );
    else if( pivoting == 1 )
        LU( A, P );
    else if( pivoting == 2 )
        LU( A, P, Q );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    const double realGFlops = 2./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    OutputFromRoot(g.Comm(),runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        if( pivoting >= 1 )
        {
            // This needs to be rewritten in light of DistPermutation
            /*
            DistMatrix<Int> P(g);
            DistMatrix<Int,STAR,STAR> p(g);
            Print( P, "rowPiv after factorization" );
            PivotsToPermutation( rowPiv, p );
            Print( p, "p after factorization");
            ExplicitPermutation( p, P );
            Print( P, "P" );
            */
        }
        if( pivoting == 2 )
        {
            // This needs to be rewritten in light of DistPermutation
            /*
            DistMatrix<Int> Q(g);
            DistMatrix<Int,STAR,STAR> q(g);
            Print( colPiv, "colPiv after factorization" );
            PivotsToPermutation( colPiv, q );
            Print( q, "q after factorization");
            ExplicitPermutation( q, Q );
            Print( Q, "Q" );
            */
        }
    }
    if( correctness )
        TestCorrectness( AOrig, A, P, Q, pivoting, print );
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
        const Int pivot = Input("--pivot","0: none, 1: partial, 2: full",1);
        const bool forceGrowth = Input
            ("--forceGrowth","force element growth?",false);
        const bool sequential = Input("--sequential","test sequential?",true);
        const bool correctness = 
          Input("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        ProcessInput();
        PrintInputReport();
        if( pivot < 0 || pivot > 2 )
            LogicError("Invalid pivot value");

#ifdef EL_HAVE_MPC
        mpfr::SetPrecision( prec );
#endif

        if( gridHeight == 0 )
            gridHeight = Grid::FindFactor( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );
        SetBlocksize( nb );
        ComplainIfDebug();
        if( pivot == 0 )
            OutputFromRoot(g.Comm(),"Testing LU with no pivoting");
        else if( pivot == 1 )
            OutputFromRoot(g.Comm(),"Testing LU with partial pivoting");
        else if( pivot == 2 )
            OutputFromRoot(g.Comm(),"Testing LU with full pivoting");

        if( sequential && mpi::Rank() == 0 )
        {
            TestLU<float>
            ( m, pivot, correctness, forceGrowth, print );
            TestLU<Complex<float>>
            ( m, pivot, correctness, forceGrowth, print );

            TestLU<double>
            ( m, pivot, correctness, forceGrowth, print );
            TestLU<Complex<double>>
            ( m, pivot, correctness, forceGrowth, print );

#ifdef EL_HAVE_QD
            TestLU<DoubleDouble>
            ( m, pivot, correctness, forceGrowth, print );
            TestLU<QuadDouble>
            ( m, pivot, correctness, forceGrowth, print );

            TestLU<Complex<DoubleDouble>>
            ( m, pivot, correctness, forceGrowth, print );
            TestLU<Complex<QuadDouble>>
            ( m, pivot, correctness, forceGrowth, print );
#endif

#ifdef EL_HAVE_QUAD
            TestLU<Quad>
            ( m, pivot, correctness, forceGrowth, print );
            TestLU<Complex<Quad>>
            ( m, pivot, correctness, forceGrowth, print );
#endif

#ifdef EL_HAVE_MPC
            TestLU<BigFloat>
            ( m, pivot, correctness, forceGrowth, print );
            TestLU<Complex<BigFloat>>
            ( m, pivot, correctness, forceGrowth, print );
#endif
        }

        TestLU<float>
        ( g, m, pivot, correctness, forceGrowth, print );
        TestLU<Complex<float>>
        ( g, m, pivot, correctness, forceGrowth, print );

        TestLU<double>
        ( g, m, pivot, correctness, forceGrowth, print );
        TestLU<Complex<double>>
        ( g, m, pivot, correctness, forceGrowth, print );

#ifdef EL_HAVE_QD
        TestLU<DoubleDouble>
        ( g, m, pivot, correctness, forceGrowth, print );
        TestLU<QuadDouble>
        ( g, m, pivot, correctness, forceGrowth, print );

        TestLU<Complex<DoubleDouble>>
        ( g, m, pivot, correctness, forceGrowth, print );
        TestLU<Complex<QuadDouble>>
        ( g, m, pivot, correctness, forceGrowth, print );
#endif

#ifdef EL_HAVE_QUAD
        TestLU<Quad>
        ( g, m, pivot, correctness, forceGrowth, print );
        TestLU<Complex<Quad>>
        ( g, m, pivot, correctness, forceGrowth, print );
#endif

#ifdef EL_HAVE_MPC
        TestLU<BigFloat>
        ( g, m, pivot, correctness, forceGrowth, print );
        TestLU<Complex<BigFloat>>
        ( g, m, pivot, correctness, forceGrowth, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
