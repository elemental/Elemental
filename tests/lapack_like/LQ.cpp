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
  const Matrix<F>& t,
  const Matrix<Base<F>>& d,
        Matrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = std::min(m,n);
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( AOrig );

    Output("Testing orthogonality of Q...");
    PushIndent();

    // Form Z := Q Q^H as an approximation to identity
    Matrix<F> Z;
    Identity( Z, m, n );
    lq::ApplyQ( RIGHT, NORMAL, A, t, d, Z );
    lq::ApplyQ( RIGHT, ADJOINT, A, t, d, Z );
    auto ZUpper = View( Z, 0, 0, minDim, minDim );

    // Form X := I - Q Q^H
    Matrix<F> X;
    Identity( X, minDim, minDim );
    X -= ZUpper;

    const Real infOrthogError = InfinityNorm( X );
    const Real relOrthogError = infOrthogError / (eps*Max(m,n));
    Output("||Q Q^H - I||_oo / (eps Max(m,n)) = ",relOrthogError);
    PopIndent();

    Output("Testing if A = LQ...");
    PushIndent();

    // Form L Q
    auto L( A );
    MakeTrapezoidal( LOWER, L );
    lq::ApplyQ( RIGHT, NORMAL, A, t, d, L );

    // Form L Q - A
    L -= AOrig;

    const Real infError = InfinityNorm( L );
    const Real relError = infError / (eps*Max(m,n)*oneNormA);
    Output("||A - LQ||_oo / (eps Max(m,n) ||A||_1) = ",relError);
    PopIndent();

    // TODO: More rigorous failure condition
    if( relOrthogError > Real(10) )
        LogicError("Relative orthogonality error was unacceptably large");
    if( relError > Real(10) ) 
        LogicError("Relative error was unacceptably large");
}

template<typename F> 
void TestCorrectness
( bool print,
  const DistMatrix<F>& A,
  const DistMatrix<F,MD,STAR>& t,
  const DistMatrix<Base<F>,MD,STAR>& d,
        DistMatrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = std::min(m,n);
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( AOrig );

    OutputFromRoot(g.Comm(),"Testing orthogonality of Q...");
    PushIndent();

    // Form Z := Q Q^H as an approximation to identity
    DistMatrix<F> Z(g);
    Identity( Z, m, n );
    lq::ApplyQ( RIGHT, NORMAL, A, t, d, Z );
    lq::ApplyQ( RIGHT, ADJOINT, A, t, d, Z );
    auto ZUpper = View( Z, 0, 0, minDim, minDim );

    // Form X := I - Q Q^H
    DistMatrix<F> X(g);
    Identity( X, minDim, minDim );
    X -= ZUpper;

    const Real infOrthogError = InfinityNorm( X );
    const Real relOrthogError = infOrthogError / (eps*Max(m,n));
    OutputFromRoot
    (g.Comm(),"||Q Q^H - I||_oo / (eps Max(m,n)) = ",relOrthogError);
    PopIndent();

    OutputFromRoot(g.Comm(),"Testing if A = LQ...");
    PushIndent();

    // Form L Q
    auto L( A );
    MakeTrapezoidal( LOWER, L );
    lq::ApplyQ( RIGHT, NORMAL, A, t, d, L );

    // Form L Q - A
    L -= AOrig;

    const Real infError = InfinityNorm( L );
    const Real relError = infError / (eps*Max(m,n)*oneNormA);
    OutputFromRoot
    (g.Comm(),"||A - LQ||_oo / (eps Max(m,n) ||A||_1) = ",relError);
    PopIndent();

    // TODO: More rigorous failure condition
    if( relOrthogError > Real(10) )
        LogicError("Relative orthogonality error was unacceptably large");
    if( relError > Real(10) ) 
        LogicError("Relative error was unacceptably large");
}

template<typename F>
void TestLQ( Int m, Int n, bool correctness, bool print )
{
    Output("Testing with ",TypeName<F>());
    PushIndent();
    Matrix<F> A, AOrig;
    Uniform( A, m, n );

    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    Matrix<F> t;
    Matrix<Base<F>> d;

    Output("Starting LQ factorization...");
    Timer timer;
    timer.Start();
    LQ( A, t, d );
    const double runTime = timer.Stop();
    const double mD = double(m);
    const double nD = double(n);
    const double realGFlops = (2.*mD*mD*nD - 2./3.*mD*mD*mD)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    Output(runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        Print( t, "phases" );
        Print( d, "diagonal" );
    }
    if( correctness )
        TestCorrectness( print, A, t, d, AOrig );
    PopIndent();
}

template<typename F>
void TestLQ( const Grid& g, Int m, Int n, bool correctness, bool print )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();
    DistMatrix<F> A(g), AOrig(g);
    Uniform( A, m, n );

    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);

    OutputFromRoot(g.Comm(),"Starting LQ factorization...");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    LQ( A, t, d );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    const double mD = double(m);
    const double nD = double(n);
    const double realGFlops = (2.*mD*mD*nD - 2./3.*mD*mD*mD)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    OutputFromRoot(g.Comm(),runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        Print( t, "phases" );
        Print( d, "diagonal" );
    }
    if( correctness )
        TestCorrectness( print, A, t, d, AOrig );
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
        const Int n = Input("--width","width of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
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
        mpc::SetPrecision( prec );
#endif

        if( gridHeight == 0 )
            gridHeight = Grid::FindFactor( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );
        SetBlocksize( nb );
        ComplainIfDebug();

        if( sequential && mpi::Rank() == 0 )
        {
            TestLQ<float>
            ( m, n, correctness, print );
            TestLQ<Complex<float>>
            ( m, n, correctness, print );

            TestLQ<double>
            ( m, n, correctness, print );
            TestLQ<Complex<double>>
            ( m, n, correctness, print );

#ifdef EL_HAVE_QD
            TestLQ<DoubleDouble>
            ( m, n, correctness, print );
            TestLQ<QuadDouble>
            ( m, n, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
            TestLQ<Quad>
            ( m, n, correctness, print );
            TestLQ<Complex<Quad>>
            ( m, n, correctness, print );
#endif

#ifdef EL_HAVE_MPC
            TestLQ<BigFloat>
            ( m, n, correctness, print );
#endif
        }

        TestLQ<float>
        ( g, m, n, correctness, print );
        TestLQ<Complex<float>>
        ( g, m, n, correctness, print );

        TestLQ<double>
        ( g, m, n, correctness, print );
        TestLQ<Complex<double>>
        ( g, m, n, correctness, print );

#ifdef EL_HAVE_QD
        TestLQ<DoubleDouble>
        ( g, m, n, correctness, print );
        TestLQ<QuadDouble>
        ( g, m, n, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
        TestLQ<Quad>
        ( g, m, n, correctness, print );
        TestLQ<Complex<Quad>>
        ( g, m, n, correctness, print );
#endif

#ifdef EL_HAVE_MPC
        TestLQ<BigFloat>
        ( g, m, n, correctness, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
