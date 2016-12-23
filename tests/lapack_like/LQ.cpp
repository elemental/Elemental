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
  const Matrix<Field>& householderScalars,
  const Matrix<Base<Field>>& signature,
        Matrix<Field>& AOrig )
{
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = std::min(m,n);
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( AOrig );

    Output("Testing orthogonality of Q...");
    PushIndent();

    // Form Z := Q Q^H as an approximation to identity
    Matrix<Field> Z;
    Identity( Z, m, n );
    lq::ApplyQ( RIGHT, NORMAL, A, householderScalars, signature, Z );
    lq::ApplyQ( RIGHT, ADJOINT, A, householderScalars, signature, Z );
    auto ZUpper = Z( IR(0,minDim), IR(0,minDim) );

    // Form X := I - Q Q^H
    Matrix<Field> X;
    Identity( X, minDim, minDim );
    X -= ZUpper;

    const Real infOrthogError = InfinityNorm( X );
    const Real relOrthogError = infOrthogError / (eps*Max(m,n));
    Output("||Q Q^H - I||_oo / (eps Max(m,n)) = ",relOrthogError);
    PopIndent();

    Output("Testing if A = LQ...");
    PushIndent();

    // Form L Q - A
    auto L( A );
    MakeTrapezoidal( LOWER, L );
    lq::ApplyQ( RIGHT, NORMAL, A, householderScalars, signature, L );
    L -= AOrig;

    const Real infError = InfinityNorm( L );
    const Real relError = infError / (eps*Max(m,n)*oneNormA);
    Output("||A - LQ||_oo / (eps Max(m,n) ||A||_1) = ",relError);
    PopIndent();

    // TODO(poulson): More rigorous failure condition
    if( relOrthogError > Real(10) )
        LogicError("Relative orthogonality error was unacceptably large");
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large");
}

template<typename Field>
void TestCorrectness
( bool print,
  const DistMatrix<Field>& A,
  const DistMatrix<Field,MD,STAR>& householderScalars,
  const DistMatrix<Base<Field>,MD,STAR>& signature,
        DistMatrix<Field>& AOrig )
{
    typedef Base<Field> Real;
    const Grid& grid = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = std::min(m,n);
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( AOrig );

    OutputFromRoot(grid.Comm(),"Testing orthogonality of Q...");
    PushIndent();

    // Form Z := Q Q^H as an approximation to identity
    DistMatrix<Field> Z(grid);
    Identity( Z, m, n );
    lq::ApplyQ( RIGHT, NORMAL, A, householderScalars, signature, Z );
    lq::ApplyQ( RIGHT, ADJOINT, A, householderScalars, signature, Z );
    auto ZUpper = Z( IR(0,minDim), IR(0,minDim) );

    // Form X := I - Q Q^H
    DistMatrix<Field> X(grid);
    Identity( X, minDim, minDim );
    X -= ZUpper;

    const Real infOrthogError = InfinityNorm( X );
    const Real relOrthogError = infOrthogError / (eps*Max(m,n));
    OutputFromRoot
    (grid.Comm(),"||Q Q^H - I||_oo / (eps Max(m,n)) = ",relOrthogError);
    PopIndent();

    OutputFromRoot(grid.Comm(),"Testing if A = LQ...");
    PushIndent();

    // Form L Q - A
    auto L( A );
    MakeTrapezoidal( LOWER, L );
    lq::ApplyQ( RIGHT, NORMAL, A, householderScalars, signature, L );
    L -= AOrig;

    const Real infError = InfinityNorm( L );
    const Real relError = infError / (eps*Max(m,n)*oneNormA);
    OutputFromRoot
    (grid.Comm(),"||A - LQ||_oo / (eps Max(m,n) ||A||_1) = ",relError);
    PopIndent();

    // TODO(poulson): More rigorous failure condition
    if( relOrthogError > Real(10) )
        LogicError("Relative orthogonality error was unacceptably large");
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large");
}

template<typename Field>
void TestLQ( Int m, Int n, bool correctness, bool print )
{
    Output("Testing with ",TypeName<Field>());
    PushIndent();
    Matrix<Field> A, AOrig;
    Uniform( A, m, n );

    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    Matrix<Field> householderScalars;
    Matrix<Base<Field>> signature;

    Output("Starting LQ factorization...");
    Timer timer;
    timer.Start();
    LQ( A, householderScalars, signature );
    const double runTime = timer.Stop();
    const double mD = double(m);
    const double nD = double(n);
    const double realGFlops = (2.*mD*mD*nD - 2./3.*mD*mD*mD)/(1.e9*runTime);
    const double gFlops = IsComplex<Field>::value ? 4*realGFlops : realGFlops;
    Output(runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        Print( householderScalars, "householderScalars" );
        Print( signature, "signature" );
    }
    if( correctness )
        TestCorrectness( print, A, householderScalars, signature, AOrig );
    PopIndent();
}

template<typename Field>
void TestLQ( const Grid& grid, Int m, Int n, bool correctness, bool print )
{
    OutputFromRoot(grid.Comm(),"Testing with ",TypeName<Field>());
    PushIndent();
    DistMatrix<Field> A(grid), AOrig(grid);
    Uniform( A, m, n );

    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    DistMatrix<Field,MD,STAR> householderScalars(grid);
    DistMatrix<Base<Field>,MD,STAR> signature(grid);

    OutputFromRoot(grid.Comm(),"Starting LQ factorization...");
    mpi::Barrier( grid.Comm() );
    Timer timer;
    timer.Start();
    LQ( A, householderScalars, signature );
    mpi::Barrier( grid.Comm() );
    const double runTime = timer.Stop();
    const double mD = double(m);
    const double nD = double(n);
    const double realGFlops = (2.*mD*mD*nD - 2./3.*mD*mD*mD)/(1.e9*runTime);
    const double gFlops = IsComplex<Field>::value ? 4*realGFlops : realGFlops;
    OutputFromRoot(grid.Comm(),runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        Print( householderScalars, "householderScalars" );
        Print( signature, "signature" );
    }
    if( correctness )
        TestCorrectness( print, A, householderScalars, signature, AOrig );
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

            TestLQ<Complex<DoubleDouble>>
            ( m, n, correctness, print );
            TestLQ<Complex<QuadDouble>>
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
            TestLQ<Complex<BigFloat>>
            ( m, n, correctness, print );
#endif
        }

        TestLQ<float>
        ( grid, m, n, correctness, print );
        TestLQ<Complex<float>>
        ( grid, m, n, correctness, print );

        TestLQ<double>
        ( grid, m, n, correctness, print );
        TestLQ<Complex<double>>
        ( grid, m, n, correctness, print );

#ifdef EL_HAVE_QD
        TestLQ<DoubleDouble>
        ( grid, m, n, correctness, print );
        TestLQ<QuadDouble>
        ( grid, m, n, correctness, print );

        TestLQ<Complex<DoubleDouble>>
        ( grid, m, n, correctness, print );
        TestLQ<Complex<QuadDouble>>
        ( grid, m, n, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
        TestLQ<Quad>
        ( grid, m, n, correctness, print );
        TestLQ<Complex<Quad>>
        ( grid, m, n, correctness, print );
#endif

#ifdef EL_HAVE_MPC
        TestLQ<BigFloat>
        ( grid, m, n, correctness, print );
        TestLQ<Complex<BigFloat>>
        ( grid, m, n, correctness, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
