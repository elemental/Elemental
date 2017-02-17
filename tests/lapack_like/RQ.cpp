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
( const Matrix<Field>& A,
  const Matrix<Field>& householderScalars,
  const Matrix<Base<Field>>& signature,
        Matrix<Field>& AOrig )
{
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int maxDim = Max(m,n);
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( AOrig );

    Output("Testing orthogonality of Q...");
    PushIndent();

    // Form Z := Q Q^H as an approximation to identity
    Matrix<Field> Z;
    Identity( Z, m, n );
    rq::ApplyQ( RIGHT, NORMAL, A, householderScalars, signature, Z );
    rq::ApplyQ( RIGHT, ADJOINT, A, householderScalars, signature, Z );
    auto ZUpper = Z( IR(0,minDim), IR(0,minDim) );

    // Form X := I - Q Q^H
    Matrix<Field> X;
    Identity( X, minDim, minDim );
    X -= ZUpper;

    const Real infOrthogError = InfinityNorm( X );
    const Real relOrthogError = infOrthogError / (eps*maxDim);
    Output("||Q^H Q - I||_oo / (eps Max(m,n)) = ",relOrthogError);
    PopIndent();

    Output("Testing if A = RQ...");
    PushIndent();

    // Form RQ - A
    auto U( A );
    MakeTrapezoidal( UPPER, U, U.Width()-U.Height() );
    rq::ApplyQ( RIGHT, NORMAL, A, householderScalars, signature, U );
    U -= AOrig;

    const Real infError = InfinityNorm( U );
    const Real relError = infError / (eps*maxDim*oneNormA);
    Output("||A - RQ||_oo / (eps Max(m,n) ||A||_1)= ",relError);

    PopIndent();

    // TODO(poulson): More rigorous failure conditions
    if( relOrthogError > Real(10) )
        LogicError("Unacceptably large relative orthogonality error");
    if( relError > Real(10) )
        LogicError("Unacceptably large relative error");
}

template<typename Field>
void TestCorrectness
( const DistMatrix<Field>& A,
  const DistMatrix<Field,MD,STAR>& householderScalars,
  const DistMatrix<Base<Field>,MD,STAR>& signature,
        DistMatrix<Field>& AOrig )
{
    typedef Base<Field> Real;
    const Grid& grid = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int maxDim = Max(m,n);
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( AOrig );

    OutputFromRoot(grid.Comm(),"Testing orthogonality of Q...");
    PushIndent();

    // Form Z := Q Q^H as an approximation to identity
    DistMatrix<Field> Z(grid);
    Identity( Z, m, n );
    rq::ApplyQ( RIGHT, NORMAL, A, householderScalars, signature, Z );
    rq::ApplyQ( RIGHT, ADJOINT, A, householderScalars, signature, Z );
    auto ZUpper = Z( IR(0,minDim), IR(0,minDim) );

    // Form X := I - Q Q^H
    DistMatrix<Field> X(grid);
    Identity( X, minDim, minDim );
    X -= ZUpper;

    const Real infOrthogError = InfinityNorm( X );
    const Real relOrthogError = infOrthogError / (eps*maxDim);
    OutputFromRoot
    (grid.Comm(),"||Q^H Q - I||_oo / (eps Max(m,n)) = ",relOrthogError);
    PopIndent();

    OutputFromRoot(grid.Comm(),"Testing if A = RQ...");
    PushIndent();

    // Form RQ - A
    auto U( A );
    MakeTrapezoidal( UPPER, U, U.Width()-U.Height() );
    rq::ApplyQ( RIGHT, NORMAL, A, householderScalars, signature, U );
    U -= AOrig;

    const Real infError = InfinityNorm( U );
    const Real relError = infError / (eps*maxDim*oneNormA);
    OutputFromRoot
    (grid.Comm(),"||A - RQ||_oo / (eps Max(m,n) ||A||_1)= ",relError);

    PopIndent();

    // TODO(poulson): More rigorous failure conditions
    if( relOrthogError > Real(10) )
        LogicError("Unacceptably large relative orthogonality error");
    if( relError > Real(10) )
        LogicError("Unacceptably large relative error");
}

template<typename Field>
void TestRQ
( Int m,
  Int n,
  bool correctness,
  bool print )
{
    Output("Testing with ",TypeName<Field>());
    PushIndent();
    Matrix<Field> A, AOrig;
    Matrix<Field> householderScalars;
    Matrix<Base<Field>> signature;

    Uniform( A, m, n );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    Output("Starting RQ factorization...");
    Timer timer;
    timer.Start();
    RQ( A, householderScalars, signature );
    const double runTime = timer.Stop();
    const double mD = double(m);
    const double nD = double(n);
    const double realGFlops = (2.*mD*nD*nD - 2./3.*nD*nD*nD)/(1.e9*runTime);
    const double gFlops = IsComplex<Field>::value ? 4*realGFlops : realGFlops;
    Output(runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        Print( householderScalars, "householderScalars" );
        Print( signature, "signature" );
    }
    if( correctness )
        TestCorrectness( A, householderScalars, signature, AOrig );
    PopIndent();
}

template<typename Field>
void TestRQ
( const Grid& g,
  Int m,
  Int n,
  bool correctness,
  bool print )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<Field>());
    PushIndent();
    DistMatrix<Field> A(g), AOrig(g);
    DistMatrix<Field,MD,STAR> householderScalars(g);
    DistMatrix<Base<Field>,MD,STAR> signature(g);

    Uniform( A, m, n );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    OutputFromRoot(g.Comm(),"Starting RQ factorization...");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    RQ( A, householderScalars, signature );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    const double mD = double(m);
    const double nD = double(n);
    const double realGFlops = (2.*mD*nD*nD - 2./3.*nD*nD*nD)/(1.e9*runTime);
    const double gFlops = IsComplex<Field>::value ? 4*realGFlops : realGFlops;
    OutputFromRoot(g.Comm(),runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        Print( householderScalars, "householderScalars" );
        Print( signature, "signature" );
    }
    if( correctness )
        TestCorrectness( A, householderScalars, signature, AOrig );
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
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid grid( comm, gridHeight, order );
        SetBlocksize( nb );
        ComplainIfDebug();

        if( sequential && mpi::Rank() == 0 )
        {
            TestRQ<float>
            ( m, n, correctness, print );
            TestRQ<Complex<float>>
            ( m, n, correctness, print );

            TestRQ<double>
            ( m, n, correctness, print );
            TestRQ<Complex<double>>
            ( m, n, correctness, print );

#ifdef EL_HAVE_QD
            TestRQ<DoubleDouble>
            ( m, n, correctness, print );
            TestRQ<QuadDouble>
            ( m, n, correctness, print );

            TestRQ<Complex<DoubleDouble>>
            ( m, n, correctness, print );
            TestRQ<Complex<QuadDouble>>
            ( m, n, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
            TestRQ<Quad>
            ( m, n, correctness, print );
            TestRQ<Complex<Quad>>
            ( m, n, correctness, print );
#endif

#ifdef EL_HAVE_MPC
            TestRQ<BigFloat>
            ( m, n, correctness, print );
            TestRQ<Complex<BigFloat>>
            ( m, n, correctness, print );
#endif
        }

        TestRQ<float>
        ( grid, m, n, correctness, print );
        TestRQ<Complex<float>>
        ( grid, m, n, correctness, print );

        TestRQ<double>
        ( grid, m, n, correctness, print );
        TestRQ<Complex<double>>
        ( grid, m, n, correctness, print );

#ifdef EL_HAVE_QD
        TestRQ<DoubleDouble>
        ( grid, m, n, correctness, print );
        TestRQ<QuadDouble>
        ( grid, m, n, correctness, print );

        TestRQ<Complex<DoubleDouble>>
        ( grid, m, n, correctness, print );
        TestRQ<Complex<QuadDouble>>
        ( grid, m, n, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
        TestRQ<Quad>
        ( grid, m, n, correctness, print );
        TestRQ<Complex<Quad>>
        ( grid, m, n, correctness, print );
#endif

#ifdef EL_HAVE_MPC
        TestRQ<BigFloat>
        ( grid, m, n, correctness, print );
        TestRQ<Complex<BigFloat>>
        ( grid, m, n, correctness, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
