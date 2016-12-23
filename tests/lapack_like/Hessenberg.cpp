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
( UpperOrLower uplo,
  const Matrix<Field>& A,
  const Matrix<Field>& householderScalars,
        Matrix<Field>& AOrig,
  bool print,
  bool display )
{
    typedef Base<Field> Real;
    const Int n = AOrig.Height();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormAOrig = OneNorm( AOrig );
    Output("Testing error...");
    PushIndent();

    // Set H to the appropriate Hessenberg portion of A
    Matrix<Field> H( A );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, H, 1 );
    else
        MakeTrapezoidal( UPPER, H, -1 );
    if( print )
        Print( H, "Hessenberg" );
    if( display )
        Display( H, "Hessenberg" );

    if( print || display )
    {
        Matrix<Field> Q;
        Identity( Q, n, n );
        hessenberg::ApplyQ( LEFT, uplo, NORMAL, A, householderScalars, Q );
        if( print )
            Print( Q, "Q" );
        if( display )
            Display( Q, "Q" );
    }

    // Reverse the accumulated Householder transforms
    hessenberg::ApplyQ( LEFT, uplo, ADJOINT, A, householderScalars, AOrig );
    hessenberg::ApplyQ( RIGHT, uplo, NORMAL, A, householderScalars, AOrig );
    if( print )
        Print( AOrig, "Manual Hessenberg" );
    if( display )
        Display( AOrig, "Manual Hessenberg" );

    // Compare the appropriate portion of AOrig and B
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, AOrig, 1 );
    else
        MakeTrapezoidal( UPPER, AOrig, -1 );
    H -= AOrig;
    if( print )
        Print( H, "Error in rotated Hessenberg" );
    if( display )
        Display( H, "Error in rotated Hessenberg" );
    const Real infError = InfinityNorm( H );
    const Real relError = infError / (n*eps*oneNormAOrig);

    Output("||H - Q^H A Q||_oo / (eps n || A ||_1) = ",relError);

    // TODO(poulson): Use a more refined failure condition
    if( relError > Real(1) )
        LogicError("Unacceptably large relative error");

    PopIndent();
}

template<typename Field>
void TestCorrectness
( UpperOrLower uplo,
  const DistMatrix<Field>& A,
  const DistMatrix<Field,STAR,STAR>& householderScalars,
        DistMatrix<Field>& AOrig,
  bool print,
  bool display )
{
    typedef Base<Field> Real;
    const Grid& grid = A.Grid();
    const Int n = AOrig.Height();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormAOrig = OneNorm( AOrig );
    OutputFromRoot(grid.Comm(),"Testing error...");
    PushIndent();

    // Set H to the appropriate Hessenberg portion of A
    DistMatrix<Field> H( A );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, H, 1 );
    else
        MakeTrapezoidal( UPPER, H, -1 );
    if( print )
        Print( H, "Hessenberg" );
    if( display )
        Display( H, "Hessenberg" );

    if( print || display )
    {
        DistMatrix<Field> Q(grid);
        Identity( Q, n, n );
        hessenberg::ApplyQ( LEFT, uplo, NORMAL, A, householderScalars, Q );
        if( print )
            Print( Q, "Q" );
        if( display )
            Display( Q, "Q" );
    }

    // Reverse the accumulated Householder transforms
    hessenberg::ApplyQ( LEFT, uplo, ADJOINT, A, householderScalars, AOrig );
    hessenberg::ApplyQ( RIGHT, uplo, NORMAL, A, householderScalars, AOrig );
    if( print )
        Print( AOrig, "Manual Hessenberg" );
    if( display )
        Display( AOrig, "Manual Hessenberg" );

    // Compare the appropriate portion of AOrig and B
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, AOrig, 1 );
    else
        MakeTrapezoidal( UPPER, AOrig, -1 );
    H -= AOrig;
    if( print )
        Print( H, "Error in rotated Hessenberg" );
    if( display )
        Display( H, "Error in rotated Hessenberg" );
    const Real infError = InfinityNorm( H );
    const Real relError = infError / (n*eps*oneNormAOrig);

    OutputFromRoot
    (grid.Comm(),"||H - Q^H A Q||_oo / (eps n || A ||_1) = ",relError);

    // TODO(poulson): Use a more refined failure condition
    if( relError > Real(1) )
        LogicError("Unacceptably large relative error");

    PopIndent();
}

template<typename Field>
void TestHessenberg
( UpperOrLower uplo,
  Int n,
  bool correctness,
  bool print,
  bool display )
{
    Matrix<Field> A, AOrig;
    Matrix<Field> householderScalars;
    Output("Testing with ",TypeName<Field>());
    PushIndent();

    Uniform( A, n, n );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    if( display )
        Display( A, "A" );

    Output("Starting reduction to Hessenberg form...");
    Timer timer;
    timer.Start();
    Hessenberg( uplo, A, householderScalars );
    const double runTime = timer.Stop();
    // TODO(poulson): Flop calculation
    Output(runTime," seconds");
    if( print )
    {
        Print( A, "A after Hessenberg" );
        Print( householderScalars, "householderScalars after Hessenberg" );
    }
    if( display )
    {
        Display( A, "A after Hessenberg" );
        Display( householderScalars, "householderScalars after Hessenberg" );
    }
    if( correctness )
        TestCorrectness( uplo, A, householderScalars, AOrig, print, display );
    PopIndent();
}

template<typename Field>
void TestHessenberg
( const Grid& grid,
  UpperOrLower uplo,
  Int n,
  bool correctness,
  bool print,
  bool display )
{
    DistMatrix<Field> A(grid), AOrig(grid);
    DistMatrix<Field,STAR,STAR> householderScalars(grid);
    OutputFromRoot(grid.Comm(),"Testing with ",TypeName<Field>());
    PushIndent();

    Uniform( A, n, n );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    if( display )
        Display( A, "A" );

    OutputFromRoot(grid.Comm(),"Starting reduction to Hessenberg form...");
    mpi::Barrier( grid.Comm() );
    Timer timer;
    timer.Start();
    Hessenberg( uplo, A, householderScalars );
    mpi::Barrier( grid.Comm() );
    const double runTime = timer.Stop();
    // TODO(poulson): Flop calculation
    OutputFromRoot(grid.Comm(),runTime," seconds");
    if( print )
    {
        Print( A, "A after Hessenberg" );
        Print( householderScalars, "householderScalars after Hessenberg" );
    }
    if( display )
    {
        Display( A, "A after Hessenberg" );
        Display( householderScalars, "householderScalars after Hessenberg" );
    }
    if( correctness )
        TestCorrectness( uplo, A, householderScalars, AOrig, print, display );
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
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int n = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool sequential = Input("--sequential","test sequential?",true);
        const bool correctness =
          Input("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool display = Input("--display","display matrices?",false);
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
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        ComplainIfDebug();

        if( sequential && mpi::Rank() == 0 )
        {
            TestHessenberg<float>
            ( uplo, n, correctness, print, display );
            TestHessenberg<Complex<float>>
            ( uplo, n, correctness, print, display );

            TestHessenberg<double>
            ( uplo, n, correctness, print, display );
            TestHessenberg<Complex<double>>
            ( uplo, n, correctness, print, display );

#ifdef EL_HAVE_QD
            TestHessenberg<DoubleDouble>
            ( uplo, n, correctness, print, display );
            TestHessenberg<QuadDouble>
            ( uplo, n, correctness, print, display );

            TestHessenberg<Complex<DoubleDouble>>
            ( uplo, n, correctness, print, display );
            TestHessenberg<Complex<QuadDouble>>
            ( uplo, n, correctness, print, display );
#endif

#ifdef EL_HAVE_QUAD
            TestHessenberg<Quad>
            ( uplo, n, correctness, print, display );
            TestHessenberg<Complex<Quad>>
            ( uplo, n, correctness, print, display );
#endif

#ifdef EL_HAVE_MPC
            TestHessenberg<BigFloat>
            ( uplo, n, correctness, print, display );
            TestHessenberg<Complex<BigFloat>>
            ( uplo, n, correctness, print, display );
#endif
        }

        TestHessenberg<float>
        ( grid, uplo, n, correctness, print, display );
        TestHessenberg<Complex<float>>
        ( grid, uplo, n, correctness, print, display );

        TestHessenberg<double>
        ( grid, uplo, n, correctness, print, display );
        TestHessenberg<Complex<double>>
        ( grid, uplo, n, correctness, print, display );

#ifdef EL_HAVE_QD
        TestHessenberg<DoubleDouble>
        ( grid, uplo, n, correctness, print, display );
        TestHessenberg<QuadDouble>
        ( grid, uplo, n, correctness, print, display );

        TestHessenberg<Complex<DoubleDouble>>
        ( grid, uplo, n, correctness, print, display );
        TestHessenberg<Complex<QuadDouble>>
        ( grid, uplo, n, correctness, print, display );
#endif

#ifdef EL_HAVE_QUAD
        TestHessenberg<Quad>
        ( grid, uplo, n, correctness, print, display );
        TestHessenberg<Complex<Quad>>
        ( grid, uplo, n, correctness, print, display );
#endif

#ifdef EL_HAVE_MPC
        TestHessenberg<BigFloat>
        ( grid, uplo, n, correctness, print, display );
        TestHessenberg<Complex<BigFloat>>
        ( grid, uplo, n, correctness, print, display );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
