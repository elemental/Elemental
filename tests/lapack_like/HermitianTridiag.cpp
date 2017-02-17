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
    const Int m = AOrig.Height();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = HermitianOneNorm( uplo, AOrig );

    Output("Testing error...");
    PushIndent();

    // Grab the diagonal and subdiagonal of the symmetric tridiagonal matrix
    Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d = GetRealPartOfDiagonal(A);
    auto e = GetRealPartOfDiagonal(A,subdiagonal);

    // Zero B and then fill its tridiagonal
    Matrix<Field> B;
    Zeros( B, m, m );
    SetRealPartOfDiagonal( B, d );
    SetRealPartOfDiagonal( B, e,  subdiagonal );
    SetRealPartOfDiagonal( B, e, -subdiagonal );
    if( print )
        Print( B, "Tridiagonal" );
    if( display )
        Display( B, "Tridiagonal" );

    // Reverse the accumulated Householder transforms, ignoring symmetry
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, householderScalars, B );
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, householderScalars, B );
    if( print )
        Print( B, "Rotated tridiagonal" );
    if( display )
        Display( B, "Rotated tridiagonal" );

    // Compare the appropriate triangle of AOrig and B
    MakeTrapezoidal( uplo, AOrig );
    MakeTrapezoidal( uplo, B );
    B -= AOrig;
    if( print )
        Print( B, "Error in rotated tridiagonal" );
    if( display )
        Display( B, "Error in rotated tridiagonal" );
    const Real infError = HermitianInfinityNorm( uplo, B );
    const Real relError = infError / (eps*m*oneNormA);
    Output("||A - Q T Q^H||_oo / (eps m ||A||_1) = ",relError);

    // Compute || I - Q Q^H ||
    MakeIdentity( B );
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, householderScalars, B );
    Matrix<Field> QHAdj;
    Adjoint( B, QHAdj );
    MakeIdentity( B );
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, householderScalars, B );
    QHAdj -= B;
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, householderScalars, B );
    ShiftDiagonal( B, Field(-1) );
    const Real infOrthogError = InfinityNorm( B );
    const Real relOrthogError = infOrthogError / (eps*m);
    Output("||I - Q^H Q||_oo / (eps m) = ",relOrthogError);

    PopIndent();

    // TODO(poulson): More rigorous failure conditions
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large");
    if( relOrthogError > Real(10) )
        LogicError("Relative orthogonality error was unacceptably large");
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
    const Int m = AOrig.Height();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = HermitianOneNorm( uplo, AOrig );

    OutputFromRoot(grid.Comm(),"Testing error...");
    PushIndent();

    // Grab the diagonal and subdiagonal of the symmetric tridiagonal matrix
    Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d = GetRealPartOfDiagonal(A);
    auto e = GetRealPartOfDiagonal(A,subdiagonal);

    // Zero B and then fill its tridiagonal
    DistMatrix<Field> B(grid);
    B.AlignWith( A );
    Zeros( B, m, m );
    SetRealPartOfDiagonal( B, d );
    SetRealPartOfDiagonal( B, e,  subdiagonal );
    SetRealPartOfDiagonal( B, e, -subdiagonal );
    if( print )
        Print( B, "Tridiagonal" );
    if( display )
        Display( B, "Tridiagonal" );

    // Reverse the accumulated Householder transforms, ignoring symmetry
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, householderScalars, B );
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, householderScalars, B );
    if( print )
        Print( B, "Rotated tridiagonal" );
    if( display )
        Display( B, "Rotated tridiagonal" );

    // Compare the appropriate triangle of AOrig and B
    MakeTrapezoidal( uplo, AOrig );
    MakeTrapezoidal( uplo, B );
    B -= AOrig;
    if( print )
        Print( B, "Error in rotated tridiagonal" );
    if( display )
        Display( B, "Error in rotated tridiagonal" );
    const Real infError = HermitianInfinityNorm( uplo, B );
    const Real relError = infError / (eps*m*oneNormA);
    OutputFromRoot
    (grid.Comm(),"||A - Q T Q^H||_oo / (eps m ||A||_1) = ",relError);

    // Compute || I - Q Q^H ||
    MakeIdentity( B );
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, householderScalars, B );
    DistMatrix<Field> QHAdj( grid );
    Adjoint( B, QHAdj );
    MakeIdentity( B );
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, householderScalars, B );
    QHAdj -= B;
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, householderScalars, B );
    ShiftDiagonal( B, Field(-1) );
    const Real infOrthogError = InfinityNorm( B );
    const Real relOrthogError = infOrthogError / (eps*m);
    OutputFromRoot(grid.Comm(),"||I - Q^H Q||_oo / (eps m) = ",relOrthogError);

    PopIndent();

    // TODO(poulson): More rigorous failure conditions
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large");
    if( relOrthogError > Real(10) )
        LogicError("Relative orthogonality error was unacceptably large");
}

template<typename Field>
void InnerTestHermitianTridiag
( UpperOrLower uplo,
        Matrix<Field>& A,
        Matrix<Field>& householderScalars,
  bool correctness,
  bool print,
  bool display )
{
    Matrix<Field> AOrig( A ), ACopy( A );
    const Int m = A.Height();
    Timer timer;

    Output("Starting tridiagonalization...");
    timer.Start();
    HermitianTridiag( uplo, A, householderScalars );
    const double runTime = timer.Stop();
    const double realGFlops = 16./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = IsComplex<Field>::value ? 4*realGFlops : realGFlops;
    Output(runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after HermitianTridiag" );
        Print
        ( householderScalars, "householderScalars after HermitianTridiag" );
    }
    if( display )
    {
        Display( A, "A after HermitianTridiag" );
        Display
        ( householderScalars, "householderScalars after HermitianTridiag" );
    }
    if( correctness )
        TestCorrectness( uplo, A, householderScalars, AOrig, print, display );
    A = ACopy;
}

template<typename Field>
void InnerTestHermitianTridiag
( UpperOrLower uplo,
        DistMatrix<Field>& A,
        DistMatrix<Field,STAR,STAR>& householderScalars,
  const HermitianTridiagCtrl<Field>& ctrl,
  bool correctness,
  bool print,
  bool display )
{
    DistMatrix<Field> AOrig( A ), ACopy( A );
    const Int m = A.Height();
    const Grid& grid = A.Grid();
    Timer timer;

    OutputFromRoot(grid.Comm(),"Starting tridiagonalization...");
    mpi::Barrier( grid.Comm() );
    timer.Start();
    HermitianTridiag( uplo, A, householderScalars, ctrl );
    mpi::Barrier( grid.Comm() );
    const double runTime = timer.Stop();
    const double realGFlops = 16./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = IsComplex<Field>::value ? 4*realGFlops : realGFlops;
    OutputFromRoot(grid.Comm(),runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after HermitianTridiag" );
        Print
        ( householderScalars, "householderScalars after HermitianTridiag" );
    }
    if( display )
    {
        Display( A, "A after HermitianTridiag" );
        Display
        ( householderScalars, "householderScalars after HermitianTridiag" );
    }
    if( correctness )
        TestCorrectness( uplo, A, householderScalars, AOrig, print, display );
    A = ACopy;
}

template<typename Field>
void TestHermitianTridiag
( UpperOrLower uplo,
  Int m,
  bool correctness,
  bool print,
  bool display )
{
    Matrix<Field> A, AOrig;
    Matrix<Field> householderScalars;
    Output("Testing with ",TypeName<Field>());
    PushIndent();

    Wigner( A, m );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    if( display )
        Display( A, "A" );

    Output("Sequential algorithm:");
    InnerTestHermitianTridiag
    ( uplo, A, householderScalars, correctness, print, display );

    PopIndent();
}

template<typename Field>
void TestHermitianTridiag
( const Grid& grid,
  UpperOrLower uplo,
  Int m,
  Int nbLocal,
  bool avoidTrmv,
  bool correctness,
  bool print,
  bool display )
{
    DistMatrix<Field> A(grid), AOrig(grid);
    DistMatrix<Field,STAR,STAR> householderScalars(grid);
    OutputFromRoot(grid.Comm(),"Testing with ",TypeName<Field>());
    PushIndent();

    HermitianTridiagCtrl<Field> ctrl;
    ctrl.symvCtrl.bsize = nbLocal;
    ctrl.symvCtrl.avoidTrmvBasedLocalSymv = avoidTrmv;

    Wigner( A, m );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    if( display )
        Display( A, "A" );

    OutputFromRoot(grid.Comm(),"Normal algorithm:");
    ctrl.approach = HERMITIAN_TRIDIAG_NORMAL;
    InnerTestHermitianTridiag
    ( uplo, A, householderScalars, ctrl, correctness, print, display );

    OutputFromRoot(grid.Comm(),"Square row-major algorithm:");
    ctrl.approach = HERMITIAN_TRIDIAG_SQUARE;
    ctrl.order = ROW_MAJOR;
    InnerTestHermitianTridiag
    ( uplo, A, householderScalars, ctrl, correctness, print, display );

    OutputFromRoot(grid.Comm(),"Square column-major algorithm:");
    ctrl.approach = HERMITIAN_TRIDIAG_SQUARE;
    ctrl.order = COLUMN_MAJOR;
    InnerTestHermitianTridiag
    ( uplo, A, householderScalars, ctrl, correctness, print, display );
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
        const Int m = Input("--height","height of matrix",75);
        const Int nb = Input("--nb","algorithmic blocksize",32);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool avoidTrmv =
          Input("--avoidTrmv","avoid Trmv local Symv",true);
        const bool sequential = Input("--sequential","test sequential?",true);
        const bool correctness =
          Input("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool display = Input("--display","display matrices?",false);
        const bool testReal = Input("--testReal","test real matrices?",true);
        const bool testCpx = Input("--testCpx","test complex matrices?",true);
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
        OutputFromRoot(grid.Comm(),"Will test HermitianTridiag",uploChar);

        if( sequential && mpi::Rank() == 0 )
        {
            if( testReal )
                TestHermitianTridiag<float>
                ( uplo, m, correctness, print, display );
            if( testCpx )
                TestHermitianTridiag<Complex<float>>
                ( uplo, m, correctness, print, display );

            if( testReal )
                TestHermitianTridiag<double>
                ( uplo, m, correctness, print, display );
            if( testCpx )
                TestHermitianTridiag<Complex<double>>
                ( uplo, m, correctness, print, display );

#ifdef EL_HAVE_QD
            if( testReal )
            {
                TestHermitianTridiag<DoubleDouble>
                ( uplo, m, correctness, print, display );
                TestHermitianTridiag<QuadDouble>
                ( uplo, m, correctness, print, display );
            }
            if( testCpx )
            {
                TestHermitianTridiag<Complex<DoubleDouble>>
                ( uplo, m, correctness, print, display );
                TestHermitianTridiag<Complex<QuadDouble>>
                ( uplo, m, correctness, print, display );
            }
#endif

#ifdef EL_HAVE_QUAD
            if( testReal )
                TestHermitianTridiag<Quad>
                ( uplo, m, correctness, print, display );
            if( testCpx )
                TestHermitianTridiag<Complex<Quad>>
                ( uplo, m, correctness, print, display );
#endif

#ifdef EL_HAVE_MPC
            if( testReal )
                TestHermitianTridiag<BigFloat>
                ( uplo, m, correctness, print, display );
#endif
        }

        if( testReal )
            TestHermitianTridiag<float>
            ( grid, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
        if( testCpx )
            TestHermitianTridiag<Complex<float>>
            ( grid, uplo, m, nbLocal, avoidTrmv, correctness, print, display );

        if( testReal )
            TestHermitianTridiag<double>
            ( grid, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
        if( testCpx )
            TestHermitianTridiag<Complex<double>>
            ( grid, uplo, m, nbLocal, avoidTrmv, correctness, print, display );

#ifdef EL_HAVE_QD
        if( testReal )
        {
            TestHermitianTridiag<DoubleDouble>
            ( grid, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
            TestHermitianTridiag<QuadDouble>
            ( grid, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
        }
        if( testCpx )
        {
            TestHermitianTridiag<Complex<DoubleDouble>>
            ( grid, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
            TestHermitianTridiag<Complex<QuadDouble>>
            ( grid, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
        }
#endif

#ifdef EL_HAVE_QUAD
        if( testReal )
            TestHermitianTridiag<Quad>
            ( grid, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
        if( testCpx )
            TestHermitianTridiag<Complex<Quad>>
            ( grid, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
#endif

#ifdef EL_HAVE_MPC
        if( testReal )
            TestHermitianTridiag<BigFloat>
            ( grid, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
