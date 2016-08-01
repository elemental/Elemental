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
( UpperOrLower uplo, 
  const Matrix<F>& A, 
  const Matrix<F>& phase,
        Matrix<F>& AOrig,
  bool print,
  bool display )
{
    typedef Base<F> Real;
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
    Matrix<F> B;
    Zeros( B, m, m );
    SetRealPartOfDiagonal( B, d );
    SetRealPartOfDiagonal( B, e,  subdiagonal );
    SetRealPartOfDiagonal( B, e, -subdiagonal );
    if( print )
        Print( B, "Tridiagonal" );
    if( display )
        Display( B, "Tridiagonal" );

    // Reverse the accumulated Householder transforms, ignoring symmetry
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, phase, B );
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, phase, B );
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
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, phase, B );
    Matrix<F> QHAdj;
    Adjoint( B, QHAdj );
    MakeIdentity( B );
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, phase, B );
    QHAdj -= B;
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, phase, B );
    ShiftDiagonal( B, F(-1) );
    const Real infOrthogError = InfinityNorm( B );
    const Real relOrthogError = infOrthogError / (eps*m);
    Output("||I - Q^H Q||_oo / (eps m) = ",relOrthogError);

    PopIndent();

    // TODO: More rigorous failure conditions
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large");
    if( relOrthogError > Real(10) )
        LogicError("Relative orthogonality error was unacceptably large");
}

template<typename F> 
void TestCorrectness
( UpperOrLower uplo, 
  const DistMatrix<F>& A, 
  const DistMatrix<F,STAR,STAR>& phase,
        DistMatrix<F>& AOrig,
  bool print,
  bool display )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = AOrig.Height();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = HermitianOneNorm( uplo, AOrig );

    OutputFromRoot(g.Comm(),"Testing error...");
    PushIndent();

    // Grab the diagonal and subdiagonal of the symmetric tridiagonal matrix
    Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d = GetRealPartOfDiagonal(A);
    auto e = GetRealPartOfDiagonal(A,subdiagonal);
     
    // Zero B and then fill its tridiagonal
    DistMatrix<F> B(g);
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
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, phase, B );
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, phase, B );
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
    (g.Comm(),"||A - Q T Q^H||_oo / (eps m ||A||_1) = ",relError);

    // Compute || I - Q Q^H ||
    MakeIdentity( B );
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, phase, B );
    DistMatrix<F> QHAdj( g );
    Adjoint( B, QHAdj );
    MakeIdentity( B );
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, phase, B );
    QHAdj -= B;
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, phase, B );
    ShiftDiagonal( B, F(-1) );
    const Real infOrthogError = InfinityNorm( B );
    const Real relOrthogError = infOrthogError / (eps*m);
    OutputFromRoot(g.Comm(),"||I - Q^H Q||_oo / (eps m) = ",relOrthogError);

    PopIndent();

    // TODO: More rigorous failure conditions
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large");
    if( relOrthogError > Real(10) )
        LogicError("Relative orthogonality error was unacceptably large");
}

template<typename F>
void InnerTestHermitianTridiag
( UpperOrLower uplo,
        Matrix<F>& A,
        Matrix<F>& phase,
  bool correctness,
  bool print,
  bool display )
{
    Matrix<F> AOrig( A ), ACopy( A );
    const Int m = A.Height();
    Timer timer;

    Output("Starting tridiagonalization...");
    timer.Start();
    HermitianTridiag( uplo, A, phase );
    const double runTime = timer.Stop();
    const double realGFlops = 16./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    Output(runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after HermitianTridiag" );
        Print( phase, "phase after HermitianTridiag" );
    }
    if( display )
    {
        Display( A, "A after HermitianTridiag" );
        Display( phase, "phase after HermitianTridiag" );
    }
    if( correctness )
        TestCorrectness( uplo, A, phase, AOrig, print, display );
    A = ACopy;
}

template<typename F>
void InnerTestHermitianTridiag
( UpperOrLower uplo,
        DistMatrix<F>& A,
        DistMatrix<F,STAR,STAR>& phase,
  const HermitianTridiagCtrl<F>& ctrl,
  bool correctness,
  bool print,
  bool display )
{
    DistMatrix<F> AOrig( A ), ACopy( A );
    const Int m = A.Height();
    const Grid& g = A.Grid();
    Timer timer;

    OutputFromRoot(g.Comm(),"Starting tridiagonalization...");
    mpi::Barrier( g.Comm() );
    timer.Start();
    HermitianTridiag( uplo, A, phase, ctrl );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    const double realGFlops = 16./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    OutputFromRoot(g.Comm(),runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after HermitianTridiag" );
        Print( phase, "phase after HermitianTridiag" );
    }
    if( display )
    {
        Display( A, "A after HermitianTridiag" );
        Display( phase, "phase after HermitianTridiag" );
    }
    if( correctness )
        TestCorrectness( uplo, A, phase, AOrig, print, display );
    A = ACopy;
}

template<typename F>
void TestHermitianTridiag
( UpperOrLower uplo,
  Int m,
  bool correctness,
  bool print,
  bool display )
{
    Matrix<F> A, AOrig;
    Matrix<F> phase;
    Output("Testing with ",TypeName<F>());
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
    ( uplo, A, phase, correctness, print, display );

    PopIndent();
}

template<typename F>
void TestHermitianTridiag
( const Grid& g,
  UpperOrLower uplo,
  Int m,
  Int nbLocal,
  bool avoidTrmv,
  bool correctness,
  bool print,
  bool display )
{
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<F,STAR,STAR> phase(g);
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();

    HermitianTridiagCtrl<F> ctrl;
    ctrl.symvCtrl.bsize = nbLocal;
    ctrl.symvCtrl.avoidTrmvBasedLocalSymv = avoidTrmv;

    Wigner( A, m );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    if( display )
        Display( A, "A" );

    OutputFromRoot(g.Comm(),"Normal algorithm:");
    ctrl.approach = HERMITIAN_TRIDIAG_NORMAL;
    InnerTestHermitianTridiag
    ( uplo, A, phase, ctrl, correctness, print, display );

    OutputFromRoot(g.Comm(),"Square row-major algorithm:");
    ctrl.approach = HERMITIAN_TRIDIAG_SQUARE;
    ctrl.order = ROW_MAJOR;
    InnerTestHermitianTridiag
    ( uplo, A, phase, ctrl, correctness, print, display );

    OutputFromRoot(g.Comm(),"Square column-major algorithm:");
    ctrl.approach = HERMITIAN_TRIDIAG_SQUARE;
    ctrl.order = COLUMN_MAJOR;
    InnerTestHermitianTridiag
    ( uplo, A, phase, ctrl, correctness, print, display );
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
            gridHeight = Grid::FindFactor( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        OutputFromRoot(g.Comm(),"Will test HermitianTridiag",uploChar);

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
            ( g, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
        if( testCpx )
            TestHermitianTridiag<Complex<float>>
            ( g, uplo, m, nbLocal, avoidTrmv, correctness, print, display );

        if( testReal )
            TestHermitianTridiag<double>
            ( g, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
        if( testCpx )
            TestHermitianTridiag<Complex<double>>
            ( g, uplo, m, nbLocal, avoidTrmv, correctness, print, display );

#ifdef EL_HAVE_QD
        if( testReal )
        {
            TestHermitianTridiag<DoubleDouble>
            ( g, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
            TestHermitianTridiag<QuadDouble>
            ( g, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
        }
        if( testCpx )
        {
            TestHermitianTridiag<Complex<DoubleDouble>>
            ( g, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
            TestHermitianTridiag<Complex<QuadDouble>>
            ( g, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
        }
#endif

#ifdef EL_HAVE_QUAD
        if( testReal )
            TestHermitianTridiag<Quad>
            ( g, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
        if( testCpx )
            TestHermitianTridiag<Complex<Quad>>
            ( g, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
#endif

#ifdef EL_HAVE_MPC
        if( testReal )
            TestHermitianTridiag<BigFloat>
            ( g, uplo, m, nbLocal, avoidTrmv, correctness, print, display );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
