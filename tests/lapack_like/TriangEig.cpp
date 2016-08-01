/*
   Copyright (c) 2009-2016, Jack Poulson and Tim Moon
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename F>
void GrcarSchurFactor( Matrix<F>& A, Int m )
{
    LogicError("Only complex Grcar is allowed for these tests");
}

template<typename Real>
void GrcarSchurFactor( Matrix<Complex<Real>>& A, Int m )
{
    Matrix<Complex<Real>> d;
    Grcar( A, m );
    Schur( A, d );
    MakeTrapezoidal( UPPER, A, -1 );
}

template<typename F>
void GrcarSchurFactor( ElementalMatrix<F>& A, Int m )
{
    LogicError("Only complex Grcar is allowed for these tests");
}

template<typename Real>
void GrcarSchurFactor( ElementalMatrix<Complex<Real>>& A, Int m )
{
    DistMatrix<Complex<Real>> d(A.Grid());
    Grcar( A, m );
    Schur( A, d );
    MakeTrapezoidal( UPPER, A, -1 );
}

template<typename F>
void FoxLiSchurFactor( Matrix<F>& A, Int m )
{
    LogicError("Fox-Li matrix is complex");
}

template<typename Real>
void FoxLiSchurFactor( Matrix<Complex<Real>>& A, Int m )
{
    Matrix<Complex<Real>> B;
    Matrix<Complex<Real>> w;
    FoxLi( B, m, Real(-0.179) );
    Schur( B, w );
    MakeTrapezoidal( UPPER, B, 0 );
    A = B;
}

template<typename F>
void FoxLiSchurFactor( ElementalMatrix<F>& A, Int m )
{
    LogicError("Fox-Li matrix is complex");
}

template<typename Real>
void FoxLiSchurFactor( ElementalMatrix<Complex<Real>>& A, Int m )
{
    const Grid& g = A.Grid();
    DistMatrix<Complex<Real>> B(g);
    DistMatrix<Complex<Real>> w(g);
    FoxLi( B, m, Real(-0.179) );
    Schur( B, w );
    MakeTrapezoidal( UPPER, B, 0 );
    A = B;
}

template<typename F>
void TestCorrectness
( const Matrix<F>& A,
  const Matrix<F>& X,
        bool print )
{
    typedef Base<F> Real;
    const Int n = A.Height();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( A );

    // Find the residual R = AX-XW
    Matrix<F> R( X );
    Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, R );
    Matrix<F> XW( X );
    Matrix<F> w;
    GetDiagonal( A, w );
    DiagonalScale( RIGHT, NORMAL, w, XW );
    R -= XW;
    const Real infError = InfinityNorm( R );
    const Real relError = infError / (eps*n*oneNormA);
    Output("||A X - X W||_oo / (eps n ||A||_1) = ",relError);

    // TODO: More rigorous failure condition
    if( relError > Real(10) || !limits::IsFinite(relError) )
        LogicError("Unacceptably large relative error");
}

template<typename F>
void TestCorrectness
( const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& X,
        bool print )
{
    typedef Base<F> Real;
    const Int n = A.Height();
    const Grid& g = A.Grid();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( A );

    // Find the residual R = AX-XW
    DistMatrix<F> R( X );
    Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, R );
    DistMatrix<F> XW( X );
    DistMatrix<F> w( g );
    GetDiagonal( A, w );
    DiagonalScale( RIGHT, NORMAL, w, XW );
    R -= XW;
    const Real infError = InfinityNorm( R );
    const Real relError = infError / (eps*n*oneNormA);
    OutputFromRoot
    (g.Comm(),"||A X - X W||_oo / (eps n ||A||_1) = ",relError);
    // TODO: More rigorous failure condition
    if( relError > Real(10) || !limits::IsFinite(relError) )
        LogicError("Unacceptably large relative error");
}

template<typename F,typename=EnableIf<IsBlasScalar<F>>>
void TestTriangEig
( Int m,
  bool correctness,
  bool print,
  Int whichMatrix )
{
    Matrix<F> A, AOrig, X;
    Matrix<F> w;

    Output("Testing with ",TypeName<F>());
    PushIndent();

    // Generate test matrix
    switch( whichMatrix )
    {
    case 0:
        {
            // LU factorization of Gaussian matrix
            Matrix<F> B;
            Gaussian( B, m, m );
            LU( B );
            Transpose( B, A );
            MakeTrapezoidal( UPPER, A, 0 );
            break;
        }
    case 1:
        {
            // Schur factorization of Fox-Li matrix
            FoxLiSchurFactor( A, m );
            break;
        }
    case 2:
        {
            // Schur factorization of (complex) Grcar matrix
            GrcarSchurFactor( A, m );
            break;
        }
    case 3:
        Jordan( A, m, F(7) );
        break;
    default: LogicError("Unknown test matrix");
    }
    
    if( correctness )
    {
        AOrig = A;
        GetDiagonal( A, w );
    }
    if( print )
        Print( A, "A" );

    Output("Starting triangular eigensolver...");
    Timer timer;
    timer.Start();
    TriangEig( A, X );
    Output("Time = ",timer.Stop()," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        Print( X, "eigenvectors:" );
    }
    if( correctness )
        TestCorrectness( AOrig, X, print );
    PopIndent();
}

template<typename F,typename=DisableIf<IsBlasScalar<F>>,typename=void>
void TestTriangEig
( Int m,
  bool correctness,
  bool print,
  Int whichMatrix )
{
    Matrix<F> A, AOrig, X;
    Matrix<F> w;

    Output("Testing with ",TypeName<F>());
    PushIndent();

    // Generate test matrix
    switch( whichMatrix )
    {
    case 0:
        {
            // LU factorization of Gaussian matrix
            Matrix<F> B;
            Gaussian( B, m, m );
            LU( B );
            Transpose( B, A );
            MakeTrapezoidal( UPPER, A, 0 );
            break;
        }
    case 3:
        Jordan( A, m, F(7) );
        break;
    default: LogicError("Schur factorization not supported for non-BLAS types");
    }
    
    if( correctness )
    {
        AOrig = A;
        GetDiagonal( A, w );
    }
    if( print )
        Print( A, "A" );

    Output("Starting triangular eigensolver...");
    Timer timer;
    timer.Start();
    TriangEig( A, X );
    Output("Time = ",timer.Stop()," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        Print( X, "eigenvectors:" );
    }
    if( correctness )
        TestCorrectness( AOrig, X, print );
    PopIndent();
}

template<typename F,Dist U=MC,Dist V=MR,Dist S=MC,
         typename=EnableIf<IsBlasScalar<F>>>
void TestTriangEig
( const Grid& g,
  Int m,
  bool correctness,
  bool print,
  Int whichMatrix )
{
    DistMatrix<F,U,V> A(g), AOrig(g), X(g);
    DistMatrix<F,S,STAR> w(g);

    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();

    // Generate test matrix
    switch( whichMatrix )
    {
    case 0:
        {
            // LU factorization of Gaussian matrix
            DistMatrix<F> B(g);
            Gaussian( B, m, m );
            LU( B );
            Transpose( B, A );
            MakeTrapezoidal( UPPER, A, 0 );
            break;
        }
    case 1:
        {
            // Schur factorization of Fox-Li matrix
            FoxLiSchurFactor( A, m );
            break;
        }
    case 2:
        {
            // Schur factorization of (complex) Grcar matrix
            GrcarSchurFactor( A, m );
            break;
        }
    case 3:
        Jordan( A, m, F(7) );
        break;
    default: LogicError("Unknown test matrix");
    }
    
    if( correctness )
    {
        AOrig = A;
        GetDiagonal( A, w );
    }
    if( print )
        Print( A, "A" );

    OutputFromRoot(g.Comm(),"Starting triangular eigensolver...");
    Timer timer;
    timer.Start();
    TriangEig( A, X );
    OutputFromRoot(g.Comm(),"Time = ",timer.Stop()," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        Print( X, "eigenvectors:" );
    }
    if( correctness )
        TestCorrectness( AOrig, X, print );
    PopIndent();
}

template<typename F,Dist U=MC,Dist V=MR,Dist S=MC,
         typename=DisableIf<IsBlasScalar<F>>,typename=void>
void TestTriangEig
( const Grid& g,
  Int m,
  bool correctness,
  bool print,
  Int whichMatrix )
{
    DistMatrix<F,U,V> A(g), AOrig(g), X(g);
    DistMatrix<F,S,STAR> w(g);

    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();

    // Generate test matrix
    switch( whichMatrix )
    {
    case 0:
        {
            // LU factorization of Gaussian matrix
            DistMatrix<F> B(g);
            Gaussian( B, m, m );
            LU( B );
            Transpose( B, A );
            MakeTrapezoidal( UPPER, A, 0 );
            break;
        }
    case 3:
        Jordan( A, m, F(7) );
        break;
    default: LogicError("Schur factorization not supported for non-BLAS types");
    }
    
    if( correctness )
    {
        AOrig = A;
        GetDiagonal( A, w );
    }
    if( print )
        Print( A, "A" );

    OutputFromRoot(g.Comm(),"Starting triangular eigensolver...");
    Timer timer;
    timer.Start();
    TriangEig( A, X );
    OutputFromRoot(g.Comm(),"Time = ",timer.Stop()," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        Print( X, "eigenvectors:" );
    }
    if( correctness )
        TestCorrectness( AOrig, X, print );
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
        const Int n = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool sequential = Input("--sequential","test sequential?",true);
        const bool correctness =
          Input("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        bool testReal = Input("--testReal","test real matrices?",true);
        const bool testCpx = Input("--testCpx","test complex matrices?",true);
        const Int whichMatrix =
          Input
          ("--whichMatrix","(0=Gaussian,1=Fox-Li,2=Grcar,3=Jordan)",0);
        ProcessInput();
        PrintInputReport();

        bool testNonstandard = true;
        if( whichMatrix == 1 || whichMatrix == 2 )
        {
            testReal = false;
            testNonstandard = false;
        }

        if( gridHeight == 0 )
            gridHeight = Grid::FindFactor( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );
        SetBlocksize( nb );
        ComplainIfDebug();

        // Test with default distributions
        OutputFromRoot(g.Comm(),"Normal algorithms:");

        if( sequential && mpi::Rank() == 0 )
        {
            // NOTE: This is not nearly enough precision for our chosen matrices
            /*
            if( testReal )
                TestTriangEig<float>
                ( n, correctness, print, whichMatrix );
            if( testCpx )
                TestTriangEig<Complex<float>>
                ( n, correctness, print, whichMatrix );
            */

            if( testReal )
                TestTriangEig<double>
                ( n, correctness, print, whichMatrix );
            if( testCpx )
                TestTriangEig<Complex<double>>
                ( n, correctness, print, whichMatrix );

#ifdef EL_HAVE_QUAD
            if( testReal && testNonstandard )
                TestTriangEig<Quad>
                ( n, correctness, print, whichMatrix );
            if( testCpx && testNonstandard )
                TestTriangEig<Complex<Quad>>
                ( n, correctness, print, whichMatrix );
#endif

#ifdef EL_HAVE_QD
            if( testReal && testNonstandard )
            {
                TestTriangEig<DoubleDouble>
                ( n, correctness, print, whichMatrix );
                TestTriangEig<QuadDouble>
                ( n, correctness, print, whichMatrix );
            }
            if( testCpx && testNonstandard )
            {
                TestTriangEig<Complex<DoubleDouble>>
                ( n, correctness, print, whichMatrix );
                TestTriangEig<Complex<QuadDouble>>
                ( n, correctness, print, whichMatrix );
            }
#endif

#ifdef EL_HAVE_MPC
            if( testReal && testNonstandard )
                TestTriangEig<BigFloat>
                ( n, correctness, print, whichMatrix );
            if( testCpx && testNonstandard )
                TestTriangEig<Complex<BigFloat>>
                ( n, correctness, print, whichMatrix );
#endif
        }

        // NOTE: This is not nearly enough precision for our chosen matrices
        /*
        if( testReal )
            TestTriangEig<float>
            ( g, n, correctness, print, whichMatrix );
        if( testCpx )
            TestTriangEig<Complex<float>>
            ( g, n, correctness, print, whichMatrix );
        */

        if( testReal )
            TestTriangEig<double>
            ( g, n, correctness, print, whichMatrix );
        if( testCpx )
            TestTriangEig<Complex<double>>
            ( g, n, correctness, print, whichMatrix );

#ifdef EL_HAVE_QUAD
        if( testReal && testNonstandard )
            TestTriangEig<Quad>
            ( g, n, correctness, print, whichMatrix );
        if( testCpx && testNonstandard )
            TestTriangEig<Complex<Quad>>
            ( g, n, correctness, print, whichMatrix );
#endif

#ifdef EL_HAVE_QD
        if( testReal && testNonstandard )
        {
            TestTriangEig<DoubleDouble>
            ( g, n, correctness, print, whichMatrix );
            TestTriangEig<QuadDouble>
            ( g, n, correctness, print, whichMatrix );
        }
        if( testCpx && testNonstandard )
        {
            TestTriangEig<Complex<DoubleDouble>>
            ( g, n, correctness, print, whichMatrix );
            TestTriangEig<Complex<QuadDouble>>
            ( g, n, correctness, print, whichMatrix );
        }
#endif

#ifdef EL_HAVE_MPC
        if( testReal && testNonstandard )
            TestTriangEig<BigFloat>
            ( g, n, correctness, print, whichMatrix );
        if( testCpx && testNonstandard )
            TestTriangEig<Complex<BigFloat>>
            ( g, n, correctness, print, whichMatrix );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
