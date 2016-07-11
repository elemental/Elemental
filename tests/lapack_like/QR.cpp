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
( const Matrix<F>& A,
  const Matrix<F>& phase,
  const Matrix<Base<F>>& d,
        Matrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int maxDim = Max(m,n);
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( AOrig );

    Output("Testing orthogonality of Q");
    PushIndent();

    // Form Z := Q^H Q as an approximation to identity
    Matrix<F> Z;
    Identity( Z, m, n );
    qr::ApplyQ( LEFT, NORMAL, A, phase, d, Z );
    qr::ApplyQ( LEFT, ADJOINT, A, phase, d, Z );
    auto ZUpper = Z( IR(0,minDim), IR(0,minDim) );

    // Form X := I - Q^H Q
    Matrix<F> X;
    Identity( X, minDim, minDim );
    X -= ZUpper;

    const Real infOrthogError = InfinityNorm( X );
    const Real relOrthogError = infOrthogError / (eps*maxDim);
    Output("||Q^H Q - I||_oo / (eps Max(m,n)) = ",relOrthogError);
    PopIndent();

    Output("Testing if A ~= QR");
    PushIndent();

    // Form Q R
    auto U( A );
    MakeTrapezoidal( UPPER, U );
    qr::ApplyQ( LEFT, NORMAL, A, phase, d, U );
    U -= AOrig;
    const Real infError = InfinityNorm( U ); 
    const Real relError = infError / (eps*maxDim*oneNormA);
    Output("||A - Q R||_oo / (eps Max(m,n) ||A||_1) = ",relError);

    PopIndent();

    // TODO: More rigorous failure condition
    if( relOrthogError > Real(10) )
        LogicError("Relative orthogonality error was unacceptably large");
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large"); 
}

template<typename F>
void TestCorrectness
( const DistMatrix<F>& A,
  const DistMatrix<F,MD,STAR>& phase,
  const DistMatrix<Base<F>,MD,STAR>& signature,
        DistMatrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int maxDim = Max(m,n);
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( AOrig );

    OutputFromRoot(g.Comm(),"Testing orthogonality of Q");
    PushIndent();

    // Form Z := Q^H Q as an approximation to identity
    DistMatrix<F> Z(g);
    Identity( Z, m, n );
    qr::ApplyQ( LEFT, NORMAL, A, phase, signature, Z );
    qr::ApplyQ( LEFT, ADJOINT, A, phase, signature, Z );
    auto ZUpper = Z( IR(0,minDim), IR(0,minDim) );

    // Form X := I - Q^H Q
    DistMatrix<F> X(g);
    Identity( X, minDim, minDim );
    X -= ZUpper;

    const Real infOrthogError = InfinityNorm( X );
    const Real relOrthogError = infOrthogError / (eps*maxDim);
    OutputFromRoot
    (g.Comm(),"||Q^H Q - I||_oo / (eps Max(m,n)) = ",relOrthogError);
    PopIndent();

    OutputFromRoot(g.Comm(),"Testing if A ~= QR");
    PushIndent();

    // Form Q R
    auto U( A );
    MakeTrapezoidal( UPPER, U );
    qr::ApplyQ( LEFT, NORMAL, A, phase, signature, U );
    U -= AOrig;
    const Real infError = InfinityNorm( U ); 
    const Real relError = infError / (eps*maxDim*oneNormA);
    OutputFromRoot
    (g.Comm(),"||A - Q R||_oo / (eps Max(m,n) ||A||_1) = ",relError);

    PopIndent();

    // TODO: More rigorous failure condition
    if( relOrthogError > Real(10) )
        LogicError("Relative orthogonality error was unacceptably large");
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large"); 
}

template<typename F>
void TestQR
( Int m,
  Int n,
  bool correctness,
  bool print )
{
    Output("Testing with ",TypeName<F>());
    PushIndent();
    Matrix<F> A, AOrig;
    Matrix<F> phase;
    Matrix<Base<F>> signature;

    Uniform( A, m, n );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    const double mD = double(m);
    const double nD = double(n);

    Timer timer;
    Output("Starting QR factorization...");
    timer.Start();
    QR( A, phase, signature );
    const double runTime = timer.Stop();
    const double realGFlops = (2.*mD*nD*nD - 2./3.*nD*nD*nD)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    Output("Elemental: ",runTime," seconds. GFlops = ",gFlops);
    if( print )
    {
        Print( A, "A after factorization" );
        Print( phase, "phase" );
        Print( signature, "signature" );
    }
    if( correctness )
        TestCorrectness( A, phase, signature, AOrig );
    PopIndent();
}

template<typename F,typename=EnableIf<IsBlasScalar<F>>>
void TestQR
( const Grid& g,
  Int m,
  Int n,
  bool correctness,
  bool print,
  bool scalapack )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<F,MD,STAR> phase(g);
    DistMatrix<Base<F>,MD,STAR> signature(g);

    Uniform( A, m, n );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    const double mD = double(m);
    const double nD = double(n);

    Timer timer;
    if( scalapack )
    {
        // TODO(poulson): Fold this interface into the standard QR via a Ctrl
        // option of 'useScaLAPACK'
        DistMatrix<F,MC,MR,BLOCK> ABlock( A );
        DistMatrix<F,MR,STAR,BLOCK> phaseBlock(g);
        mpi::Barrier( g.Comm() );
        timer.Start();
        QR( ABlock, phaseBlock ); 
        const double runTime = timer.Stop();
        const double realGFlops = (2.*mD*nD*nD - 2./3.*nD*nD*nD)/(1.e9*runTime);
        const double gFlops =
          ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
        OutputFromRoot
        (g.Comm(),"ScaLAPACK: ",runTime," seconds. GFlops = ",gFlops);
    }

    OutputFromRoot(g.Comm(),"Starting QR factorization...");
    mpi::Barrier( g.Comm() );
    timer.Start();
    QR( A, phase, signature );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    const double realGFlops = (2.*mD*nD*nD - 2./3.*nD*nD*nD)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    OutputFromRoot(g.Comm(),"Elemental: ",runTime," seconds. GFlops = ",gFlops);
    if( print )
    {
        Print( A, "A after factorization" );
        Print( phase, "phase" );
        Print( signature, "signature" );
    }
    if( correctness )
        TestCorrectness( A, phase, signature, AOrig );
    PopIndent();
}

template<typename F,typename=DisableIf<IsBlasScalar<F>>,typename=void>
void TestQR
( const Grid& g,
  Int m,
  Int n,
  bool correctness,
  bool print )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<F,MD,STAR> phase(g);
    DistMatrix<Base<F>,MD,STAR> signature(g);

    Uniform( A, m, n );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    const double mD = double(m);
    const double nD = double(n);

    OutputFromRoot(g.Comm(),"Starting QR factorization...");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    QR( A, phase, signature );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = (2.*mD*nD*nD - 2./3.*nD*nD*nD)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    OutputFromRoot(g.Comm(),"Elemental: ",runTime," seconds. GFlops = ",gFlops);
    if( print )
    {
        Print( A, "A after factorization" );
        Print( phase, "phase" );
        Print( signature, "signature" );
    }
    if( correctness )
        TestCorrectness( A, phase, signature, AOrig );
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
        const Int nb = Input("--nb","algorithmic blocksize",64);
        const bool sequential = Input("--sequential","test sequential?",true);
        const bool correctness =
          Input("--correctness","test correctness?",true);
#ifdef EL_HAVE_SCALAPACK
        const bool scalapack = Input("--scalapack","test ScaLAPACK?",true);
#else
        const bool scalapack = false;
#endif
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        const bool print = Input("--print","print matrices?",false);
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
            TestQR<float>
            ( m, n, correctness, print );
            TestQR<Complex<float>>
            ( m, n, correctness, print );

            TestQR<double>
            ( m, n, correctness, print );
            TestQR<Complex<double>>
            ( m, n, correctness, print );

#ifdef EL_HAVE_QD
            TestQR<DoubleDouble>
            ( m, n, correctness, print );
            TestQR<QuadDouble>
            ( m, n, correctness, print );

            TestQR<Complex<DoubleDouble>>
            ( m, n, correctness, print );
            TestQR<Complex<QuadDouble>>
            ( m, n, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
            TestQR<Quad>
            ( m, n, correctness, print );
            TestQR<Complex<Quad>>
            ( m, n, correctness, print );
#endif

#ifdef EL_HAVE_MPC
            TestQR<BigFloat>
            ( m, n, correctness, print );
            TestQR<Complex<BigFloat>>
            ( m, n, correctness, print );
#endif
        }

        TestQR<float>
        ( g, m, n, correctness, print, scalapack );
        TestQR<Complex<float>>
        ( g, m, n, correctness, print, scalapack );

        TestQR<double>
        ( g, m, n, correctness, print, scalapack );
        TestQR<Complex<double>>
        ( g, m, n, correctness, print, scalapack );

#ifdef EL_HAVE_QD
        TestQR<DoubleDouble>
        ( g, m, n, correctness, print );
        TestQR<QuadDouble>
        ( g, m, n, correctness, print );

        TestQR<Complex<DoubleDouble>>
        ( g, m, n, correctness, print );
        TestQR<Complex<QuadDouble>>
        ( g, m, n, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
        TestQR<Quad>
        ( g, m, n, correctness, print );
        TestQR<Complex<Quad>>
        ( g, m, n, correctness, print );
#endif

#ifdef EL_HAVE_MPC
        TestQR<BigFloat>
        ( g, m, n, correctness, print );
        TestQR<Complex<BigFloat>>
        ( g, m, n, correctness, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
