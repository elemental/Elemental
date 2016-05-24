/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace std;
using namespace El;

template<typename F> 
void TestCorrectness
( const DistMatrix<F,VC,  STAR>& Q,
  const DistMatrix<F,STAR,STAR>& R,
        DistMatrix<F,VC,  STAR>& A )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( A );

    // Form I - Q^H Q
    OutputFromRoot(g.Comm(),"Testing orthogonality of Q");
    PushIndent();
    DistMatrix<F> Z(g);
    Identity( Z, n, n );
    DistMatrix<F> Q_MC_MR( Q );
    Herk( UPPER, ADJOINT, Base<F>(-1), Q_MC_MR, Base<F>(1), Z );
    const Real infOrthogError = HermitianInfinityNorm( UPPER, Z );
    const Real relOrthogError = infOrthogError / (eps*Max(m,n)*oneNormA);
    OutputFromRoot
    (g.Comm(),"||Q^H Q - I||_oo / (eps Max(m,n) || A ||_1) = ",relOrthogError);
    PopIndent();

    // Form A - Q R
    OutputFromRoot(g.Comm(),"Testing if A = QR");
    PushIndent();
    LocalGemm( NORMAL, NORMAL, F(-1), Q, R, F(1), A );
    const Real infNormError = InfinityNorm( A );
    const Real relError = infNormError / (eps*Max(m,n)*oneNormA);
    OutputFromRoot
    (g.Comm(),"||A - QR||_oo / (eps Max(m,n) || A ||_1) = ",relError);
    PopIndent();

    // TODO: More refined failure conditions (especially in this case...)
    if( relOrthogError > Real(10) )
        LogicError("Relative orthog error was unacceptably large");
    if( relError > Real(1) )
        LogicError("Relative error was unacceptably large");
}

template<typename F>
void TestQR
( const Grid& g,
  Int m, 
  Int n,
  bool testCorrectness,
  bool print )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();
    DistMatrix<F,VC,STAR> A(g), Q(g);
    DistMatrix<F,STAR,STAR> R(g);

    Uniform( A, m, n );
    if( print )
        Print( A, "A" );
    Q = A;

    OutputFromRoot(g.Comm(),"Starting Cholesky QR factorization");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    qr::Cholesky( Q, R );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    const double mD = double(m);
    const double nD = double(n);
    const double gFlops = (2.*mD*nD*nD + 1./3.*nD*nD*nD)/(1.e9*runTime);
    OutputFromRoot(g.Comm(),"Time: ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( Q, "Q" );
        Print( R, "R" );
    }
    if( testCorrectness )
        TestCorrectness( Q, R, A );
    PopIndent();
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        ProcessInput();
        PrintInputReport();

#ifdef EL_HAVE_MPC
        mpfr::SetPrecision( prec );
#endif

        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, order );
        SetBlocksize( nb );
        ComplainIfDebug();

        TestQR<float>( g, m, n, testCorrectness, print );
        TestQR<Complex<float>>( g, m, n, testCorrectness, print );

        TestQR<double>( g, m, n, testCorrectness, print );
        TestQR<Complex<double>>( g, m, n, testCorrectness, print );

#ifdef EL_HAVE_QD
        TestQR<DoubleDouble>( g, m, n, testCorrectness, print );
        TestQR<QuadDouble>( g, m, n, testCorrectness, print );
#endif

#ifdef EL_HAVE_QUAD
        TestQR<Quad>( g, m, n, testCorrectness, print );
        TestQR<Complex<Quad>>( g, m, n, testCorrectness, print );
#endif

#ifdef EL_HAVE_MPC
        TestQR<BigFloat>( g, m, n, testCorrectness, print );
        TestQR<Complex<BigFloat>>( g, m, n, testCorrectness, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
