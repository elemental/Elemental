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
( const DistMatrix<F,VC,  STAR>& Q,
  const DistMatrix<F,STAR,STAR>& R,
        DistMatrix<F,VC,  STAR>& A )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Int maxDim = Max(m,n);
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( A );

    // Form I - Q^H Q
    OutputFromRoot(g.Comm(),"Testing orthogonality of Q...");
    PushIndent();
    DistMatrix<F> Z(g);
    Identity( Z, n, n );
    Herk( UPPER, ADJOINT, Real(-1), Q, Real(1), Z );
    const Real infOrthogError = HermitianInfinityNorm( UPPER, Z );
    const Real relOrthogError = infOrthogError / (eps*maxDim);
    OutputFromRoot
    (g.Comm(),
     "||Q^H Q - I||_oo / (eps Max(m,n)) = ",relOrthogError);
    PopIndent();

    // Form A - Q R
    OutputFromRoot(g.Comm(),"Testing if A ~= QR...");
    PushIndent();
    LocalGemm( NORMAL, NORMAL, F(-1), Q, R, F(1), A );
    const Real infError = InfinityNorm( A );
    const Real relError = infError / (eps*maxDim*oneNormA);
    OutputFromRoot
    (g.Comm(),"||A - QR||_oo / (eps Max(m,n) ||A||_1) = ",relError);

    PopIndent();

    // TODO: More rigorous failure conditions
    if( relOrthogError > Real(10) )
        LogicError("Unacceptably large relative orthogonality error");
    if( relError > Real(10) )
        LogicError("Unacceptably large relative error");
}

template<typename F>
void TestQR
( const Grid& g,
  Int m,
  Int n,
  bool correctness,
  bool print )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();

    DistMatrix<F,VC,STAR> A(g), AFact(g);
    DistMatrix<F,STAR,STAR> R(g);

    Uniform( A, m, n );
    if( print )
        Print( A, "A" );
    AFact = A;

    Timer timer;

    OutputFromRoot(g.Comm(),"Starting TSQR factorization...");
    mpi::Barrier( g.Comm() );
    timer.Start();
    qr::ExplicitTS( AFact, R );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    const double mD = double(m);
    const double nD = double(n);
    const double gFlops = (2.*mD*nD*nD + 1./3.*nD*nD*nD)/(1.e9*runTime);
    OutputFromRoot(g.Comm(),"Time = ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( AFact, "Q" );
        Print( R, "R" );
    }
    if( correctness )
        TestCorrectness( AFact, R, A );
    PopIndent();
    OutputFromRoot(g.Comm(),"");
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

        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, order );
        SetBlocksize( nb );
        ComplainIfDebug();
        OutputFromRoot(comm,"Will test TSQR");

        TestQR<float>
        ( g, m, n, correctness, print );
        TestQR<Complex<float>>
        ( g, m, n, correctness, print );

        TestQR<double>
        ( g, m, n, correctness, print );
        TestQR<Complex<double>>
        ( g, m, n, correctness, print );

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
