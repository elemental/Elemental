/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template<typename F> 
void TestCorrectness
( const DistMatrix<F,VC,  STAR>& Q,
  const DistMatrix<F,STAR,STAR>& R,
        DistMatrix<F,VC,  STAR>& A )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = A.Width();

    // Form I - Q^H Q
    if( g.Rank() == 0 )
        Output("  Testing orthogonality of Q...");
    DistMatrix<F> Z(g);
    Identity( Z, n, n );
    Herk( UPPER, ADJOINT, Real(-1), Q, Real(1), Z );
    Real oneNormError = HermitianOneNorm( UPPER, Z );
    Real infNormError = HermitianInfinityNorm( UPPER, Z );
    Real frobNormError = HermitianFrobeniusNorm( UPPER, Z );
    if( g.Rank() == 0 )
        Output
        ("    ||Q^H Q - I||_1  = ",oneNormError,"\n",
         "    ||Q^H Q - I||_oo = ",infNormError,"\n",
         "    ||Q^H Q - I||_F  = ",frobNormError);

    // Form A - Q R
    if( g.Rank() == 0 )
        Output("  Testing if A = QR...");
    const Real oneNormA = OneNorm( A );
    const Real infNormA = InfinityNorm( A );
    const Real frobNormA = FrobeniusNorm( A );
    LocalGemm( NORMAL, NORMAL, F(-1), Q, R, F(1), A );
    oneNormError = OneNorm( A );
    infNormError = InfinityNorm( A );
    frobNormError = FrobeniusNorm( A );
    if( g.Rank() == 0 )
        Output
        ("    ||A||_1       = ",oneNormA,"\n",
         "    ||A||_oo      = ",infNormA,"\n",
         "    ||A||_F       = ",frobNormA,"\n",
         "    ||A - QR||_1  = ",oneNormError,"\n",
         "    ||A - QR||_oo = ",infNormError,"\n",
         "    ||A - QR||_F  = ",frobNormError);
}

template<typename F>
void TestQR
( const Grid& g,
  Int m,
  Int n,
  bool testCorrectness,
  bool print )
{
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());
    DistMatrix<F,VC,STAR> A(g), AFact(g);
    DistMatrix<F,STAR,STAR> R(g);

    Uniform( A, m, n );
    if( print )
        Print( A, "A" );
    AFact = A;

    if( g.Rank() == 0 )
        Output("  Starting TSQR factorization...");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    qr::ExplicitTS( AFact, R );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double mD = double(m);
    const double nD = double(n);
    const double gFlops = (2.*mD*nD*nD + 1./3.*nD*nD*nD)/(1.e9*runTime);
    if( g.Rank() == 0 )
        Output("  Time = ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( AFact, "Q" );
        Print( R, "R" );
    }
    if( testCorrectness )
        TestCorrectness( AFact, R, A );
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );

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
        mpc::SetPrecision( prec );
#endif

        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, order );
        SetBlocksize( nb );
        ComplainIfDebug();
        if( commRank == 0 )
            Output("Will test TSQR");

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
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
