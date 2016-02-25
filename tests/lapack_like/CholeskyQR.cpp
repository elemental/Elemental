/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
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
    const Int n = A.Width();

    // Form I - Q^H Q
    if( g.Rank() == 0 )
        Output("  Testing orthogonality of Q");
    DistMatrix<F> Z(g);
    Identity( Z, n, n );
    DistMatrix<F> Q_MC_MR( Q );
    Herk( UPPER, ADJOINT, Base<F>(-1), Q_MC_MR, Base<F>(1), Z );
    Real oneNormOfError = HermitianOneNorm( UPPER, Z );
    Real infNormOfError = HermitianInfinityNorm( UPPER, Z );
    Real frobNormOfError = HermitianFrobeniusNorm( UPPER, Z );
    if( g.Rank() == 0 )
    {
        Output("    ||Q^H Q - I||_1  = ",oneNormOfError);
        Output("    ||Q^H Q - I||_oo = ",infNormOfError);
        Output("    ||Q^H Q - I||_F  = ",frobNormOfError);
    }

    // Form A - Q R
    if( g.Rank() == 0 )
        Output("  Testing if A = QR");
    const Real oneNormOfA = OneNorm( A );
    const Real infNormOfA = InfinityNorm( A );
    const Real frobNormOfA = FrobeniusNorm( A );
    LocalGemm( NORMAL, NORMAL, F(-1), Q, R, F(1), A );
    oneNormOfError = OneNorm( A );
    infNormOfError = InfinityNorm( A );
    frobNormOfError = FrobeniusNorm( A );
    if( g.Rank() == 0 )
    {
        Output("    ||A||_1       = ",oneNormOfA);
        Output("    ||A||_oo      = ",infNormOfA);
        Output("    ||A||_F       = ",frobNormOfA);
        Output("    ||A - QR||_1  = ",oneNormOfError);
        Output("    ||A - QR||_oo = ",infNormOfError);
        Output("    ||A - QR||_F  = ",frobNormOfError);
    }
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
    DistMatrix<F,VC,STAR> A(g), Q(g);
    DistMatrix<F,STAR,STAR> R(g);

    Uniform( A, m, n );
    if( print )
        Print( A, "A" );
    Q = A;

    if( g.Rank() == 0 )
        Output("  Starting Cholesky QR factorization");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    qr::Cholesky( Q, R );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double mD = double(m);
    const double nD = double(n);
    const double gFlops = (2.*mD*nD*nD + 1./3.*nD*nD*nD)/(1.e9*runTime);
    if( g.Rank() == 0 )
        Output("  Time: ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( Q, "Q" );
        Print( R, "R" );
    }
    if( testCorrectness )
        TestCorrectness( Q, R, A );
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
        mpc::SetPrecision( prec );
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
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
