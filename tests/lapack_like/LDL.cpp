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
( bool conjugated,
  bool print, 
  const DistMatrix<F>& A,
  const DistMatrix<F,MD,STAR>& dSub,
  const DistPermutation& p,
  const DistMatrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = AOrig.Height();

    DistMatrix<F> X(g), Y(g);
    Uniform( X, m, 100 );
    Y = X;

    // Test correctness by comparing the application of AOrig against a 
    // random set of 100 vectors to the application of tril(A) tril(A)^H
    if( print )
        Print( X, "X" );
    ldl::MultiplyAfter( A, dSub, p, Y, conjugated );
    if( print )
        Print( Y, "P' L B L' P X" );
    Symm( LEFT, LOWER, F(-1), AOrig, X, F(1), Y, conjugated );
    if( print )
        Print( Y, "P' L B L' P X - A X" );
    const Real infNormOfError = InfinityNorm( Y );
    const Real frobNormOfError = FrobeniusNorm( Y );
    const Real infNormOfA = HermitianInfinityNorm( LOWER, AOrig );
    const Real frobNormOfA = HermitianFrobeniusNorm( LOWER, AOrig );
    const Real infNormOfX = InfinityNorm( X );
    const Real frobNormOfX = FrobeniusNorm( X );
    if( g.Rank() == 0 )
        Output
        ("||A||_oo   = ",infNormOfA,"\n",
         "||A||_F    = ",frobNormOfA,"\n",
         "||X||_oo   = ",infNormOfX,"\n",
         "||X||_F    = ",frobNormOfX,"\n",
         "||A X - L D L^[T/H] X||_oo = ",infNormOfError,"\n",
         "||A X - L D L^[T/H] X||_F  = ",frobNormOfError);
}

template<typename F> 
void TestLDL
( const Grid& g,
  Int m,
  bool conjugated,
  Int nbLocal,
  bool testCorrectness,
  bool print )
{
    DistMatrix<F> A(g), AOrig(g);
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());

    SetLocalTrrkBlocksize<F>( nbLocal );

    if( conjugated )
        HermitianUniformSpectrum( A, m, -100, 100 );
    else
        Uniform( A, m, m );
    if( testCorrectness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    if( g.Rank() == 0 )
        Output("  Starting LDL^[T/H] factorization...");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    DistMatrix<F,MD,STAR> dSub(g);
    DistPermutation p(g);
    LDL( A, dSub, p, conjugated );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 1./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
        Output("  ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        DistMatrix<Int> P(g);
        p.ExplicitMatrix( P ); 
        Print( P, "P" );
    }
    if( testCorrectness )
        TestCorrectness( conjugated, print, A, dSub, p, AOrig );
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--gridHeight","process grid height",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int m = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool conjugated = Input("--conjugate","conjugate LDL?",false);
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

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        SetBlocksize( nb );
        ComplainIfDebug();

        TestLDL<float>
        ( g, m, conjugated, nbLocal, testCorrectness, print );
        TestLDL<Complex<float>>
        ( g, m, conjugated, nbLocal, testCorrectness, print );

        TestLDL<double>
        ( g, m, conjugated, nbLocal, testCorrectness, print );
        TestLDL<Complex<double>>
        ( g, m, conjugated, nbLocal, testCorrectness, print );

#ifdef EL_HAVE_QD
        TestLDL<DoubleDouble>
        ( g, m, conjugated, nbLocal, testCorrectness, print );
        TestLDL<QuadDouble>
        ( g, m, conjugated, nbLocal, testCorrectness, print );
#endif

#ifdef EL_HAVE_QUAD
        TestLDL<Quad>
        ( g, m, conjugated, nbLocal, testCorrectness, print );
        TestLDL<Complex<Quad>>
        ( g, m, conjugated, nbLocal, testCorrectness, print );
#endif

#ifdef EL_HAVE_MPC
        TestLDL<BigFloat>
        ( g, m, conjugated, nbLocal, testCorrectness, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
