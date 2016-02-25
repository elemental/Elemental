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
( bool pivoted,
  bool print, 
  const Matrix<F>& A,
  const Permutation& P,
  const Matrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Int m = AOrig.Height();

    Output("Testing error...");

    // Generate random right-hand sides
    Matrix<F> X, Y;
    Uniform( X, m, 100 );
    Y = X;
    if( pivoted )
        P.PermuteRows( Y );

    // Solve against the (pivoted) right-hand sides
    Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, Y );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, Y );

    // Now investigate the residual, ||AOrig Y - X||_oo
    const Real infNormX = InfinityNorm( X );
    const Real frobNormX = FrobeniusNorm( X );
    Gemm( NORMAL, NORMAL, F(-1), AOrig, Y, F(1), X );
    const Real infNormError = InfinityNorm( X );
    const Real frobNormError = FrobeniusNorm( X );
    const Real infNormA = InfinityNorm( AOrig );
    const Real frobNormA = FrobeniusNorm( AOrig );

    Output
    ("||A||_oo                 = ",infNormA,"\n",
     "||A||_F                  = ",frobNormA,"\n",
     "||X||_oo                 = ",infNormX,"\n",
     "||X||_F                  = ",frobNormX,"\n",
     "||A U^-1 L^-1 X - X||_oo = ",infNormError,"\n",
     "||A U^-1 L^-1 X - X||_F  = ",frobNormError);
}

template<typename F> 
void TestLU( Int m, bool pivot, bool testCorrectness, bool print )
{
    if( mpi::Rank() == 0 )
        Output("Testing with ",TypeName<F>());
    Matrix<F> A, ARef;
    Permutation P;

    Uniform( A, m, m );
    if( testCorrectness )
        ARef = A;
    if( print )
        Print( A, "A" );

    Output("  Starting LU factorization...");
    const double startTime = mpi::Time();
    if( pivot )
        LU( A, P );
    else
        LU( A );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 2./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    Output("  ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        if( pivot )
        {
            // TODO: Print permutation
            /*
            PivotsToPermutation( rowPiv, p );
            Print( p, "p after factorization" );
            Matrix<Int> P;
            ExplicitPermutation( p, P );
            Print( P, "P after factorization" );
            */
        }
    }
    if( testCorrectness )
        TestCorrectness( pivot, print, A, P, ARef );
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );

    try
    {
        const Int m = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool pivot = Input("--pivot","pivoted LU?",true);
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

        SetBlocksize( nb );
        ComplainIfDebug();
        if( commRank == 0 )
            Output("Will test LU",( pivot ? " with partial pivoting" : " " ));

        TestLU<float>( m, pivot, testCorrectness, print );
        TestLU<Complex<float>>( m, pivot, testCorrectness, print );

        TestLU<double>( m, pivot, testCorrectness, print );
        TestLU<Complex<double>>( m, pivot, testCorrectness, print );

#ifdef EL_HAVE_QD
        TestLU<DoubleDouble>( m, pivot, testCorrectness, print );
        TestLU<QuadDouble>( m, pivot, testCorrectness, print );
#endif

#ifdef EL_HAVE_QUAD
        TestLU<Quad>( m, pivot, testCorrectness, print );
        TestLU<Complex<Quad>>( m, pivot, testCorrectness, print );
#endif

#ifdef EL_HAVE_MPC
        TestLU<BigFloat>( m, pivot, testCorrectness, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
