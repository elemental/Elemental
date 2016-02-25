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
( UpperOrLower uplo,
  const Matrix<F>& T,
  Base<F> alpha,
  const Matrix<F>& V,
  const Matrix<F>& A )
{
    typedef Base<F> Real;
    const Int m = V.Height();

    Matrix<F> B( A );
    Herk( uplo, NORMAL, alpha, V, Real(1), B );

    // Test correctness by multiplying a random set of vectors by 
    // A + alpha V V^H, then using the Cholesky factorization to solve.
    Matrix<F> X, Y;
    Uniform( X, m, 100 );
    Zeros( Y, m, 100 );
    Hemm( LEFT, uplo, F(1), B, X, F(0), Y );
    const Real maxNormT = MaxNorm( T );
    const Real maxNormB = HermitianMaxNorm( uplo, B );
    const Real frobNormB = HermitianFrobeniusNorm( uplo, B );
    const Real frobNormY = FrobeniusNorm( Y );

    cholesky::SolveAfter( uplo, NORMAL, T, Y );
    X -= Y;
    const Real frobNormE = FrobeniusNorm( X );

    if( mpi::Rank() == 0 )
        Output
        ("||T||_max = ",maxNormT,"\n",
         "||B||_max = ",maxNormB,"\n",
         "||B||_F   = ",frobNormB,"\n",
         "||Y||_F   = ",frobNormY,"\n",
         "||X - inv(B) X||_F  = ",frobNormE);
}

template<typename F> 
void TestCholeskyMod
( UpperOrLower uplo,
  Int m,
  Int n, 
  Base<F> alpha,
  bool testCorrectness,
  bool print )
{
    if( mpi::Rank() == 0 )
        Output("Testing with ",TypeName<F>());
    Matrix<F> T, A;
    HermitianUniformSpectrum( T, m, 1e-9, 10 );
    if( testCorrectness )
        A = T;
    if( print && mpi::Rank() == 0 )
        Print( T, "A" );

    if( mpi::Rank() == 0 )
        Output("  Starting Cholesky...");
    double startTime = mpi::Time();
    Cholesky( uplo, T );
    double runTime = mpi::Time() - startTime;
    double realGFlops = 1./3.*Pow(double(m),3.)/(1.e9*runTime);
    double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    if( mpi::Rank() == 0 )
        Output("  Time = ",runTime," seconds (",gFlops," GFlop/s)");
    MakeTrapezoidal( uplo, T );
    if( print && mpi::Rank() == 0 )
        Print( T, "Cholesky factor" );

    Matrix<F> V, VMod;
    Uniform( V, m, n );
    V *= F(1)/Sqrt(F(m)*F(n));
    VMod = V;
    if( print && mpi::Rank() == 0 )
        Print( V, "V" );

    if( mpi::Rank() == 0 )
        Output("  Starting Cholesky mod...");
    startTime = mpi::Time();
    CholeskyMod( uplo, T, alpha, VMod );
    runTime = mpi::Time() - startTime;
    if( mpi::Rank() == 0 )
        Output("  Time = ",runTime," seconds");
    if( print && mpi::Rank() == 0 )
        Print( T, "Modified Cholesky factor" );

    if( testCorrectness )
        TestCorrectness( uplo, T, alpha, V, A );
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int m = Input("--m","height of matrix",100);
        const Int n = Input("--n","rank of update",5);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const double alpha = Input("--alpha","update scaling",3.);
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

        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        ComplainIfDebug();

        TestCholeskyMod<float>
        ( uplo, m, n, alpha, testCorrectness, print );
        TestCholeskyMod<Complex<float>>
        ( uplo, m, n, alpha, testCorrectness, print );

        TestCholeskyMod<double>
        ( uplo, m, n, alpha, testCorrectness, print );
        TestCholeskyMod<Complex<double>>
        ( uplo, m, n, alpha, testCorrectness, print );

#ifdef EL_HAVE_QD
        TestCholeskyMod<DoubleDouble>
        ( uplo, m, n, alpha, testCorrectness, print );
        TestCholeskyMod<QuadDouble>
        ( uplo, m, n, alpha, testCorrectness, print );
#endif

#ifdef EL_HAVE_QUAD
        TestCholeskyMod<Quad>
        ( uplo, m, n, alpha, testCorrectness, print );
        TestCholeskyMod<Complex<Quad>>
        ( uplo, m, n, alpha, testCorrectness, print );
#endif

#ifdef EL_HAVE_MPC
        TestCholeskyMod<BigFloat>
        ( uplo, m, n, alpha, testCorrectness, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
