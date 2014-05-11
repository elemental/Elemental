/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_HEMM_INC
#include ELEM_TRMM_INC
#include ELEM_CHOLESKY_INC
#include ELEM_HERMITIANUNIFORMSPECTRUM_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_INFINITYNORM_INC
#include ELEM_ONENORM_INC
using namespace std;
using namespace elem;

template<typename F>
void TestCorrectness
( UpperOrLower uplo, const Matrix<F>& T, Base<F> alpha, const Matrix<F>& V,
  const Matrix<F>& A )
{
    typedef Base<F> Real;
    const Int m = V.Height();

    Matrix<F> B( A );
    Herk( uplo, NORMAL, F(alpha), V, F(1), B );

    // Test correctness by multiplying a random set of vectors by 
    // A + alpha V V^H, then using the Cholesky factorization to solve.
    Matrix<F> X, Y;
    Uniform( X, m, 100 );
    Zeros( Y, m, 100 );
    Hemm( LEFT, uplo, F(1), B, X, F(0), Y );
    const Real maxNormT = MaxNorm( T );
    const Real maxNormB = HermitianMaxNorm( uplo, B );
    const Real infNormB = HermitianInfinityNorm( uplo, B );
    const Real frobNormB = HermitianFrobeniusNorm( uplo, B );
    const Real oneNormY = OneNorm( Y );
    const Real infNormY = InfinityNorm( Y );
    const Real frobNormY = FrobeniusNorm( Y );

    cholesky::SolveAfter( uplo, NORMAL, T, Y );
    Axpy( F(-1), Y, X );
    const Real oneNormE = OneNorm( X );
    const Real infNormE = InfinityNorm( X );
    const Real frobNormE = FrobeniusNorm( X );

    if( mpi::WorldRank() == 0 )
    {
        cout << "||T||_max = " << maxNormT << "\n"
             << "||B||_max = " << maxNormB << "\n"
             << "||B||_1   = " << infNormB << "\n"
             << "||B||_F   = " << frobNormB << "\n"
             << "||Y||_1   = " << oneNormY << "\n"
             << "||Y||_oo  = " << infNormY << "\n"
             << "||Y||_F   = " << frobNormY << "\n"
             << "||X - inv(B) X||_1  = " << oneNormE << "\n"
             << "||X - inv(B) X||_oo = " << infNormE << "\n"
             << "||X - inv(B) X||_F  = " << frobNormE << endl;
    }
}

template<typename F> 
void TestCholeskyMod
( bool testCorrectness, bool print, UpperOrLower uplo, Int m, Int n, 
  Base<F> alpha )
{
    Matrix<F> T, A;
    HermitianUniformSpectrum( T, m, 1e-9, 10 );
    if( testCorrectness )
    {
        if( mpi::WorldRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        A = T;
        if( mpi::WorldRank() == 0 )
            cout << "DONE" << endl;
    }
    if( print && mpi::WorldRank() == 0 )
        Print( T, "A" );

    if( mpi::WorldRank() == 0 )
    {
        cout << "  Starting Cholesky factorization...";
        cout.flush();
    }
    double startTime = mpi::Time();
    Cholesky( uplo, T );
    double runTime = mpi::Time() - startTime;
    double realGFlops = 1./3.*Pow(double(m),3.)/(1.e9*runTime);
    double gFlops = ( IsComplex<F>::val ? 4*realGFlops : realGFlops );
    if( mpi::WorldRank() == 0 )
    {
        cout << "DONE\n"
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    MakeTriangular( uplo, T );
    if( print && mpi::WorldRank() == 0 )
        Print( T, "Cholesky factor" );

    if( mpi::WorldRank() == 0 )
    {
        cout << "  Generating random update vectors...";
        cout.flush();
    }
    Matrix<F> V, VMod;
    Uniform( V, m, n );
    Scale( F(1)/Sqrt(F(m)*F(n)), V );
    VMod = V;
    if( mpi::WorldRank() == 0 )
        cout << "DONE" << endl;
    if( print && mpi::WorldRank() == 0 )
        Print( V, "V" );

    if( mpi::WorldRank() == 0 )
    {
        cout << "  Starting Cholesky modification...";
        cout.flush();
    }
    startTime = mpi::Time();
    CholeskyMod( uplo, T, alpha, VMod );
    runTime = mpi::Time() - startTime;
    if( mpi::WorldRank() == 0 )
        cout << "DONE\n"
             << "  Time = " << runTime << " seconds." << endl;
    if( print && mpi::WorldRank() == 0 )
        Print( T, "Modified Cholesky factor" );

    if( testCorrectness )
        TestCorrectness( uplo, T, alpha, V, A );
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );

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
        ProcessInput();
        PrintInputReport();

        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test CholeskyMod" << uploChar << endl;

        if( commRank == 0 )
            cout << "Testing with doubles:" << endl;
        TestCholeskyMod<double>
        ( testCorrectness, print, uplo, m, n, alpha );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestCholeskyMod<Complex<double>>
        ( testCorrectness, print, uplo, m, n, alpha );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
