/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_MAKETRAPEZOIDAL_INC
#include ELEM_MAKETRIANGULAR_INC
#include ELEM_BIDIAG_INC
#include ELEM_INFINITYNORM_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_IDENTITY_INC
#include ELEM_UNIFORM_INC
using namespace std;
using namespace elem;

template<typename F> 
void TestCorrectness
( const Matrix<F>& A, const Matrix<F>& tP, const Matrix<F>& tQ,
        Matrix<F>& AOrig,
  bool print, bool display )
{
    typedef Base<F> Real;
    const Int m = AOrig.Height();
    const Int n = AOrig.Width();
    const Real infNormAOrig = InfinityNorm( AOrig );
    const Real frobNormAOrig = FrobeniusNorm( AOrig );
    if( mpi::WorldRank() == 0 )
        cout << "Testing error..." << endl;

    // Grab the diagonal and superdiagonal of the bidiagonal matrix
    auto d = A.GetDiagonal( 0 );
    auto e = A.GetDiagonal( (m>=n ? 1 : -1) );

    // Zero B and then fill its bidiagonal
    Matrix<F> B;
    Zeros( B, m, n );
    B.SetDiagonal( d, 0  );
    B.SetDiagonal( e, (m>=n ? 1 : -1) );
    if( print && mpi::WorldRank() == 0 )
        Print( B, "Bidiagonal" );
    if( display && mpi::WorldRank() == 0 )
        Display( B, "Bidiagonal" );

    if( print || display )
    {
        Matrix<F> Q, P;
        Identity( Q, m, m );
        Identity( P, n, n );
        bidiag::ApplyQ( LEFT,  NORMAL, A, tQ, Q );
        bidiag::ApplyP( RIGHT, NORMAL, A, tP, P );
        if( print && mpi::WorldRank() == 0 )
        {
            Print( Q, "Q" );
            Print( P, "P" );
        }
        if( display && mpi::WorldRank() == 0 )
        {
            Display( Q, "Q" );
            Display( P, "P" );
        }
    }

    // Reverse the accumulated Householder transforms
    bidiag::ApplyQ( LEFT,  ADJOINT, A, tQ, AOrig );
    bidiag::ApplyP( RIGHT, NORMAL,  A, tP, AOrig );
    if( print && mpi::WorldRank() == 0 )
        Print( AOrig, "Manual bidiagonal" );
    if( display && mpi::WorldRank() == 0 )
        Display( AOrig, "Manual bidiagonal" );

    // Compare the appropriate portion of AOrig and B
    if( m >= n )
    {
        MakeTriangular( UPPER, AOrig );
        MakeTrapezoidal( LOWER, AOrig, 1 );
    }
    else
    {
        MakeTriangular( LOWER, AOrig );
        MakeTrapezoidal( UPPER, AOrig, -1 );
    }
    Axpy( F(-1), AOrig, B );
    if( print && mpi::WorldRank() == 0 )
        Print( B, "Error in rotated bidiagonal" );
    if( display && mpi::WorldRank() == 0 )
        Display( B, "Error in rotated bidiagonal" );
    const Real infNormError = InfinityNorm( B );
    const Real frobNormError = FrobeniusNorm( B );

    if( mpi::WorldRank() == 0 )
    {
        cout << "    ||A||_oo = " << infNormAOrig << "\n"
             << "    ||A||_F  = " << frobNormAOrig << "\n"
             << "    ||B - Q^H A P||_oo = " << infNormError << "\n"
             << "    ||B - Q^H A P||_F  = " << frobNormError << endl;
    }
}

template<typename F>
void TestBidiag( Int m, Int n, bool testCorrectness, bool print, bool display )
{
    Matrix<F> A, AOrig;
    Matrix<F> tP, tQ;

    Uniform( A, m, n );
    if( testCorrectness )
    {
        if( mpi::WorldRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( mpi::WorldRank() == 0 )
            cout << "DONE" << endl;
    }
    if( print && mpi::WorldRank() == 0 )
        Print( A, "A" );
    if( display && mpi::WorldRank() == 0 )
        Display( A, "A" );

    if( mpi::WorldRank() == 0 )
    {
        cout << "  Starting bidiagonalization...";
        cout.flush();
    }
    mpi::Barrier( mpi::COMM_WORLD );
    const double startTime = mpi::Time();
    Bidiag( A, tP, tQ );
    mpi::Barrier( mpi::COMM_WORLD );
    const double runTime = mpi::Time() - startTime;
    // TODO: Flop calculation
    if( mpi::WorldRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds." << std::endl;
    }
    if( print && mpi::WorldRank() == 0 )
    {
        Print( A, "A after Bidiag" );
        Print( tP, "tP after Bidiag" );
        Print( tQ, "tQ after Bidiag" );
    }
    if( display && mpi::WorldRank() == 0 )
    {
        Display( A, "A after Bidiag" );
        Display( tP, "tP after Bidiag" );
        Display( tQ, "tQ after Bidiag" );
    }
    if( testCorrectness )
        TestCorrectness( A, tP, tQ, AOrig, print, display );
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );

    try
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool display = Input("--display","display matrices?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );
        ComplainIfDebug();

        if( commRank == 0 )
            cout << "Double-precision:" << endl;
        TestBidiag<double>( m, n, testCorrectness, print, display );

        if( commRank == 0 )
            cout << "Double-precision complex:" << endl;
        TestBidiag<Complex<double>>( m, n, testCorrectness, print, display );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
