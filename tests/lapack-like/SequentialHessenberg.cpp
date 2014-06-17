/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"
#include EL_IDENTITY_INC

using namespace std;
using namespace El;

template<typename F> 
void TestCorrectness
( UpperOrLower uplo, const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& AOrig,
  bool print, bool display )
{
    typedef Base<F> Real;
    const Int n = AOrig.Height();
    const Real infNormAOrig = InfinityNorm( AOrig );
    const Real frobNormAOrig = FrobeniusNorm( AOrig );
    if( mpi::WorldRank() == 0 )
        cout << "Testing error..." << endl;

    // Set H to the appropriate Hessenberg portion of A
    Matrix<F> H( A );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, H, 1 );
    else
        MakeTrapezoidal( UPPER, H, -1 );
    if( print && mpi::WorldRank() == 0 )
        Print( H, "Hessenberg" );
    if( display && mpi::WorldRank() == 0 )
        Display( H, "Bidiagonal" );

    if( (print || display) && mpi::WorldRank() == 0 )
    {
        Matrix<F> Q;
        Identity( Q, n, n );
        hessenberg::ApplyQ( uplo, LEFT, NORMAL, A, t, Q );
        if( print )
            Print( Q, "Q" );
        if( display )
            Display( Q, "Q" );
    }

    // Reverse the accumulated Householder transforms
    hessenberg::ApplyQ( uplo, LEFT, ADJOINT, A, t, AOrig );
    hessenberg::ApplyQ( uplo, RIGHT, NORMAL, A, t, AOrig );
    if( print && mpi::WorldRank() == 0 )
        Print( AOrig, "Manual Hessenberg" );
    if( display && mpi::WorldRank() == 0 )
        Display( AOrig, "Manual Hessenberg" );

    // Compare the appropriate portion of AOrig and B
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, AOrig, 1 );
    else
        MakeTrapezoidal( UPPER, AOrig, -1 );
    Axpy( F(-1), AOrig, H );
    if( print && mpi::WorldRank() == 0 )
        Print( H, "Error in rotated Hessenberg" );
    if( display && mpi::WorldRank() == 0 )
        Display( H, "Error in rotated Hessenberg" );
    const Real infNormError = InfinityNorm( H );
    const Real frobNormError = FrobeniusNorm( H );

    if( mpi::WorldRank() == 0 )
    {
        cout << "    ||A||_oo = " << infNormAOrig << "\n"
             << "    ||A||_F  = " << frobNormAOrig << "\n"
             << "    ||H - Q^H A Q||_oo = " << infNormError << "\n"
             << "    ||H - Q^H A Q||_F  = " << frobNormError << endl;
    }
}

template<typename F>
void TestHessenberg
( UpperOrLower uplo, Int n, bool testCorrectness, bool print, bool display )
{
    Matrix<F> A, AOrig;
    Matrix<F> t;

    Uniform( A, n, n );
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
        cout << "  Starting reduction to Hessenberg form...";
        cout.flush();
    }
    mpi::Barrier( mpi::COMM_WORLD );
    const double startTime = mpi::Time();
    Hessenberg( uplo, A, t );
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
        Print( A, "A after Hessenberg" );
        Print( t, "t after Hessenberg" );
    }
    if( display && mpi::WorldRank() == 0 )
    {
        Display( A, "A after Hessenberg" );
        Display( t, "t after Hessenberg" );
    }
    if( testCorrectness )
        TestCorrectness( uplo, A, t, AOrig, print, display );
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
        const Int n = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool display = Input("--display","display matrices?",false);
        ProcessInput();
        PrintInputReport();

        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        ComplainIfDebug();

        if( commRank == 0 )
            cout << "Double-precision:" << endl;
        TestHessenberg<double>( uplo, n, testCorrectness, print, display );

        if( commRank == 0 )
            cout << "Double-precision complex:" << endl;
        TestHessenberg<Complex<double>>
        ( uplo, n, testCorrectness, print, display );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
