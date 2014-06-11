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
#include EL_UNIFORM_INC
using namespace std;
using namespace El;

template<typename F> 
void TestCorrectness
( UpperOrLower uplo, 
  const DistMatrix<F>& A, 
  const DistMatrix<F,STAR,STAR>& t,
        DistMatrix<F>& AOrig,
  bool print, bool display )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = AOrig.Height();
    const Real infNormAOrig = InfinityNorm( AOrig );
    const Real frobNormAOrig = FrobeniusNorm( AOrig );
    if( g.Rank() == 0 )
        cout << "Testing error..." << endl;

    // Set H to the appropriate Hessenberg portion of A
    DistMatrix<F> H( A );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, H, 1 );
    else
        MakeTrapezoidal( UPPER, H, -1 );
    if( print )
        Print( H, "Hessenberg" );
    if( display )
        Display( H, "Bidiagonal" );

    if( print || display )
    {
        DistMatrix<F> Q(g);
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
    if( print )
        Print( AOrig, "Manual Hessenberg" );
    if( display )
        Display( AOrig, "Manual Hessenberg" );

    // Compare the appropriate portion of AOrig and B
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, AOrig, 1 );
    else
        MakeTrapezoidal( UPPER, AOrig, -1 );
    Axpy( F(-1), AOrig, H );
    if( print )
        Print( H, "Error in rotated Hessenberg" );
    if( display )
        Display( H, "Error in rotated Hessenberg" );
    const Real infNormError = InfinityNorm( H );
    const Real frobNormError = FrobeniusNorm( H );

    if( g.Rank() == 0 )
    {
        cout << "    ||A||_oo = " << infNormAOrig << "\n"
             << "    ||A||_F  = " << frobNormAOrig << "\n"
             << "    ||H - Q^H A Q||_oo = " << infNormError << "\n"
             << "    ||H - Q^H A Q||_F  = " << frobNormError << endl;
    }
}

template<typename F>
void TestHessenberg
( UpperOrLower uplo, Int n, const Grid& g, bool testCorrectness, 
  bool print, bool display )
{
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<F,STAR,STAR> t(g);

    Uniform( A, n, n );
    if( testCorrectness )
    {
        if( g.Rank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.Rank() == 0 )
            cout << "DONE" << endl;
    }
    if( print )
        Print( A, "A" );
    if( display )
        Display( A, "A" );

    if( g.Rank() == 0 )
    {
        cout << "  Starting reduction to Hessenberg form...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Hessenberg( uplo, A, t );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    // TODO: Flop calculation
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds." << std::endl;
    }
    if( print )
    {
        Print( A, "A after Hessenberg" );
        Print( t, "t after Hessenberg" );
    }
    if( display )
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
    const Int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int n = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool display = Input("--display","display matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        ComplainIfDebug();

        if( commRank == 0 )
            cout << "Double-precision:" << endl;
        TestHessenberg<double>( uplo, n, g, testCorrectness, print, display );

        if( commRank == 0 )
            cout << "Double-precision complex:" << endl;
        TestHessenberg<Complex<double>>
        ( uplo, n, g, testCorrectness, print, display );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
