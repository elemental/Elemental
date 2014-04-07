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
( const DistMatrix<F>& A, 
  const DistMatrix<F,STAR,STAR>& tP,
  const DistMatrix<F,STAR,STAR>& tQ,
        DistMatrix<F>& AOrig,
  bool print, bool display )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = AOrig.Height();
    const Int n = AOrig.Width();
    const Real infNormAOrig = InfinityNorm( AOrig );
    const Real frobNormAOrig = FrobeniusNorm( AOrig );
    if( g.Rank() == 0 )
        cout << "Testing error..." << endl;

    // Grab the diagonal and superdiagonal of the bidiagonal matrix
    auto d = A.GetDiagonal( 0 );
    auto e = A.GetDiagonal( (m>=n ? 1 : -1) );

    // Zero B and then fill its bidiagonal
    DistMatrix<F> B(g);
    B.AlignWith( A );
    Zeros( B, m, n );
    B.SetDiagonal( d, 0  );
    B.SetDiagonal( e, (m>=n ? 1 : -1) );
    if( print )
        Print( B, "Bidiagonal" );
    if( display )
        Display( B, "Bidiagonal" );

    if( print || display )
    {
        DistMatrix<F> Q(g), P(g);
        Identity( Q, m, m );
        Identity( P, n, n );
        bidiag::ApplyQ( LEFT,  NORMAL, A, tQ, Q );
        bidiag::ApplyP( RIGHT, NORMAL, A, tP, P );
        if( print )
        {
            Print( Q, "Q" );
            Print( P, "P" );
        }
        if( display )
        {
            Display( Q, "Q" );
            Display( P, "P" );
        }
    }

    // Reverse the accumulated Householder transforms
    bidiag::ApplyQ( LEFT,  ADJOINT, A, tQ, AOrig );
    bidiag::ApplyP( RIGHT, NORMAL,  A, tP, AOrig );
    if( print )
        Print( AOrig, "Manual bidiagonal" );
    if( display )
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
    if( print )
        Print( B, "Error in rotated bidiagonal" );
    if( display )
        Display( B, "Error in rotated bidiagonal" );
    const Real infNormError = InfinityNorm( B );
    const Real frobNormError = FrobeniusNorm( B );

    if( g.Rank() == 0 )
    {
        cout << "    ||A||_oo = " << infNormAOrig << "\n"
             << "    ||A||_F  = " << frobNormAOrig << "\n"
             << "    ||B - Q^H A P||_oo = " << infNormError << "\n"
             << "    ||B - Q^H A P||_F  = " << frobNormError << endl;
    }
}

template<typename F>
void TestBidiag
( Int m, Int n, const Grid& g, bool testCorrectness, bool print, bool display )
{
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<F,STAR,STAR> tP(g), tQ(g);

    Uniform( A, m, n );
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
        cout << "  Starting bidiagonalization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Bidiag( A, tP, tQ );
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
        Print( A, "A after Bidiag" );
        Print( tP, "tP after Bidiag" );
        Print( tQ, "tQ after Bidiag" );
    }
    if( display )
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
    const Int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
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
        SetBlocksize( nb );
        ComplainIfDebug();

        if( commRank == 0 )
            cout << "Double-precision:" << endl;
        TestBidiag<double>( m, n, g, testCorrectness, print, display );

        if( commRank == 0 )
            cout << "Double-precision complex:" << endl;
        TestBidiag<Complex<double>>( m, n, g, testCorrectness, print, display );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
