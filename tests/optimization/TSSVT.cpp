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
( const DistMatrix<F,VC,  STAR>& A,
        DistMatrix<F,VC,  STAR>& B, Base<F> tau )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();

    DistMatrix<F> BNormal( A );
    SVT( BNormal, tau );
    Axpy( F(-1), B, BNormal );
    const Real AFrob = FrobeniusNorm( A );
    const Real BFrob = FrobeniusNorm( B );
    const Real errorFrob = FrobeniusNorm( BNormal );
    if( g.Rank() == 0 )
        cout << "  || A ||_F = " << AFrob << "\n"
             << "  || B ||_F = " << BFrob << "\n"
             << "  || E ||_F = " << errorFrob << "\n"
             << endl;
}

template<typename F>
void TestSVT
( bool testCorrectness, bool print,
  Int m, Int n, const Grid& g, Base<F> tau )
{
    DistMatrix<F,VC,STAR> A(g), B(g);
    Uniform( A, m, n );
    if( print )
        Print( A, "A" );
    B = A;

    if( g.Rank() == 0 )
    {
        cout << "  Starting TS-SVT factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    SVT( B, tau );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    if( g.Rank() == 0 )
    {
        cout << "DONE. Time = " << runTime << " seconds." << endl;
    }
    if( print )
        Print( B, "B" );
    if( testCorrectness )
        TestCorrectness( A, B, tau );
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );

    try
    {
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const double tau = Input("--tau","soft-threshold parameter",0.5);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, order );
        SetBlocksize( nb );
        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test TSQR" << endl;

        if( commRank == 0 )
            cout << "Testing with doubles:" << endl;
        TestSVT<double>( testCorrectness, print, m, n, g, tau );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestSVT<double>( testCorrectness, print, m, n, g, tau );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
