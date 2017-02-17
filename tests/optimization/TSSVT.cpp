/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename F>
void TestCorrectness
( const DistMatrix<F,VC,STAR>& A,
        DistMatrix<F,VC,STAR>& B,
  Base<F> tau,
  bool print )
{
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Real eps = limits::Epsilon<Real>();
    const Grid& grid = A.Grid();

    DistMatrix<F> BNormal( A );
    SVT( BNormal, tau );
    if( print )
        Print( BNormal, "BNormal" );
    BNormal -= B;
    if( print )
        Print( BNormal, "E" );
    const Real AFrob = FrobeniusNorm( A );
    const Real BFrob = FrobeniusNorm( B );
    const Real errorFrob = FrobeniusNorm( BNormal );
    if( grid.Rank() == 0 )
        Output
        ("  || A ||_F = ",AFrob,"\n",
         "  || B ||_F = ",BFrob,"\n",
         "  || E ||_F = ",errorFrob);

    // TODO(poulson): A more careful failure condition
    if( errorFrob/BFrob > 100*eps*Max(m,n) )
        LogicError("The error between the two approaches was too high");
}

template<typename F>
void TestSVT
( bool testCorrectness,
  bool print,
  Int m,
  Int n,
  const Grid& grid,
  Base<F> tau )
{
    DistMatrix<F,VC,STAR> A(grid), B(grid);
    Uniform( A, m, n );
    if( print )
        Print( A, "A" );
    B = A;

    if( grid.Rank() == 0 )
        Output("  Starting TS-SVT factorization...");
    mpi::Barrier( grid.Comm() );
    const double startTime = mpi::Time();
    SVT( B, tau );
    mpi::Barrier( grid.Comm() );
    const double runTime = mpi::Time() - startTime;
    if( grid.Rank() == 0 )
        Output("  ",runTime," seconds");
    if( print )
        Print( B, "B" );
    if( testCorrectness )
        TestCorrectness( A, B, tau, print );
}

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
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
            Output("Testing with doubles:");
        TestSVT<double>( testCorrectness, print, m, n, g, tau );

        if( commRank == 0 )
            Output("Testing with double-precision complex:");
        TestSVT<double>( testCorrectness, print, m, n, g, tau );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
