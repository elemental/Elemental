/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

typedef double Real;
typedef Complex<Real> F;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commSize = mpi::Size( comm );

    try 
    {
        const Int m = Input("--m","height of A",100);
        const Int n = Input("--n","width of A",100);
        const Int p = Input("--p","number of constraints",20);
        const Int numRhs = Input("--numRhs","# of right-hand sides",5);
        const Int blocksize = Input("--blocksize","algorithmic blocksize",64);
        const bool print = Input("--print","print matrices?",false);
        Int gridHeight = Input("--gridHeight","grid height",0);
        ProcessInput();
        PrintInputReport();

        // Set the algorithmic blocksize
        SetBlocksize( blocksize );

        // If the grid height wasn't specified, then we should attempt to build
        // a nearly-square process grid
        if( gridHeight == 0 )
            gridHeight = Grid::FindFactor( commSize );
        Grid g( comm, gridHeight );

        DistMatrix<F> A(g), B(g), C(g), D(g), X(g);
        Uniform( A, m, n );
        Uniform( B, p, n );
        Uniform( C, m, numRhs );
        Uniform( D, p, numRhs );

        const Real CFrob = FrobeniusNorm( C );
        const Real DFrob = FrobeniusNorm( D );
        if( print )
        {
            Print( A, "A" );
            Print( B, "B" );
            Print( C, "C" );
            Print( D, "D" );
        }

        Timer timer;
        if( mpi::Rank() == 0 )
            timer.Start();
        LSE( A, B, C, D, X );
        if( mpi::Rank() == 0 )
            timer.Stop();
        if( print ) 
            Print( X, "X" );
        
        Gemm( NORMAL, NORMAL, F(-1), A, X, F(1), C );
        Gemm( NORMAL, NORMAL, F(-1), B, X, F(1), D );
        const Real resid = FrobeniusNorm( C );
        const Real EFrob = FrobeniusNorm( D );
        if( print )
            Print( D, "D - B X" );
        if( mpi::Rank() == 0 )
        {
            Output("LSE time: ",timer.Total()," secs");
            Output
            ("|| C - A X ||_F / || C ||_F = ",resid/CFrob,"\n",
             "|| D - B X ||_F / || D ||_F = ",EFrob/DFrob,"\n");
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
