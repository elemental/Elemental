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
    mpi::Comm comm;
    const Int commSize = mpi::Size( comm );

    try 
    {
        const Int m = Input("--m","height of A",100);
        const Int n = Input("--n","width of A",100);
        const Int p = Input("--p","width of B",20);
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

        DistMatrix<F> A(g), B(g), D(g), X(g), Y(g);
        Uniform( A, m, n );
        Uniform( B, m, p );
        Uniform( D, m, numRhs );
        if( print )
        {
            Print( A, "A" );
            Print( B, "B" );
            Print( D, "D" );
        }

        Timer timer;
        if( mpi::Rank() == 0 )
            timer.Start();
        GLM( A, B, D, X, Y );
        if( mpi::Rank() == 0 )
            timer.Stop();
        if( print ) 
        {
            Print( X, "X" );
            Print( Y, "Y" );
        }
        
        const Real DFrob = FrobeniusNorm( D ); 
        Gemm( NORMAL, NORMAL, F(-1), A, D, F(1), D );
        Gemm( NORMAL, NORMAL, F(-1), B, Y, F(1), D );
        const Real EFrob = FrobeniusNorm( D );
        if( print )
            Print( D, "D - A X - B Y" );
        if( mpi::Rank() == 0 )
        {
            Output("GLM time: ",timer.Total()," secs");
            Output
            ("|| D             ||_F = ",DFrob,"\n",
             "|| A X + B Y - D ||_F = ",EFrob,"\n");
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
