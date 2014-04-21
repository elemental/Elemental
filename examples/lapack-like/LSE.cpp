/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_FROBENIUSNORM_INC
#include ELEM_LSE_INC
#include ELEM_UNIFORM_INC
using namespace elem;

typedef double Real;
typedef Complex<Real> F;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );

    try 
    {
        const Int m = Input("--m","height of A",100);
        const Int n = Input("--n","width of A",100);
        const Int p = Input("--p","number of constraints",20);
        const Int numRhs = Input("--numRhs","# of right-hand sides",5);
        const Int blocksize = Input("--blocksize","algorithmic blocksize",64);
        const bool resid = Input("--resid","compute residual?",true);
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
        auto BCopy( B );
        auto DCopy( D );

        const Real AFrob = FrobeniusNorm( A );
        const Real BFrob = FrobeniusNorm( B );
        const Real CFrob = FrobeniusNorm( C );
        const Real DFrob = FrobeniusNorm( D );
        if( print )
        {
            Print( A, "A" );
            Print( B, "B" );
            Print( C, "C" );
            Print( D, "D" );
        }

        LSE( A, B, C, D, X, resid );

        if( print ) 
        { 
            Print( X, "X" );
            if( resid )
                Print( C, "rotated residuals" );
        }
        
        Gemm( NORMAL, NORMAL, F(-1), BCopy, X, F(1), DCopy );
        const Real EFrob = FrobeniusNorm( DCopy );
        if( print )
            Print( DCopy, "D - B X" );
        if( commRank == 0 )
            std::cout << "|| A       ||_F = " << AFrob << "\n"
                      << "|| B       ||_F = " << BFrob << "\n"
                      << "|| C       ||_F = " << CFrob << "\n"
                      << "|| D       ||_F = " << DFrob << "\n"
                      << "|| B X - D ||_F = " << EFrob << "\n"
                      << std::endl;
        if( resid )
        {
            const Real residFrob = FrobeniusNorm( C );
            if( commRank == 0 )
                std::cout << "|| A X - C ||_F = " << residFrob << std::endl;
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
