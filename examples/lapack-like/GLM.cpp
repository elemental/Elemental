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
#include ELEM_GLM_INC
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

        DistMatrix<F> A(g), B(g), D(g), Y(g);
        Uniform( A, m, n );
        Uniform( B, m, p );
        Uniform( D, m, numRhs );
        auto ACopy( A );
        auto BCopy( B );
        auto DCopy( D );

        if( print )
        {
            Print( A, "A" );
            Print( B, "B" );
            Print( D, "D" );
        }

        GLM( A, B, D, Y );

        if( print ) 
        {
            Print( D, "X" );
            Print( Y, "Y" );
        }
        
        const Real DFrob = FrobeniusNorm( DCopy ); 
        Gemm( NORMAL, NORMAL, F(-1), ACopy, D, F(1), DCopy );
        Gemm( NORMAL, NORMAL, F(-1), BCopy, Y, F(1), DCopy );
        const Real EFrob = FrobeniusNorm( DCopy );
        if( print )
            Print( DCopy, "D - A X - B Y" );
        if( commRank == 0 )
            std::cout << "|| D             ||_F = " << DFrob << "\n"
                      << "|| A X + B Y - D ||_F = " << EFrob << "\n"
                      << std::endl;
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
