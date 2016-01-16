/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

// Typedef our real and complex types to 'Real' and 'C' for convenience
typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> A;
        Uniform( A, m, n );

        Timer timer;
        // Compute the pseudoinverseof A (but do not overwrite A)
        DistMatrix<C> pinvA( A );
        if( mpi::Rank() == 0 )
            timer.Start();
        Pseudoinverse( pinvA );
        if( mpi::Rank() == 0 )
            timer.Stop();
        if( print )
        {
            Print( A, "A" );
            Print( pinvA, "pinv(A)" );
        }

        const Real frobA = FrobeniusNorm( A );
        const Real frobPinvA = FrobeniusNorm( pinvA );

        if( mpi::Rank() == 0 )
        {
            Output("PseudoInverse time: ",timer.Total()," secs");
            Output
            ("||   A     ||_F = ",frobA,"\n",
             "|| pinv(A) ||_F = ",frobPinvA,"\n");
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
