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
        const Int n = Input("--size","size of HPSD matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> L;
        Uniform( L, n, n );
        MakeTrapezoidal( LOWER, L, -1 );
        DistMatrix<C> A;
        Zeros( A, n, n );
        Herk( LOWER, NORMAL, Real(1), L, Real(0), A );
        if( print )
            Print( A, "A" );

        Timer timer;
        if( mpi::Rank() == 0 )
            timer.Start();
        // Replace A with its matrix square root
        HPSDSquareRoot( LOWER, A );
        if( mpi::Rank() == 0 )
            timer.Stop();
        if( print )
        {
            MakeHermitian( LOWER, A );
            Print( A, "sqrt(A)" );
        }
        if( mpi::Rank() == 0 )
            Output("HPSDSquareRoot time: ",timer.Total()," secs");
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
