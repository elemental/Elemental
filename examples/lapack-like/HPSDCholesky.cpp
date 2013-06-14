/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/HPSDCholesky.hpp"
#include "elemental/matrices/Uniform.hpp"
using namespace std;
using namespace elem;

// Typedef our real and complex types to 'R' and 'C' for convenience
typedef double R;
typedef Complex<R> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try 
    {
        const int n = Input("--size","size of HPSD matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        Grid g( mpi::COMM_WORLD );
        DistMatrix<C> L(g), A(g);
        Uniform( L, n, n );
        MakeTrapezoidal( LOWER, L, -1 );
        Zeros( A, n, n );
        Herk( LOWER, NORMAL, C(1), L, C(0), A );

        if( print )
            Print( A, "A" );

        // Replace A with its Cholesky factor
        HPSDCholesky( LOWER, A );

        if( print )
        {
            MakeTriangular( LOWER, A );
            Print( A, "chol(A)" );
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
