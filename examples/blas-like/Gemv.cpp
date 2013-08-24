/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/matrices/Uniform.hpp"
using namespace std;
using namespace elem;

// Typedef our real and complex types to 'Real' and 'C' for convenience
typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try 
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const bool adjoint = Input("--adjoint","apply adjoint?",false);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        const Orientation orientation = ( adjoint ? ADJOINT : NORMAL );

        DistMatrix<C> A;
        Uniform( A, m, n );

        // Draw the entries of the original x and y from uniform distributions 
        // over the complex unit ball
        DistMatrix<C,VC,STAR> x, y;
        if( orientation == NORMAL )
        {
            Uniform( x, n, 1 );
            Uniform( y, m, 1 );
        }
        else
        {
            Uniform( x, m, 1 );
            Uniform( y, n, 1 );
        }

        if( print )
        {
            Print( A, "A" );
            Print( x, "x" );
            Print( y, "y" );
        }

        // Run the matrix-vector product
        Gemv( orientation, C(3), A, x, C(4), y );

        if( print )
        {
            if( orientation == NORMAL )
                Print( y, "y := 3 A x + 4 y" );
            else
                Print( y, "y := 3 A^H x + 4 y" );
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
