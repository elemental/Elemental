/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level1/MakeHermitian.hpp"
#include "elemental/lapack-like/HermitianFunction.hpp"
using namespace std;
using namespace elem;

// Create a typedef for convenience
typedef double R;

// A functor for returning the exponential of a real number
class ExpFunctor {
public:
    R operator()( R alpha ) const { return std::exp(alpha); }
};

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try 
    {
        const Int n = Input("--size","size of matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        Grid g( mpi::COMM_WORLD );
        DistMatrix<R> H( n, n, g );

        // Fill the matrix since we did not pass in a buffer. 
        //
        // We will fill entry (i,j) with the value i+j so that
        // the global matrix is symmetric. However, only one triangle of the 
        // matrix actually needs to be filled, the symmetry can be implicit.
        //
        const Int colShift = H.ColShift(); // first row we own
        const Int rowShift = H.RowShift(); // first col we own
        const Int colStride = H.ColStride();
        const Int rowStride = H.RowStride();
        const Int localHeight = H.LocalHeight();
        const Int localWidth = H.LocalWidth();
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
            {
                // Our process owns the rows colShift:colStride:n,
                //           and the columns rowShift:rowStride:n
                const Int i = colShift + iLocal*colStride;
                const Int j = rowShift + jLocal*rowStride;
                H.SetLocal( iLocal, jLocal, R(i+j) );
            }
        }

        if( print )
            Print( H, "H" );

        // Reform the matrix with the exponentials of the original eigenvalues
        RealHermitianFunction( LOWER, H, ExpFunctor() );

        if( print )
        {
            MakeHermitian( LOWER, H );
            Print( H, "exp(H)" );
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
