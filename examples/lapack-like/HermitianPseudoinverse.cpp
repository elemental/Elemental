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
#include "elemental/lapack-like/Pseudoinverse.hpp"
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
        const Int n = Input("--size","size of matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> H( n, n );

        // Fill the matrix since we did not pass in a buffer. 
        //
        // We will fill entry (i,j) with the complex value (i+j,i-j) so that 
        // the global matrix is Hermitian. However, only one triangle of the 
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
                H.SetLocal( iLocal, jLocal, C(i+j,i-j) );
            }
        }

        if( print )
            Print( H, "H" );

        // Replace H with its pseudoinverse
        HermitianPseudoinverse( LOWER, H );
        MakeHermitian( LOWER, H );

        if( print )
            Print( H, "pinv(H)" );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
