/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

// Create a typedef for convenience
typedef double Real;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int n = Input("--size","size of matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<Real> H( n, n );

        // We will fill entry (i,j) with the value i+j so that
        // the global matrix is symmetric. However, only one triangle of the 
        // matrix actually needs to be filled, the symmetry can be implicit.
        const Int localHeight = H.LocalHeight();
        const Int localWidth = H.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            // Our process owns the rows colShift:colStride:n,
            //           and the columns rowShift:rowStride:n
            const Int j = H.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = H.GlobalRow(iLoc);
                H.SetLocal( iLoc, jLoc, Real(i+j) );
            }
        }
        if( print )
            Print( H, "H" );

        // Reform the matrix with the exponentials of the original eigenvalues
        auto expFunc = []( Real alpha ) { return Exp(alpha); };
        HermitianFunction( LOWER, H, function<Real(Real)>(expFunc) );
        if( print )
        {
            MakeHermitian( LOWER, H );
            Print( H, "exp(H)" );
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
