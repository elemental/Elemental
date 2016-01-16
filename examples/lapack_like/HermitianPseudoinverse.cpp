/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace std;
using namespace El;

// Typedef our real and complex types to 'R' and 'C' for convenience
typedef double R;
typedef Complex<R> C;

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

        DistMatrix<C> H( n, n );

        // Fill the matrix since we did not pass in a buffer. 
        //
        // We will fill entry (i,j) with the complex value (i+j,i-j) so that 
        // the global matrix is Hermitian. However, only one triangle of the 
        // matrix actually needs to be filled, the symmetry can be implicit.
        //
        const Int localHeight = H.LocalHeight();
        const Int localWidth = H.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = H.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = H.GlobalRow(iLoc);
                H.SetLocal( iLoc, jLoc, C(i+j,i-j) );
            }
        }
        if( print )
            Print( H, "H" );

        Timer timer;
        if( mpi::Rank() == 0 )
            timer.Start();
        // Replace H with its pseudoinverse
        HermitianPseudoinverse( LOWER, H );
        MakeHermitian( LOWER, H );
        if( mpi::Rank() == 0 )
            timer.Stop();

        if( print )
            Print( H, "pinv(H)" );
        if( mpi::Rank() == 0 )
            Output("HermitianPseudoinverse time: ",timer.Total()," secs");
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
