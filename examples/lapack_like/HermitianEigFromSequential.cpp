/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;
using namespace std;
 
// Typedef our real and complex types to 'Real' and 'C' for convenience
typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    // This detects whether or not you have already initialized MPI and 
    // does so if necessary. 
    Environment env( argc, argv );

    // Surround the Elemental calls with try/catch statements in order to 
    // safely handle any exceptions that were thrown during execution.
    try 
    {
        const Int n = Input("--size","size of matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        // Create an n x n complex matrix residing on a single process.
        DistMatrix<C,CIRC,CIRC> H( n, n );
        if( mpi::Rank() == 0 )
        {
            // Set entry (i,j) to (i+j,i-j)
            for( Int j=0; j<n; ++j )
                for( Int i=0; i<n; ++i )
                    H.SetLocal( i, j, C(i+j,i-j) );
        }
        if( print )
            Print( H, "H on process 0" );

        // Call the eigensolver.
        DistMatrix<Real,CIRC,CIRC> w;
        DistMatrix<C,CIRC,CIRC> Q;
        // Optional: set blocksizes and algorithmic choices here. See the 
        //           'Tuning' section of the README for details.
        Timer timer;
        if( mpi::Rank() == 0 )
            timer.Start();
        HermitianEig( LOWER, H, w, Q );
        if( mpi::Rank() == 0 )
            timer.Stop();
        if( print )
        {
            Print( w, "Eigenvalues of H" );
            Print( Q, "Eigenvectors of H" );
        }
        if( mpi::Rank() == 0 )
            Output("HermitianEig time: ",timer.Total()," secs");
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
